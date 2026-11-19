//! Geology layer: geologic provinces as a partition, plus the geology timeline.
//!
//! A port of geoforge v1 Stage 2 (same 20 province types, same painting order, same
//! constants), on geodesic distances instead of pixel steps. Changes from v1, each because
//! v1's shortcut broke on the sphere or on authored data:
//!
//! - Crust is tested per cell, so hand-painted continents work; a plate's "kind" is its
//!   majority crust.
//! - Widths are geodesic km (v1's latitude factor widened features toward the poles).
//! - Shields and basins are spherical caps around the continental centroid, not pixel-space
//!   ellipses around a bounding-box centre (which broke across the antimeridian).
//! - Paleo-orogens follow great circles; hotspot tracks follow the plate's Euler rotation.
//! - Large igneous provinces are contiguous blobs (v1 never placed them: its "unassigned
//!   pixel" filter excluded every pixel, since the base layers cover the whole map).
//! - In ocean–ocean subduction the older crust subducts (v1 used the lower plate id).

pub mod lithology;
pub mod provinces;

use crate::tectonics::boundaries::{self, BoundaryParams};
use geoforge_core::grid::{distance_field, Cell, Grid};
use geoforge_core::key::Key;
use geoforge_core::product::{Generator, LayerId, Precondition, ProductData, ProductId, Products, RunContext};
use geoforge_core::shapes::{Event, EventKind, Field, Partition, Place, Table, Timeline};
use geoforge_core::surface::{BoundaryKind, BoundarySegment, CrustType, GeologicProvince, Plate, Province};
use geoforge_core::{BodyContext, Error, UnitVec3};
use provinces as pc;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::collections::{BTreeMap, BinaryHeap, HashMap};

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct GeologyParams {
    pub lip_probability: f64,
    pub boundaries: BoundaryParams,
    pub oceanic: bool,
    pub stable_regions: bool,
    pub arc_systems: bool,
    pub extensional: bool,
}

impl Default for GeologyParams {
    fn default() -> Self {
        Self {
            lip_probability: pc::LIP_PROBABILITY,
            boundaries: BoundaryParams::default(),
            oceanic: true,
            stable_regions: true,
            arc_systems: true,
            extensional: true,
        }
    }
}

pub struct Geology;

impl Generator for Geology {
    type Params = GeologyParams;
    const NAME: &'static str = "geology";
    const VERSION: u32 = 1;
    const LAYER: LayerId = LayerId::Geology;

    fn needs() -> &'static [ProductId] {
        &[ProductId::Plates, ProductId::CrustType]
    }
    fn wants() -> &'static [ProductId] {
        &[ProductId::CrustAge]
    }
    fn applies(body: &BodyContext) -> Result<(), Precondition> {
        crate::preconditions::plate_tectonics(body)
    }

    fn run(ctx: &RunContext<'_>, inputs: &Products, p: &GeologyParams) -> Result<Products, Error> {
        let plates = inputs.plates().ok_or(Error::MissingInput(ProductId::Plates))?;
        let crust = inputs.crust_type().ok_or(Error::MissingInput(ProductId::CrustType))?;
        let mut out = Products::new();
        let provinces = if inputs.provinces().is_some() && inputs.geology_events().is_some() {
            None
        } else {
            let segments = boundaries::analyze(ctx.grid, plates, ctx.radius_km(), &p.boundaries);
            let (provinces, timeline) = build(
                ctx.grid,
                plates,
                &crust.values,
                inputs.crust_age().map(|f| f.values.as_slice()),
                &segments,
                ctx.radius_km(),
                ctx.seed(ProductId::Provinces).value(),
                p,
            );
            if inputs.geology_events().is_none() {
                out.insert(ProductData::GeologyEvents(timeline));
            }
            Some(provinces)
        };
        let provinces = match (inputs.provinces(), provinces) {
            (Some(pv), _) => pv.clone(),
            (None, Some(pv)) => {
                out.insert(ProductData::Provinces(pv.clone()));
                pv
            }
            (None, None) => unreachable!("built above when provinces are missing"),
        };
        if inputs.lithology().is_none() || inputs.rock_age().is_none() {
            let (lith, age) = lithology::build(
                ctx.grid,
                &provinces,
                &crust.values,
                inputs.crust_age().map(|f| f.values.as_slice()),
                ctx.seed(ProductId::Lithology).value(),
            );
            if inputs.lithology().is_none() {
                out.insert(ProductData::Lithology(Field::from_values(ctx.grid.spec(), lith)));
            }
            if inputs.rock_age().is_none() {
                out.insert(ProductData::RockAge(Field::from_values(ctx.grid.spec(), age)));
            }
        }
        Ok(out)
    }
}

#[derive(PartialEq)]
struct Item(f64, Cell);
impl Eq for Item {}
impl Ord for Item {
    fn cmp(&self, o: &Self) -> Ordering {
        o.0.partial_cmp(&self.0).unwrap_or(Ordering::Equal).then_with(|| o.1.cmp(&self.1))
    }
}
impl PartialOrd for Item {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}

struct PlateInfo {
    area_km2: f64,
    cont: Vec<Cell>,
    ocean: Vec<Cell>,
    cont_area_km2: f64,
    kind: CrustType,
}

#[derive(Clone, Copy)]
struct Side {
    plate: u32,
    continental: bool,
    mean_age: f64,
}

struct B<'a> {
    grid: &'a dyn Grid,
    owner: &'a [u32],
    crust: &'a [CrustType],
    age: Option<&'a [f32]>,
    r: f64,
    plates: &'a Partition<Plate>,
    info: BTreeMap<u32, PlateInfo>,
    edge_dist_km: Vec<f64>,
    regions: Vec<(Vec<Cell>, Province)>,
    events: Vec<Event>,
    rng: StdRng,
    /// Smooth noise for varying belt widths along their length.
    width_noise: crate::noise::TerrainNoise,
}

impl B<'_> {
    fn km(&self, a: Cell, b: Cell) -> f64 {
        self.grid.center(a).angle_to(self.grid.center(b)) * self.r
    }

    /// Cells within `max_km` (geodesic, along the grid) of `sources`, with distances.
    fn reach(&self, sources: &[Cell], max_km: f64, allow: impl Fn(Cell) -> bool) -> HashMap<Cell, f64> {
        let mut best: HashMap<Cell, f64> = HashMap::new();
        let mut heap = BinaryHeap::new();
        for &s in sources {
            if best.insert(s, 0.0).is_none() {
                heap.push(Item(0.0, s));
            }
        }
        while let Some(Item(d, c)) = heap.pop() {
            if best.get(&c).is_some_and(|&b| d > b) {
                continue;
            }
            for n in self.grid.neighbors(c).iter() {
                if !allow(n) {
                    continue;
                }
                let nd = d + self.km(c, n);
                if nd <= max_km && best.get(&n).is_none_or(|&b| nd < b) {
                    best.insert(n, nd);
                    heap.push(Item(nd, n));
                }
            }
        }
        best
    }

    fn cells_of(m: &HashMap<Cell, f64>, keep: impl Fn(Cell, f64) -> bool) -> Vec<Cell> {
        let mut v: Vec<Cell> = m.iter().filter(|(c, d)| keep(**c, **d)).map(|(c, _)| *c).collect();
        v.sort_unstable();
        v
    }

    fn push(&mut self, cells: Vec<Cell>, mut prov: Province, boundary: Option<usize>, age_ma: f64) {
        if cells.is_empty() {
            return;
        }
        prov.boundary = boundary;
        prov.age_ma = age_ma;
        self.regions.push((cells, prov));
    }

    fn area(&self, cells: &[Cell]) -> f64 {
        cells.iter().map(|&c| self.grid.area_sr(c)).sum::<f64>() * self.r * self.r
    }

    fn sides(&self, s: &BoundarySegment) -> (Side, Side) {
        let mk = |plate: u32| {
            let cells: Vec<Cell> = s.cells.iter().copied().filter(|&c| self.owner[c as usize] == plate).collect();
            let n = cells.len().max(1) as f64;
            let cont = cells.iter().filter(|&&c| self.crust[c as usize] == CrustType::Continental).count() as f64;
            let mean_age = self.age.map_or(0.0, |a| cells.iter().map(|&c| f64::from(a[c as usize])).sum::<f64>() / n);
            Side { plate, continental: cont / n >= 0.5, mean_age }
        };
        (mk(s.a.0 as u32), mk(s.b.0 as u32))
    }

    fn plate_area(&self, plate: u32) -> f64 {
        self.info.get(&plate).map_or(0.0, |i| i.area_km2)
    }

    // ---------------------------------------------------------------- foundation layers

    fn oceanic_base(&mut self) {
        let items: Vec<(Vec<Cell>, f64)> = self.info.values().filter(|i| !i.ocean.is_empty()).map(|i| (i.ocean.clone(), i.area_km2)).collect();
        for (cells, area) in items {
            let age = self.age.map_or(80.0, |a| cells.iter().map(|&c| f64::from(a[c as usize])).sum::<f64>() / cells.len() as f64);
            self.push(cells, pc::abyssal_plain(area), None, age);
        }
    }

    fn stable_continental_base(&mut self) {
        let plates: Vec<u32> = self.info.iter().filter(|(_, i)| !i.cont.is_empty()).map(|(&k, _)| k).collect();
        for k in plates {
            let (cont, cont_area) = {
                let i = &self.info[&k];
                (i.cont.clone(), i.cont_area_km2)
            };
            self.push(cont.clone(), pc::platform(cont_area), None, 500.0);
            // Shield: a noisy spherical cap around the continental centroid.
            let Some(centroid) = UnitVec3::centroid(cont.iter().map(|&c| self.grid.center(c))) else { continue };
            let center = *cont
                .iter()
                .min_by(|&&a, &&b| centroid.angle_to(self.grid.center(a)).total_cmp(&centroid.angle_to(self.grid.center(b))))
                .unwrap_or(&cont[0]);
            let pc0 = self.grid.center(center);
            let target = cont_area * pc::SHIELD_AREA_FRACTION;
            let theta = (1.0 - target / (2.0 * std::f64::consts::PI * self.r * self.r)).clamp(-1.0, 1.0).acos();
            let shield: Vec<Cell> = cont
                .iter()
                .copied()
                .filter(|&c| {
                    let q = self.grid.center(c);
                    let b = pc0.bearing_to(q).to_radians();
                    let noise = (b * 3.0).sin() * 0.2 + (b * 7.0).cos() * 0.1;
                    pc0.angle_to(q) < theta * (1.0 + noise)
                })
                .collect();
            let max_age = self.age.map_or(2500.0, |a| cont.iter().map(|&c| f64::from(a[c as usize])).fold(0.0, f64::max));
            let sa = self.area(&shield);
            let mut prov = pc::craton(sa);
            prov.anchor = Some(pc0);
            self.push(shield, prov, None, max_age);
            // Intracratonic basin.
            if cont_area > pc::INTRACRATONIC_BASIN_MIN_AREA_KM2 && self.rng.r#gen::<f64>() < pc::INTRACRATONIC_BASIN_PROBABILITY {
                let c0 = cont[self.rng.gen_range(0..cont.len())];
                let basin = Self::cells_of(&self.reach(&[c0], pc::INTRACRATONIC_BASIN_RADIUS_KM, |c| self.owner[c as usize] == k && self.crust[c as usize] == CrustType::Continental), |_, _| true);
                let ba = self.area(&basin);
                self.push(basin, pc::intracratonic_basin(ba), None, 400.0);
            }
        }
    }

    fn passive_margins(&mut self, segments: &[BoundarySegment]) {
        for (i, s) in segments.iter().enumerate() {
            let (a, b) = self.sides(s);
            if !a.continental && !b.continental {
                continue;
            }
            let w = if s.kind == BoundaryKind::Divergent { pc::PASSIVE_MARGIN_DIVERGENT_WIDTH_KM } else { pc::PASSIVE_MARGIN_OTHER_WIDTH_KM };
            let m = self.reach(&s.cells, w, |_| true);
            let cells = Self::cells_of(&m, |c, _| self.crust[c as usize] == CrustType::Continental);
            self.push(cells, pc::extended_crust(s.length_km), Some(i), 100.0);
        }
        // Coastlines inside a plate (authored or mixed-crust plates) are passive margins too.
        let interior_edges: Vec<Cell> = self
            .grid
            .cells()
            .filter(|&c| {
                self.crust[c as usize] == CrustType::Continental
                    && self.grid.neighbors(c).iter().any(|n| self.owner[n as usize] == self.owner[c as usize] && self.crust[n as usize] == CrustType::Oceanic)
            })
            .collect();
        if !interior_edges.is_empty() {
            let m = self.reach(&interior_edges, 250.0, |c| self.crust[c as usize] == CrustType::Continental);
            let cells = Self::cells_of(&m, |_, _| true);
            let len = cells.len() as f64;
            self.push(cells, pc::extended_crust(len), None, 100.0);
        }
    }

    fn paleo_orogens(&mut self) {
        let plates: Vec<(u32, f64, f64)> = self.info.iter().map(|(&k, i)| (k, i.cont_area_km2, i.area_km2)).collect();
        for (k, cont_area, area) in plates {
            if cont_area < pc::PALEO_OROGEN_MIN_AREA_KM2 {
                continue;
            }
            let n = if area > pc::HOTSPOT_LARGE_PLATE_THRESHOLD_KM2 { pc::PALEO_OROGEN_COUNT_LARGE } else { pc::PALEO_OROGEN_COUNT_MEDIUM };
            let cont = self.info[&k].cont.clone();
            for _ in 0..n {
                let mut ends = None;
                for _ in 0..10 {
                    let a = cont[self.rng.gen_range(0..cont.len())];
                    let b = cont[self.rng.gen_range(0..cont.len())];
                    if self.km(a, b) > pc::PALEO_OROGEN_MIN_LENGTH_KM {
                        ends = Some((a, b));
                        break;
                    }
                }
                let width = pc::PALEO_OROGEN_MIN_WIDTH_KM + self.rng.r#gen::<f64>() * pc::PALEO_OROGEN_WIDTH_RANGE_KM;
                let age = 300.0 + self.rng.r#gen::<f64>() * 900.0;
                let Some((a, b)) = ends else { continue };
                let bow = (self.rng.r#gen::<f64>() - 0.5) * 0.5;
                let line = self.bowed_path(a, b, bow);
                let m = self.reach(&line, width / 2.0, |c| self.owner[c as usize] == k && self.crust[c as usize] == CrustType::Continental);
                let cells = Self::cells_of(&m, |c, _| self.owner[c as usize] == k && self.crust[c as usize] == CrustType::Continental);
                if self.area(&cells) < 50_000.0 {
                    continue;
                }
                let mid = UnitVec3::centroid(cells.iter().map(|&c| self.grid.center(c))).unwrap_or(UnitVec3::NORTH);
                self.events.push(Event { t_ma: age.round(), kind: EventKind::PaleoOrogeny, at: Place::Point { at: mid }, magnitude: width, causes: vec![], note: "ancient mountain belt, now eroded".into() });
                self.push(cells, pc::paleo_orogen(width), None, age);
            }
        }
    }

    /// A great-circle path bowed sideways by `bow` × its length at the middle (belts curve).
    fn bowed_path(&self, a: Cell, b: Cell, bow: f64) -> Vec<Cell> {
        let pa = self.grid.center(a);
        let pb = self.grid.center(b);
        let ang = pa.angle_to(pb);
        let axis = pa.v().cross(pb.v()).normalized().unwrap_or(UnitVec3::NORTH);
        let steps = ((ang / (self.grid.spacing_rad() * 0.5)).ceil() as usize).max(1);
        let mut v: Vec<Cell> = (0..=steps)
            .map(|i| {
                let t = i as f64 / steps as f64;
                let s = ang.sin().max(1e-12);
                let q = (pa.v() * (((1.0 - t) * ang).sin() / s) + pb.v() * ((t * ang).sin() / s)).normalized().unwrap_or(pa);
                let off = bow * ang * (std::f64::consts::PI * t).sin();
                let q = (q.v() * off.cos() + axis.v() * off.sin()).normalized().unwrap_or(q);
                self.grid.locate(q)
            })
            .collect();
        v.dedup();
        v
    }

    #[allow(dead_code)]
    fn great_circle(&self, a: Cell, b: Cell) -> Vec<Cell> {
        let pa = self.grid.center(a);
        let pb = self.grid.center(b);
        let ang = pa.angle_to(pb);
        let steps = ((ang / (self.grid.spacing_rad() * 0.5)).ceil() as usize).max(1);
        let mut v: Vec<Cell> = (0..=steps)
            .map(|i| {
                let t = i as f64 / steps as f64;
                let s = ang.sin().max(1e-12);
                let wa = ((1.0 - t) * ang).sin() / s;
                let wb = (t * ang).sin() / s;
                let p = (pa.v() * wa + pb.v() * wb).normalized().unwrap_or(pa);
                self.grid.locate(p)
            })
            .collect();
        v.dedup();
        v
    }

    // ---------------------------------------------------------------- active features

    fn collision_orogens(&mut self, segments: &[BoundarySegment]) {
        for (i, s) in segments.iter().enumerate() {
            if s.kind != BoundaryKind::Convergent || s.rate_cm_yr < pc::MIN_CONVERGENCE_RATE_CM_PER_YEAR {
                continue;
            }
            let (a, b) = self.sides(s);
            if !(a.continental && b.continental) {
                continue;
            }
            let width = pc::rate_scaled(pc::COLLISION_BASE_WIDTH_KM, pc::COLLISION_WIDTH_FACTOR, pc::COLLISION_MAX_MULTIPLIER, s.rate_cm_yr);
            // Belts swell and pinch along their length (half to one and a half times the
            // nominal width), so their edges are not parallel to a straight boundary.
            let m = self.reach(&s.cells, width * 0.75, |c| self.crust[c as usize] == CrustType::Continental);
            let wn = &self.width_noise;
            let grid = self.grid;
            let cells = Self::cells_of(&m, |c, d| d <= width / 2.0 * (1.0 + 1.6 * wn.fbm(grid.center(c), 6.0, 3)).clamp(0.5, 1.5));
            let age = (1500.0 / (s.rate_cm_yr * 10.0).max(5.0)).clamp(5.0, 80.0);
            self.push(cells, pc::collision_orogen(s.rate_cm_yr, width), Some(i), age);
        }
    }

    fn arc_systems(&mut self, segments: &[BoundarySegment]) {
        for (i, s) in segments.iter().enumerate() {
            if s.kind != BoundaryKind::Convergent {
                continue;
            }
            let (a, b) = self.sides(s);
            let (sub, over) = match (a.continental, b.continental) {
                (true, true) => continue,
                (false, true) => (a, b),
                (true, false) => (b, a),
                (false, false) => {
                    // The older (denser) oceanic crust sinks.
                    if a.mean_age > b.mean_age || (a.mean_age == b.mean_age && a.plate < b.plate) { (a, b) } else { (b, a) }
                }
            };
            let rate = s.rate_cm_yr;
            let age = (1500.0 / (rate * 10.0).max(5.0)).clamp(5.0, 120.0);
            // Trench: on the subducting plate, at the boundary.
            let t = self.reach(&s.cells, pc::OCEAN_TRENCH_WIDTH_KM / 2.0, |_| true);
            let trench = Self::cells_of(&t, |c, _| self.owner[c as usize] == sub.plate);
            self.push(trench, pc::ocean_trench(rate, s.length_km), Some(i), age);
            // Wedge → forearc → arc → backarc, stepping into the overriding plate.
            let w1 = pc::rate_scaled(pc::ACCRETIONARY_WEDGE_BASE_WIDTH_KM, pc::ACCRETIONARY_WEDGE_WIDTH_FACTOR, pc::ACCRETIONARY_WEDGE_MAX_MULTIPLIER, rate);
            let w2 = pc::rate_scaled(pc::FOREARC_BASIN_BASE_WIDTH_KM, pc::FOREARC_BASIN_WIDTH_FACTOR, pc::FOREARC_BASIN_MAX_MULTIPLIER, rate);
            let w3 = pc::rate_scaled(pc::VOLCANIC_ARC_BASE_WIDTH_KM, pc::VOLCANIC_ARC_WIDTH_FACTOR, pc::VOLCANIC_ARC_MAX_MULTIPLIER, rate);
            let w4 = pc::rate_scaled(pc::BACKARC_BASIN_BASE_WIDTH_KM, pc::BACKARC_BASIN_WIDTH_FACTOR, pc::BACKARC_BASIN_MAX_MULTIPLIER, rate);
            let backarc = self.plate_area(a.plate) > pc::BACKARC_BASIN_MIN_PLATE_AREA_KM2 || self.plate_area(b.plate) > pc::BACKARC_BASIN_MIN_PLATE_AREA_KM2;
            let total = w1 + w2 + w3 + if backarc { w4 } else { 0.0 };
            let sources: Vec<Cell> = s.cells.iter().copied().filter(|&c| self.owner[c as usize] == over.plate).collect();
            let m = self.reach(&sources, total, |c| self.owner[c as usize] == over.plate);
            let band = |lo: f64, hi: f64| Self::cells_of(&m, |_, d| d > lo && d <= hi);
            let wedge = Self::cells_of(&m, |_, d| d <= w1);
            let fore = band(w1, w1 + w2);
            let arc = band(w1 + w2, w1 + w2 + w3);
            self.push(wedge, pc::accretionary_wedge(rate, w1), Some(i), age);
            self.push(fore, pc::forearc_basin(s.length_km), Some(i), age);
            self.push(arc, pc::volcanic_arc(rate, s.length_km), Some(i), age);
            if backarc {
                let ba = band(w1 + w2 + w3, total);
                self.push(ba, pc::backarc_basin(s.length_km), Some(i), age * 0.5);
            }
        }
    }

    fn oceanic_overlays(&mut self, segments: &[BoundarySegment]) {
        for (i, s) in segments.iter().enumerate() {
            let (a, b) = self.sides(s);
            if s.kind == BoundaryKind::Divergent && !a.continental && !b.continental {
                let r = s.rate_cm_yr;
                let half = if r > pc::SPREADING_RATE_FAST {
                    pc::RIDGE_HALF_WIDTH_FAST_KM
                } else if r > pc::SPREADING_RATE_MEDIUM {
                    pc::RIDGE_HALF_WIDTH_MEDIUM_KM
                } else if r > pc::SPREADING_RATE_SLOW {
                    pc::RIDGE_HALF_WIDTH_SLOW_KM
                } else {
                    pc::RIDGE_HALF_WIDTH_ULTRASLOW_KM
                };
                let cells = Self::cells_of(&self.reach(&s.cells, half, |_| true), |_, _| true);
                self.push(cells, pc::mid_ocean_ridge(r, s.length_km), Some(i), 0.0);
            }
        }
        for (i, s) in segments.iter().enumerate() {
            let (a, b) = self.sides(s);
            if s.kind == BoundaryKind::Transform && (!a.continental || !b.continental) {
                let cells = Self::cells_of(&self.reach(&s.cells, pc::FRACTURE_ZONE_HALF_WIDTH_KM, |_| true), |_, _| true);
                let age = (s.length_km * 0.3 / (s.rate_cm_yr * 10.0).max(5.0)).clamp(2.0, 60.0);
                self.push(cells, pc::oceanic_fracture_zone(s.length_km), Some(i), age);
            }
        }
    }

    fn continental_rifts(&mut self, segments: &[BoundarySegment]) {
        for (i, s) in segments.iter().enumerate() {
            if s.kind != BoundaryKind::Divergent {
                continue;
            }
            let (a, b) = self.sides(s);
            if a.continental && b.continental {
                let half = pc::RIFT_MIN_HALF_WIDTH_KM + self.rng.r#gen::<f64>() * pc::RIFT_HALF_WIDTH_RANGE_KM;
                let cells = Self::cells_of(&self.reach(&s.cells, half, |_| true), |_, _| true);
                let age = (300.0 / (s.rate_cm_yr * 10.0).max(5.0)).clamp(2.0, 40.0);
                self.push(cells, pc::continental_rift(s.rate_cm_yr, s.length_km), Some(i), age);
            }
        }
    }

    fn large_igneous_provinces(&mut self, probability: f64) {
        let keys: Vec<u32> = self.info.keys().copied().collect();
        for k in keys {
            if self.rng.r#gen::<f64>() >= probability {
                continue;
            }
            let (kind, area, pool) = {
                let i = &self.info[&k];
                if i.kind == CrustType::Continental {
                    (GeologicProvince::ContinentalFloodBasalt, i.cont_area_km2 * pc::CONTINENTAL_FLOOD_BASALT_AREA_FRACTION, i.cont.clone())
                } else {
                    (GeologicProvince::OceanicPlateau, (i.area_km2 - i.cont_area_km2) * pc::OCEANIC_PLATEAU_AREA_FRACTION, i.ocean.clone())
                }
            };
            let area = area.min(pc::LIP_MAX_AREA_KM2);
            let interior: Vec<Cell> = pool.iter().copied().filter(|&c| self.edge_dist_km[c as usize] > 300.0).collect();
            let age = 15.0 + self.rng.r#gen::<f64>() * 235.0;
            if interior.is_empty() || area < 50_000.0 {
                continue;
            }
            let c0 = interior[self.rng.gen_range(0..interior.len())];
            let radius = (area / std::f64::consts::PI).sqrt();
            let crust0 = self.crust[c0 as usize];
            let m = self.reach(&[c0], radius * 1.3, |c| self.owner[c as usize] == k && self.crust[c as usize] == crust0);
            let p0 = self.grid.center(c0);
            let cells = Self::cells_of(&m, |c, d| {
                let b = p0.bearing_to(self.grid.center(c)).to_radians();
                d < radius * (1.0 + 0.25 * (b * 3.0 + 1.0).sin() + 0.1 * (b * 5.0).cos())
            });
            let a = self.area(&cells);
            let mut prov = if kind == GeologicProvince::ContinentalFloodBasalt { pc::continental_flood_basalt(a) } else { pc::oceanic_plateau(a) };
            prov.anchor = Some(p0);
            self.events.push(Event { t_ma: age.round(), kind: EventKind::LargeIgneousProvince, at: Place::Point { at: p0 }, magnitude: a, causes: vec![], note: format!("{} erupts", kind_name(kind)) });
            self.push(cells, prov, None, age);
        }
    }

    fn hotspot_tracks(&mut self) {
        let keys: Vec<(u32, f64)> = self.info.iter().map(|(&k, i)| (k, i.area_km2)).collect();
        for (k, area) in keys {
            if area < pc::HOTSPOT_MIN_PLATE_AREA_KM2 {
                continue;
            }
            let prob = if area > pc::HOTSPOT_LARGE_PLATE_THRESHOLD_KM2 { pc::HOTSPOT_PROBABILITY_LARGE_PLATE } else { pc::HOTSPOT_PROBABILITY_MEDIUM_PLATE };
            if self.rng.r#gen::<f64>() > prob {
                continue;
            }
            let min_d = if area > pc::HOTSPOT_LARGE_PLATE_THRESHOLD_KM2 { pc::HOTSPOT_MIN_INTERIOR_DISTANCE_LARGE_KM } else { pc::HOTSPOT_MIN_INTERIOR_DISTANCE_MEDIUM_KM };
            let interior: Vec<Cell> = self.grid.cells().filter(|&c| self.owner[c as usize] == k && self.edge_dist_km[c as usize] >= min_d).collect();
            if interior.is_empty() {
                continue;
            }
            let h = interior[self.rng.gen_range(0..interior.len())];
            let oceanic = self.crust[h as usize] == CrustType::Oceanic;
            let (t_span, max_len, half_w) = if oceanic {
                (pc::HOTSPOT_OCEANIC_MIN_AGE_MA + self.rng.r#gen::<f64>() * pc::HOTSPOT_OCEANIC_AGE_RANGE_MA, pc::HOTSPOT_OCEANIC_MAX_LENGTH_KM, pc::HOTSPOT_OCEANIC_HALF_WIDTH_KM)
            } else {
                (pc::HOTSPOT_CONTINENTAL_MIN_AGE_MA + self.rng.r#gen::<f64>() * pc::HOTSPOT_CONTINENTAL_AGE_RANGE_MA, pc::HOTSPOT_CONTINENTAL_MAX_LENGTH_KM, pc::HOTSPOT_CONTINENTAL_HALF_WIDTH_KM)
            };
            let Some(plate) = self.plates.table.get(Key(u64::from(k))) else { continue };
            let ph = self.grid.center(h);
            let speed_km_myr = plate.velocity_km_per_myr(ph, self.r).norm();
            let length = (speed_km_myr * t_span).min(max_len);
            if length < pc::HOTSPOT_MIN_VISIBLE_LENGTH_KM || speed_km_myr <= 0.0 {
                continue;
            }
            // Volcanoes formed t Myr ago have since ridden the plate's rotation away from
            // the fixed hotspot, so the track is the hotspot carried forward by the plate.
            let t_len = length / speed_km_myr;
            let step_myr = (self.grid.spacing_rad() * 0.5 * self.r / speed_km_myr).max(1e-3);
            let mut track = vec![h];
            let mut t = step_myr;
            while t <= t_len {
                let q = ph.rotate_about(plate.euler_pole, (plate.rate_deg_per_myr * t).to_radians());
                let c = self.grid.locate(q);
                if self.owner[c as usize] != k {
                    break;
                }
                if track.last() != Some(&c) {
                    track.push(c);
                }
                t += step_myr;
            }
            let cells = Self::cells_of(&self.reach(&track, half_w, |_| true), |_, _| true);
            let mut prov = if oceanic { pc::oceanic_hotspot_track(length) } else { pc::continental_hotspot_track(length) };
            prov.anchor = Some(ph);
            self.events.push(Event { t_ma: (t_len * 10.0).round() / 10.0, kind: EventKind::HotspotTrack, at: Place::Point { at: ph }, magnitude: length, causes: vec![], note: "hotspot track begins; active end is the youngest".into() });
            self.push(cells, prov, None, t_len);
        }
    }
}

fn kind_name(k: GeologicProvince) -> &'static str {
    use geoforge_core::surface::Category;
    k.name()
}

/// Build provinces and the geology timeline (v1's layer order: oldest first, overlays last).
#[allow(clippy::too_many_arguments)]
#[must_use]
pub fn build(
    grid: &dyn Grid,
    plates: &Partition<Plate>,
    crust: &[CrustType],
    age: Option<&[f32]>,
    segments: &[BoundarySegment],
    radius_km: f64,
    seed: u64,
    p: &GeologyParams,
) -> (Partition<Province>, Timeline) {
    let owner = &plates.owner.values;
    let r2 = radius_km * radius_km;
    let mut info: BTreeMap<u32, PlateInfo> = BTreeMap::new();
    for c in grid.cells() {
        let o = owner[c as usize];
        if o == 0 {
            continue;
        }
        let e = info.entry(o).or_insert(PlateInfo { area_km2: 0.0, cont: vec![], ocean: vec![], cont_area_km2: 0.0, kind: CrustType::Oceanic });
        let a = grid.area_sr(c) * r2;
        e.area_km2 += a;
        if crust[c as usize] == CrustType::Continental {
            e.cont.push(c);
            e.cont_area_km2 += a;
        } else {
            e.ocean.push(c);
        }
    }
    for i in info.values_mut() {
        i.kind = if i.cont_area_km2 * 2.0 >= i.area_km2 { CrustType::Continental } else { CrustType::Oceanic };
    }
    let edges: Vec<Cell> = grid.cells().filter(|&c| grid.neighbors(c).iter().any(|n| owner[n as usize] != owner[c as usize])).collect();
    let edge_dist_km: Vec<f64> = distance_field(grid, edges, std::f64::consts::PI, |_| true).into_iter().map(|d| d * radius_km).collect();

    let mut b = B { grid, owner, crust, age, r: radius_km, plates, info, edge_dist_km, regions: vec![], events: vec![], rng: StdRng::seed_from_u64(seed), width_noise: crate::noise::TerrainNoise::new(seed ^ 0x0B1D_7A11) };

    if p.oceanic {
        b.oceanic_base();
    }
    if p.stable_regions {
        b.stable_continental_base();
    }
    if p.extensional {
        b.passive_margins(segments);
    }
    if p.stable_regions {
        b.paleo_orogens();
    }
    b.collision_orogens(segments);
    if p.arc_systems {
        b.arc_systems(segments);
    }
    if p.oceanic {
        b.oceanic_overlays(segments);
    }
    if p.extensional {
        b.continental_rifts(segments);
    }
    b.large_igneous_provinces(p.lip_probability);
    if p.oceanic {
        b.hotspot_tracks();
    }

    // Paint in order: later regions overwrite earlier ones.
    let mut owner_out = vec![0u32; grid.len()];
    let mut table = Table::new();
    for (i, (cells, prov)) in b.regions.into_iter().enumerate() {
        let key = i as u32 + 1;
        for c in cells {
            owner_out[c as usize] = key;
        }
        table.insert(Key(u64::from(key)), prov);
    }
    let mut part = Partition { owner: Field::from_values(grid.spec(), owner_out), table };
    part.prune();

    let mut events = Table::new();
    let mut evs = b.events;
    evs.sort_by(|a, b| b.t_ma.total_cmp(&a.t_ma).then(a.kind.cmp(&b.kind)));
    for (i, e) in evs.into_iter().enumerate() {
        events.insert(Key(i as u64 + 1), e);
    }
    (part, Timeline { events })
}
