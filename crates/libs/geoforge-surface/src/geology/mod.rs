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
//! - Continental collisions are asymmetric (v1 drew a ridge centred on the boundary): the
//!   smaller continent goes under (unless the other rides over an ocean slab nearby), a
//!   front range and a plateau grow on the upper plate, and a foreland basin bends down in
//!   front (see [`provinces::Collision`]).
//! - Ocean–continent margins are not built like island arcs: a narrow wedge under the
//!   slope, an onshore forearc, a broad cordillera, a plateau behind it over an old slab,
//!   and a retroarc foreland instead of an extensional backarc (see
//!   [`provinces::ANDEAN_WEDGE_BASE_WIDTH_KM`]).
//! - Rifts do not cut across belts already built at convergent margins.

pub mod lithology;
pub mod provinces;

use crate::elevation::smoothstep;
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
    /// Rocks added to Earth's catalog (see [`lithology::RockDef`]).
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub extra_rocks: Vec<lithology::RockDef>,
}

impl Default for GeologyParams {
    fn default() -> Self {
        Self { lip_probability: pc::LIP_PROBABILITY, boundaries: BoundaryParams::default(), oceanic: true, stable_regions: true, arc_systems: true, extensional: true, extra_rocks: Vec::new() }
    }
}

pub struct Geology;

impl Generator for Geology {
    type Params = GeologyParams;
    const NAME: &'static str = "geology";
    const VERSION: u32 = 6;
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
            let (provinces, timeline) = build(ctx.grid, plates, &crust.values, inputs.crust_age().map(|f| f.values.as_slice()), &segments, ctx.radius_km(), ctx.seed(ProductId::Provinces).value(), p);
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
            let catalog = lithology::Catalog::new(&p.extra_rocks)?;
            let (lith, age) = lithology::build(ctx.grid, &provinces, &crust.values, inputs.crust_age().map(|f| f.values.as_slice()), ctx.seed(ProductId::Lithology).value(), &catalog);
            if inputs.lithology().is_none() {
                out.insert(ProductData::Lithology(lith));
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

/// How far from an arc's trench the upper plate's crust counts as hot and weak (the
/// central Andes are ~600 km across from the trench to the Subandean front).
const ARC_WEAKENED_KM: f64 = 600.0;

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
    /// Each plate's craton: its centre and cells.
    cratons: BTreeMap<u32, (UnitVec3, Vec<Cell>)>,
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

    /// Like [`Self::reach`], but also returns each cell's nearest source.
    fn nearest(&self, sources: &[Cell], max_km: f64, allow: impl Fn(Cell) -> bool) -> HashMap<Cell, (f64, Cell)> {
        let mut best: HashMap<Cell, (f64, Cell)> = HashMap::new();
        let mut heap = BinaryHeap::new();
        for &s in sources {
            if best.insert(s, (0.0, s)).is_none() {
                heap.push(Item(0.0, s));
            }
        }
        while let Some(Item(d, c)) = heap.pop() {
            let Some(&(b, src)) = best.get(&c) else { continue };
            if d > b {
                continue;
            }
            for n in self.grid.neighbors(c).iter() {
                if !allow(n) {
                    continue;
                }
                let nd = d + self.km(c, n);
                if nd <= max_km && best.get(&n).is_none_or(|&(b, _)| nd < b) {
                    best.insert(n, (nd, src));
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

    /// Each plate's side of a boundary, judged over the first 200 km into the plate: at an
    /// active margin the trench lies offshore of a steep continental slope, so the cells on
    /// the boundary itself can be deep water on a continental plate (Peru, Chile).
    fn sides(&self, s: &BoundarySegment) -> (Side, Side) {
        const SIDE_KM: f64 = 200.0;
        let mk = |plate: u32| {
            let edge: Vec<Cell> = s.cells.iter().copied().filter(|&c| self.owner[c as usize] == plate).collect();
            let cells = Self::cells_of(&self.reach(&edge, SIDE_KM, |c| self.owner[c as usize] == plate), |_, _| true);
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
            let center = *cont.iter().min_by(|&&a, &&b| centroid.angle_to(self.grid.center(a)).total_cmp(&centroid.angle_to(self.grid.center(b)))).unwrap_or(&cont[0]);
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
            self.cratons.insert(k, (pc0, shield.clone()));
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

    /// Old mountain belts, worn down, where blocks were stitched onto a craton. A continent
    /// is its craton plus what accreted to it, so the belts run through that accreted
    /// ground, between the craton and the coast and parallel to both (the Appalachians
    /// between the Canadian Shield and the Atlantic, the Tasmanides east of Australia's
    /// cratons): each at a share of the way out, on the side where the most ground lies
    /// beyond the craton; a large continent's second on another side, older and closer in.
    /// Belts wander, swell and pinch along their length, thin out at their ends, break where
    /// there is no room, and never cross a craton or each other.
    fn paleo_orogens(&mut self) {
        const BINS: usize = 72;
        let plates: Vec<(u32, f64)> = self.info.iter().filter(|(_, i)| i.cont_area_km2 >= pc::PALEO_OROGEN_MIN_AREA_KM2).map(|(&k, i)| (k, i.area_km2)).collect();
        let mut taken: std::collections::HashSet<Cell> = std::collections::HashSet::new();
        for (k, area) in plates {
            let Some((centre, craton)) = self.cratons.get(&k).cloned() else { continue };
            let in_craton: std::collections::HashSet<Cell> = craton.iter().copied().collect();
            let (owner, crust, grid) = (self.owner, self.crust, self.grid);
            let continental = |c: Cell| owner[c as usize] == k && crust[c as usize] == CrustType::Continental;
            let rim: Vec<Cell> = craton.iter().copied().filter(|&c| grid.neighbors(c).iter().any(|n| continental(n) && !in_craton.contains(&n))).collect();
            if rim.is_empty() {
                continue;
            }
            let out = self.reach(&rim, f64::INFINITY, |c| continental(c) && !in_craton.contains(&c));
            let bin = |c: Cell| ((centre.bearing_to(grid.center(c)).rem_euclid(360.0) / 360.0 * BINS as f64) as usize).min(BINS - 1);
            // How far the accreted ground reaches beyond the craton, by direction, smoothed.
            // Ground running out to an active margin is that margin's foreland (the Amazon
            // behind the Andes) and counts for little; old belts face passive coasts.
            let mut reach = [(0.0f64, Cell::MAX); BINS];
            for (&c, &d) in &out {
                let r = &mut reach[bin(c)];
                if d > r.0 || (d == r.0 && c < r.1) {
                    *r = (d, c);
                }
            }
            let reach: Vec<f64> = reach.iter().map(|&(d, c)| if c != Cell::MAX && self.edge_dist_km[c as usize] < pc::PALEO_OROGEN_ACTIVE_MARGIN_KM { d * 0.25 } else { d }).collect();
            let room: Vec<f64> = (0..BINS).map(|b| (0..5).map(|o| reach[(b + BINS + o - 2) % BINS]).sum::<f64>() / 5.0).collect();
            let n = if area > pc::HOTSPOT_LARGE_PLATE_THRESHOLD_KM2 { pc::PALEO_OROGEN_COUNT_LARGE } else { pc::PALEO_OROGEN_COUNT_MEDIUM };
            let mut used: Vec<usize> = Vec::new();
            for j in 0..n {
                let half = pc::PALEO_OROGEN_MIN_HALF_ARC_DEG + self.rng.r#gen::<f64>() * pc::PALEO_OROGEN_HALF_ARC_RANGE_DEG;
                let hb = (half / 360.0 * BINS as f64).round() as usize;
                let apart = |b: usize| used.iter().all(|&u| (b + BINS - u) % BINS > BINS / 4 && (u + BINS - b) % BINS > BINS / 4);
                let sector_room = |b: usize| (0..=2 * hb).map(|o| room[(b + BINS + o - hb) % BINS].min(pc::PALEO_OROGEN_MAX_ROOM_KM)).sum::<f64>();
                let Some(toward) = (0..BINS).filter(|&b| apart(b)).max_by(|&a, &b| sector_room(a).total_cmp(&sector_room(b)).then(b.cmp(&a))) else { break };
                used.push(toward);
                let toward_deg = (toward as f64 + 0.5) * 360.0 / BINS as f64;
                let (lo, hi) = if j == 0 { pc::PALEO_OROGEN_OUTER_SHARE } else { pc::PALEO_OROGEN_INNER_SHARE };
                let share = lo + self.rng.r#gen::<f64>() * (hi - lo);
                let width = pc::PALEO_OROGEN_MIN_WIDTH_KM + self.rng.r#gen::<f64>() * pc::PALEO_OROGEN_WIDTH_RANGE_KM;
                let age = if j == 0 { 300.0 + self.rng.r#gen::<f64>() * 400.0 } else { 700.0 + self.rng.r#gen::<f64>() * 500.0 };
                // Two smooth fields along the belt: where its axis lies, and how wide it is.
                let wander = crate::noise::TerrainNoise::new(self.rng.r#gen());
                let swell = crate::noise::TerrainNoise::new(self.rng.r#gen());
                let mut cells: Vec<Cell> = out
                    .iter()
                    .filter(|(c, _)| !taken.contains(*c))
                    .filter_map(|(&c, &d)| {
                        let q = grid.center(c);
                        let off = ((centre.bearing_to(q) - toward_deg + 540.0).rem_euclid(360.0) - 180.0).abs();
                        let r = room[bin(c)];
                        if off >= half || r < pc::PALEO_OROGEN_MIN_ROOM_KM {
                            return None;
                        }
                        let ends = smoothstep(0.0, 0.25, 1.0 - off / half).sqrt();
                        let axis = share * r * (1.0 + 0.4 * wander.fbm(q, 4.0, 2));
                        // Swells to half again its width, pinches to half (never wider than the
                        // ground allows), and breaks where the noise is lowest.
                        let sw = swell.fbm(q, 6.0, 2);
                        let w = if sw < -0.35 { 0.0 } else { (width * (1.0 + 1.5 * sw).clamp(0.5, 1.5)).min(0.6 * r) * ends };
                        ((d - axis).abs() <= w / 2.0).then_some(c)
                    })
                    .collect();
                cells.sort_unstable();
                if self.area(&cells) < 50_000.0 {
                    continue;
                }
                taken.extend(cells.iter().copied());
                let mid = UnitVec3::centroid(cells.iter().map(|&c| self.grid.center(c))).unwrap_or(UnitVec3::NORTH);
                self.events.push(Event {
                    t_ma: age.round(),
                    kind: EventKind::PaleoOrogeny,
                    at: Place::Point { at: mid },
                    magnitude: width,
                    causes: vec![],
                    note: "ancient mountain belt, now eroded".into(),
                });
                self.push(cells, pc::paleo_orogen(width), None, age);
            }
        }
    }

    // ---------------------------------------------------------------- active features

    /// Of two colliding continents, the one that goes under: the smaller, driven into the
    /// larger, as India and Arabia are into Asia. Polarity is inherited: the small continent
    /// rode in behind ocean floor that was already subducting beneath the big one, and once
    /// its buoyant upper crust is scraped off, the rest is dense enough to follow.
    ///
    /// The exception is a plate already riding over an ocean slab next to the collision:
    /// its crust is hot and weak above the slab, so the other plate's colder lithosphere
    /// drives under it whatever their sizes (in Taiwan, China's margin goes under the Luzon
    /// arc).
    fn lower_and_upper(&self, s: &BoundarySegment, a: Side, b: Side, arc_fronts: &BTreeMap<u32, Vec<Cell>>) -> (Side, Side) {
        let weak = |side: Side| {
            let Some(front) = arc_fronts.get(&side.plate) else { return 0.0 };
            let near = self.reach(front, ARC_WEAKENED_KM, |c| self.owner[c as usize] == side.plate);
            let edge: Vec<Cell> = s.cells.iter().copied().filter(|&c| self.owner[c as usize] == side.plate).collect();
            edge.iter().filter(|c| near.contains_key(c)).count() as f64 / edge.len().max(1) as f64
        };
        let (wa, wb) = (weak(a), weak(b));
        if (wa - wb).abs() > 0.5 {
            return if wa > wb { (b, a) } else { (a, b) };
        }
        let size = |side: Side| (self.info.get(&side.plate).map_or(0.0, |i| i.cont_area_km2), side.plate);
        if size(a) < size(b) {
            (a, b)
        } else {
            (b, a)
        }
    }

    /// The overriding plate and its boundary cells at every subduction under it.
    fn arc_fronts(&self, segments: &[BoundarySegment]) -> BTreeMap<u32, Vec<Cell>> {
        let mut fronts: BTreeMap<u32, Vec<Cell>> = BTreeMap::new();
        for s in segments.iter().filter(|s| s.kind == BoundaryKind::Convergent) {
            if let Some((_, over)) = self.subduction(s) {
                fronts.entry(over.plate).or_default().extend(s.cells.iter().copied().filter(|&c| self.owner[c as usize] == over.plate));
            }
        }
        fronts
    }

    /// The subducting and overriding sides of a convergent boundary with ocean floor on at
    /// least one side (`None` for two continents).
    fn subduction(&self, s: &BoundarySegment) -> Option<(Side, Side)> {
        let (a, b) = self.sides(s);
        match (a.continental, b.continental) {
            (true, true) => None,
            (false, true) => Some((a, b)),
            (true, false) => Some((b, a)),
            // The older (denser) oceanic crust sinks.
            (false, false) => Some(if a.mean_age > b.mean_age || (a.mean_age == b.mean_age && a.plate < b.plate) { (a, b) } else { (b, a) }),
        }
    }

    fn collision_orogens(&mut self, segments: &[BoundarySegment]) {
        let arc_fronts = self.arc_fronts(segments);
        for (i, s) in segments.iter().enumerate() {
            if s.kind != BoundaryKind::Convergent || s.rate_cm_yr < pc::COLLISION_MIN_RATE_CM_PER_YEAR {
                continue;
            }
            let (a, b) = self.sides(s);
            if !(a.continental && b.continental) {
                continue;
            }
            let (lower, upper) = self.lower_and_upper(s, a, b, &arc_fronts);
            let shape = pc::Collision::at_rate(s.rate_cm_yr);
            let age = (1500.0 / (s.rate_cm_yr * 10.0).max(5.0)).clamp(5.0, 80.0);
            let (owner, crust) = (self.owner, self.crust);
            let continental_on = |plate: u32| move |c: Cell| owner[c as usize] == plate && crust[c as usize] == CrustType::Continental;
            let edge = |plate: u32| -> Vec<Cell> { s.cells.iter().copied().filter(|&c| owner[c as usize] == plate).collect() };
            // Belts swell and pinch along their length, so their edges are not parallel to a
            // straight boundary.
            let wn = &self.width_noise;
            let grid = self.grid;
            let swell = |c: Cell, lo: f64, hi: f64| (1.0 + 1.6 * wn.fbm(grid.center(c), 6.0, 3)).clamp(lo, hi);

            // Upper plate: front range, then plateau.
            let depth = shape.range_km + shape.plateau_km;
            let m = self.reach(&edge(upper.plate), depth * 1.2, continental_on(upper.plate));
            let range = Self::cells_of(&m, |c, d| d <= shape.range_km * swell(c, 0.8, 1.2));
            let plateau = Self::cells_of(&m, |c, d| d > shape.range_km * swell(c, 0.8, 1.2) && d <= shape.range_km + shape.plateau_km * swell(c, 0.8, 1.2));
            // Lower plate: foreland basin.
            let f = self.reach(&edge(lower.plate), shape.foreland_km, continental_on(lower.plate));
            let foreland = Self::cells_of(&f, |_, _| true);

            self.push(foreland, pc::foreland_basin(s.rate_cm_yr, shape.foreland_km), Some(i), age * 0.5);
            if shape.plateau_km > 0.0 {
                self.push(plateau, pc::orogenic_plateau(s.rate_cm_yr, shape.plateau_km), Some(i), age);
            }
            self.push(range, pc::collision_orogen(s.rate_cm_yr, shape.range_km), Some(i), age);
        }
    }

    fn arc_systems(&mut self, segments: &[BoundarySegment]) {
        for (i, s) in segments.iter().enumerate() {
            if s.kind != BoundaryKind::Convergent {
                continue;
            }
            let Some((sub, over)) = self.subduction(s) else { continue };
            let rate = s.rate_cm_yr;
            let age = (1500.0 / (rate * 10.0).max(5.0)).clamp(5.0, 120.0);
            // Trench: on the subducting plate, at the boundary.
            let t = self.reach(&s.cells, pc::OCEAN_TRENCH_WIDTH_KM / 2.0, |_| true);
            let trench = Self::cells_of(&t, |c, _| self.owner[c as usize] == sub.plate);
            self.push(trench, pc::ocean_trench(rate, s.length_km), Some(i), age);
            if over.continental {
                self.andean_margin(i, s, sub.plate, over.plate, age);
                continue;
            }
            // Wedge → forearc → arc → backarc, stepping into the overriding plate.
            let w1 = pc::rate_scaled(pc::ACCRETIONARY_WEDGE_BASE_WIDTH_KM, pc::ACCRETIONARY_WEDGE_WIDTH_FACTOR, pc::ACCRETIONARY_WEDGE_MAX_MULTIPLIER, rate);
            let w2 = pc::rate_scaled(pc::FOREARC_BASIN_BASE_WIDTH_KM, pc::FOREARC_BASIN_WIDTH_FACTOR, pc::FOREARC_BASIN_MAX_MULTIPLIER, rate);
            let w3 = pc::rate_scaled(pc::VOLCANIC_ARC_BASE_WIDTH_KM, pc::VOLCANIC_ARC_WIDTH_FACTOR, pc::VOLCANIC_ARC_MAX_MULTIPLIER, rate);
            let w4 = pc::rate_scaled(pc::BACKARC_BASIN_BASE_WIDTH_KM, pc::BACKARC_BASIN_WIDTH_FACTOR, pc::BACKARC_BASIN_MAX_MULTIPLIER, rate);
            let backarc = self.plate_area(sub.plate) > pc::BACKARC_BASIN_MIN_PLATE_AREA_KM2 || self.plate_area(over.plate) > pc::BACKARC_BASIN_MIN_PLATE_AREA_KM2;
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

    /// Wedge → forearc → cordillera → (plateau) → retroarc foreland, stepping into a
    /// continental upper plate (see [`pc::ANDEAN_WEDGE_BASE_WIDTH_KM`]). The plateau stands
    /// only where the slab entering the nearest stretch of trench is old, but not very old.
    /// Plateau and foreland are built on continental crust only.
    fn andean_margin(&mut self, i: usize, s: &BoundarySegment, sub: u32, over: u32, age: f64) {
        let rate = s.rate_cm_yr;
        let scaled = |base: f64| pc::rate_scaled(base, pc::ACCRETIONARY_WEDGE_WIDTH_FACTOR, pc::ACCRETIONARY_WEDGE_MAX_MULTIPLIER, rate);
        let (w1, w2, w3, w5) = (scaled(pc::ANDEAN_WEDGE_BASE_WIDTH_KM), scaled(pc::ANDEAN_FOREARC_BASE_WIDTH_KM), scaled(pc::ANDEAN_ARC_BASE_WIDTH_KM), scaled(pc::RETROARC_FORELAND_BASE_WIDTH_KM));
        let w4 = if rate >= pc::ANDEAN_PLATEAU_MIN_RATE_CM_PER_YEAR { scaled(pc::ANDEAN_PLATEAU_BASE_WIDTH_KM) } else { 0.0 };
        let arc_end = w1 + w2 + w3;
        let total = arc_end + 1.3 * w4 + w5;
        let (owner, crust) = (self.owner, self.crust);
        let on = |plate: u32| move |c: Cell| owner[c as usize] == plate;
        let edge = |plate: u32| -> Vec<Cell> { s.cells.iter().copied().filter(|&c| owner[c as usize] == plate).collect() };
        let continental = |c: Cell| crust[c as usize] == CrustType::Continental;
        let m = self.reach(&edge(over), total, on(over));
        // Slab age where it enters the trench, carried inland from the nearest trench cell.
        let slab: Vec<Cell> = edge(sub).into_iter().filter(|&c| !continental(c)).collect();
        let near = self.nearest(&slab, total + self.r * self.grid.spacing_rad(), on(over));
        let slab_age = |c: Cell| near.get(&c).zip(self.age).map(|(&(_, src), a)| f64::from(a[src as usize]));
        let old = |c: Cell| w4 > 0.0 && slab_age(c).is_some_and(|t| (pc::ANDEAN_PLATEAU_MIN_SLAB_AGE_MA..=pc::ANDEAN_PLATEAU_MAX_SLAB_AGE_MA).contains(&t));
        let (wn, grid) = (&self.width_noise, self.grid);
        let swell = |c: Cell| (1.0 + 1.6 * wn.fbm(grid.center(c), 6.0, 3)).clamp(0.7, 1.3);
        let back = |c: Cell| if old(c) { arc_end + w4 * swell(c) } else { arc_end };
        let band = |lo: f64, hi: f64| Self::cells_of(&m, |_, d| d > lo && d <= hi);
        let wedge = Self::cells_of(&m, |_, d| d <= w1);
        let fore = band(w1, w1 + w2);
        let arc = band(w1 + w2, arc_end);
        let plateau = Self::cells_of(&m, |c, d| d > arc_end && d <= back(c) && continental(c));
        let foreland = Self::cells_of(&m, |c, d| d > back(c) && d <= back(c) + w5 && continental(c));
        self.push(foreland, pc::foreland_basin(rate, w5), Some(i), age * 0.5);
        self.push(plateau, pc::andean_plateau(rate, w4), Some(i), age);
        self.push(wedge, pc::accretionary_wedge(rate, w1), Some(i), age);
        self.push(fore, pc::forearc_basin(w2), Some(i), age);
        self.push(arc, pc::volcanic_arc(rate, s.length_km), Some(i), age);
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

    /// Rifts open between continents, but do not cut across belts already built at a
    /// convergent margin: a short divergent stretch on an arc microplate's edge (PB2002's
    /// Altiplano) is not the East African rift.
    fn continental_rifts(&mut self, segments: &[BoundarySegment]) {
        use GeologicProvince as G;
        let active: std::collections::HashSet<Cell> = self
            .regions
            .iter()
            .filter(|(_, p)| matches!(p.kind, G::CollisionOrogen | G::OrogenicPlateau | G::ForelandBasin | G::VolcanicArc | G::ForearcBasin | G::AccretionaryWedge | G::OceanTrench))
            .flat_map(|(cells, _)| cells.iter().copied())
            .collect();
        for (i, s) in segments.iter().enumerate() {
            if s.kind != BoundaryKind::Divergent {
                continue;
            }
            let (a, b) = self.sides(s);
            if a.continental && b.continental {
                let half = pc::RIFT_MIN_HALF_WIDTH_KM + self.rng.r#gen::<f64>() * pc::RIFT_HALF_WIDTH_RANGE_KM;
                let cells = Self::cells_of(&self.reach(&s.cells, half, |_| true), |c, _| !active.contains(&c));
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
            self.events.push(Event {
                t_ma: age.round(),
                kind: EventKind::LargeIgneousProvince,
                at: Place::Point { at: p0 },
                magnitude: a,
                causes: vec![],
                note: format!("{} erupts", kind_name(kind)),
            });
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
            self.events.push(Event {
                t_ma: (t_len * 10.0).round() / 10.0,
                kind: EventKind::HotspotTrack,
                at: Place::Point { at: ph },
                magnitude: length,
                causes: vec![],
                note: "hotspot track begins; active end is the youngest".into(),
            });
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

    let mut b = B {
        grid,
        owner,
        crust,
        age,
        r: radius_km,
        plates,
        info,
        edge_dist_km,
        regions: vec![],
        events: vec![],
        rng: StdRng::seed_from_u64(seed),
        width_noise: crate::noise::TerrainNoise::new(seed ^ 0x0B1D_7A11),
        cratons: BTreeMap::new(),
    };

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
