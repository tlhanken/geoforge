//! Geologic provinces of a rigid lid.
//!
//! The base is the crust's split: cratered highlands and volcanic plains (tesserae instead of
//! highlands on an episodic lid, cratered ice on a frozen body). On top, in order of age: impact basins (the largest
//! craters, more on older ground), volcanic rises with shield volcanoes, rifts radiating from
//! the largest rise, and, on an episodic lid, coronae. Each feature is a province row with
//! its centre (`anchor`), size (`width_km`) and strength (`intensity`), so the elevation
//! layer can draw it and an editor can move or resize it.

use super::landforms::{craters, BASIN_MIN_KM};
use crate::geology::lithology::{self, RockDef};
use crate::preconditions;
use geoforge_core::body::{Species, TectonicRegime, EARTH_MASS_KG};
use geoforge_core::grid::{cells_within, components, Cell, Grid};
use geoforge_core::key::Key;
use geoforge_core::product::{Generator, LayerId, Precondition, ProductData, ProductId, Products, RunContext};
use geoforge_core::shapes::{Event, EventKind, Field, Partition, Place, Table, Timeline};
use geoforge_core::surface::{CrustType, GeologicProvince as G, Province};
use geoforge_core::{BodyContext, Error, UnitVec3};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde::{Deserialize, Serialize};

/// Counts `None` = chosen from the body's mass and regime.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct LidGeologyParams {
    /// Broad volcanic rises (Tharsis, Elysium; Beta and Atla Regio).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub volcanic_rises: Option<u32>,
    pub shields_per_rise: u32,
    /// Rift valleys (Valles Marineris; Venus's chasmata).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub rifts: Option<u32>,
    /// Coronae (episodic lids only by default).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub coronae: Option<u32>,
    /// Rocks added to the built-in catalog (see [`RockDef`]).
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub extra_rocks: Vec<RockDef>,
}

impl Default for LidGeologyParams {
    fn default() -> Self {
        Self { volcanic_rises: None, shields_per_rise: 3, rifts: None, coronae: None, extra_rocks: Vec::new() }
    }
}

/// Feature counts for a body.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Counts {
    rises: u32,
    rifts: u32,
    coronae: u32,
}

impl LidGeologyParams {
    fn counts(&self, body: &BodyContext) -> Counts {
        // Long-lived mantle plumes need a big body: Mars has two rises; Mercury and the Moon
        // cooled too fast and have none.
        let mass = body.exterior.mass_kg / EARTH_MASS_KG;
        let episodic = body.interior.regime == TectonicRegime::Episodic;
        let rises = if episodic {
            (10.0 * mass).round().min(12.0)
        } else if mass < 0.07 {
            0.0
        } else {
            (20.0 * mass).round().min(4.0)
        } as u32;
        let rises = self.volcanic_rises.unwrap_or(rises);
        let rifts = self.rifts.unwrap_or(if episodic { rises.min(3) } else { rises.min(1) });
        let area_venus = 4.6e8;
        let coronae = self.coronae.unwrap_or(if episodic { (40.0 * body.surface_area_km2() / area_venus).round() as u32 } else { 0 });
        Counts { rises, rifts, coronae }
    }
}

pub struct LidGeology;

/// Provinces being painted, later ones over earlier ones.
struct Painter<'a> {
    grid: &'a dyn Grid,
    radius_km: f64,
    owner: Vec<u32>,
    table: Table<Province>,
    events: Table<Event>,
}

impl Painter<'_> {
    fn province(kind: G, width_km: f64, intensity: f64, age_ma: f64, anchor: Option<UnitVec3>) -> Province {
        Province { kind, elevation_intensity: 0.0, roughness: 0.0, width_km, intensity, convergence_rate_cm_yr: 0.0, boundary: None, age_ma, anchor }
    }

    fn add(&mut self, prov: Province, cells: &[Cell], event: Option<EventKind>) -> Key {
        let key = Key(self.table.rows.len() as u64 + 1);
        for &c in cells {
            self.owner[c as usize] = key.0 as u32;
        }
        if let Some(kind) = event {
            let t_ma = prov.age_ma;
            let magnitude = prov.width_km;
            self.events.insert(Key(self.events.rows.len() as u64 + 1), Event { t_ma, kind, at: Place::Entity { product: ProductId::Provinces, key }, magnitude, causes: Vec::new(), note: String::new() });
        }
        self.table.insert(key, prov);
        key
    }

    /// Cells within `radius_km` of `at`.
    fn disk(&self, at: UnitVec3, radius_km: f64) -> Vec<Cell> {
        let r = (radius_km / self.radius_km).max(self.grid.spacing_rad() * 0.5);
        cells_within(self.grid, &[self.grid.locate(at)], r, |_| true)
    }

    /// A circular feature.
    fn circle(&mut self, kind: G, at: UnitVec3, radius_km: f64, intensity: f64, age_ma: f64, event: Option<EventKind>) {
        let cells = self.disk(at, radius_km);
        self.add(Self::province(kind, 2.0 * radius_km, intensity, age_ma, Some(at)), &cells, event);
    }
}

/// A random point at least `min_sep_rad` from every point in `taken` (best of a few tries).
fn spaced_point(rng: &mut StdRng, taken: &[(UnitVec3, f64)], radius_rad: f64) -> UnitVec3 {
    let mut best = (UnitVec3::from_lat_lon(0.0, 0.0), f64::NEG_INFINITY);
    for _ in 0..40 {
        let p = UnitVec3::from_lat_lon(rng.gen::<f64>().mul_add(2.0, -1.0).asin().to_degrees(), rng.gen_range(-180.0..180.0));
        let gap = taken.iter().map(|(q, r)| p.angle_to(*q) - r - radius_rad).fold(f64::INFINITY, f64::min);
        if gap > 0.0 {
            return p;
        }
        if gap > best.1 {
            best = (p, gap);
        }
    }
    best.0
}

impl Generator for LidGeology {
    type Params = LidGeologyParams;
    const NAME: &'static str = "lid_geology";
    const VERSION: u32 = 1;
    const LAYER: LayerId = LayerId::Geology;

    fn needs() -> &'static [ProductId] {
        &[ProductId::CrustType, ProductId::CrustAge]
    }
    fn applies(body: &BodyContext) -> Result<(), Precondition> {
        preconditions::rigid_lid(body)
    }

    fn run(ctx: &RunContext<'_>, inputs: &Products, p: &LidGeologyParams) -> Result<Products, Error> {
        let grid = ctx.grid;
        let crust = &inputs.crust_type().ok_or(Error::MissingInput(ProductId::CrustType))?.values;
        let age = &inputs.crust_age().ok_or(Error::MissingInput(ProductId::CrustAge))?.values;
        let mut out = Products::new();

        let provinces = match inputs.provinces() {
            Some(pv) => pv.clone(),
            None => {
                let (pv, events) = build(ctx, crust, age, p);
                if inputs.geology_events().is_none() {
                    out.insert(ProductData::GeologyEvents(events));
                }
                out.insert(ProductData::Provinces(pv.clone()));
                pv
            }
        };
        if inputs.lithology().is_none() || inputs.rock_age().is_none() {
            let catalog = lithology::Catalog::new(&p.extra_rocks)?;
            let (lith, rock_age) = lithology::build(grid, &provinces, crust, Some(age), ctx.seed(ProductId::Lithology).value(), &catalog);
            if inputs.lithology().is_none() {
                out.insert(ProductData::Lithology(lith));
            }
            if inputs.rock_age().is_none() {
                out.insert(ProductData::RockAge(Field::from_values(grid.spec(), rock_age)));
            }
        }
        Ok(out)
    }
}

/// Paint the provinces and record their events.
fn build(ctx: &RunContext<'_>, crust: &[CrustType], age: &[f32], p: &LidGeologyParams) -> (Partition<Province>, Timeline) {
    let grid = ctx.grid;
    let body = ctx.body;
    let r = ctx.radius_km();
    let episodic = body.interior.regime == TectonicRegime::Episodic;
    let icy = body.interior.surface_shell().is_some_and(|s| s.composition.fraction(Species::Ice) > 0.5);
    let counts = p.counts(body);
    let mut rng = StdRng::seed_from_u64(ctx.seed(ProductId::Provinces).value());
    let mut paint = Painter { grid, radius_km: r, owner: vec![0; grid.len()], table: Table::new(), events: Table::new() };
    let mean_age = |cells: &[Cell]| cells.iter().map(|&c| f64::from(age[c as usize])).sum::<f64>() / cells.len().max(1) as f64;

    // 1. The crust's own regions.
    let all: Vec<Cell> = grid.cells().collect();
    for high in [true, false] {
        let member = |c: Cell| (crust[c as usize] == CrustType::Continental) == high;
        for region in components(grid, &all, member) {
            let kind = match (high, episodic, icy) {
                (true, _, true) => G::CrateredIce,
                (true, true, false) => G::Tessera,
                (true, false, false) => G::CrateredHighlands,
                (false, ..) => G::VolcanicPlains,
            };
            let a = mean_age(&region);
            let width = (region.iter().map(|&c| grid.area_sr(c)).sum::<f64>()).sqrt() * r;
            paint.add(Painter::province(kind, width, 0.5, a, None), &region, None);
        }
    }

    // 2. Impact basins, oldest (largest) first.
    // The largest known (the Moon's South Pole–Aitken) is about 1.5 body radii across.
    let basins = craters(grid, age, r, BASIN_MIN_KM, 1.5 * r, ctx.seed(ProductId::Provinces).child(b"basins").value());
    for b in &basins {
        let a = f64::from(age[grid.locate(b.at) as usize]) * rng.gen_range(0.9..1.0);
        paint.circle(G::ImpactBasin, b.at, b.diameter_km / 2.0, 1.0, a, Some(EventKind::ImpactBasin));
    }

    // 3. Volcanic rises, each with shield volcanoes; the first is the largest.
    let mut rises: Vec<(UnitVec3, f64)> = Vec::new();
    for i in 0..counts.rises {
        let radius_km = if episodic {
            r * 0.12 * rng.gen_range(0.6..1.4)
        } else if i == 0 {
            r * 0.45 * rng.gen_range(0.8..1.1)
        } else {
            r * 0.2 * rng.gen_range(0.7..1.3)
        };
        let at = spaced_point(&mut rng, &rises, radius_km / r);
        let a = f64::from(age[grid.locate(at) as usize]).min(3700.0) * rng.gen_range(0.3..0.9);
        paint.circle(G::VolcanicRise, at, radius_km, rng.gen_range(0.6..1.0), a, Some(EventKind::VolcanicRise));
        rises.push((at, radius_km / r));
    }

    // 4. Rifts: from a rise's flank outward along a great circle.
    for i in 0..counts.rifts as usize {
        let Some(&(centre, rise_rad)) = rises.get(i % rises.len().max(1)) else { break };
        let bearing = rng.gen_range(0.0..360.0);
        let start = centre.travel(bearing, rise_rad * 0.5);
        let length = rng.gen_range(0.6..1.0); // radians: 0.6–1.0 body radii
        let half_km = rng.gen_range(80.0..130.0);
        let steps = (length / grid.spacing_rad()).ceil().max(1.0) as usize;
        let spine: Vec<Cell> = (0..=steps).map(|k| grid.locate(start.travel(bearing, length * k as f64 / steps as f64))).collect();
        let cells = cells_within(grid, &spine, half_km / r, |_| true);
        let a = rng.gen_range(1500.0..3500.0_f64).min(f64::from(age[spine[0] as usize]));
        paint.add(Painter::province(G::ContinentalRift, 2.0 * half_km, rng.gen_range(0.6..1.0), a, Some(start)), &cells, None);
    }

    // 5. Coronae.
    for _ in 0..counts.coronae {
        let radius_km = rng.gen_range(80.0..300.0);
        let at = spaced_point(&mut rng, &[], 0.0);
        let a = f64::from(age[grid.locate(at) as usize]) * rng.gen_range(0.2..0.9);
        paint.circle(G::Corona, at, radius_km, rng.gen_range(0.5..1.0), a, None);
    }

    // 6. Shields: on their rises, plus a few scattered.
    let scattered = counts.rises / 2;
    for i in 0..counts.rises * p.shields_per_rise + scattered {
        let (at, rise_age) = match rises.get((i / p.shields_per_rise.max(1)) as usize) {
            Some(&(c, rr)) if i < counts.rises * p.shields_per_rise => (c.travel(rng.gen_range(0.0..360.0), rr * 0.7 * rng.gen::<f64>().sqrt()), 3000.0),
            _ => (spaced_point(&mut rng, &[], 0.0), 2000.0),
        };
        let radius_km = rng.gen_range(80.0..300.0) * if episodic { 0.7 } else { 1.0 };
        let a = rng.gen_range(0.05..0.6) * rise_age;
        paint.circle(G::ShieldVolcano, at, radius_km, rng.gen_range(0.4..1.0), a, None);
    }

    let mut pv = Partition { owner: Field::from_values(grid.spec(), paint.owner), table: paint.table };
    pv.prune();
    (pv, Timeline { events: paint.events })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lid::tectonics::{LidTectonics, LidTectonicsParams};
    use geoforge_core::grid::GridSpec;
    use geoforge_core::{Archetype, Seed};

    fn run(a: Archetype) -> Products {
        let body = a.context(Key(1), a.slug());
        let grid = GridSpec::Equirect { width: 360, height: 180 }.build().unwrap();
        let ctx = RunContext { root_seed: Seed(5), body: &body, grid: &grid };
        let mut all = LidTectonics::run(&ctx, &Products::new(), &LidTectonicsParams::default()).unwrap();
        all.extend(LidGeology::run(&ctx, &all, &LidGeologyParams::default()).unwrap());
        all
    }

    fn kinds(pr: &Products) -> Vec<G> {
        pr.provinces().unwrap().table.rows.values().map(|p| p.kind).collect()
    }

    #[test]
    fn mars_moon_and_venus_differ() {
        let mars = run(Archetype::Mars);
        let k = kinds(&mars);
        assert!(k.contains(&G::CrateredHighlands) && k.contains(&G::VolcanicPlains));
        assert_eq!(k.iter().filter(|x| **x == G::VolcanicRise).count(), 2, "Tharsis and Elysium");
        assert!(k.contains(&G::ShieldVolcano) && k.contains(&G::ContinentalRift));
        assert!(k.contains(&G::ImpactBasin), "old ground keeps its basins");

        let moon = kinds(&run(Archetype::Moon));
        assert!(!moon.contains(&G::VolcanicRise) && !moon.contains(&G::ShieldVolcano));
        assert!(moon.iter().filter(|x| **x == G::ImpactBasin).count() >= 2);

        let venus = kinds(&run(Archetype::Venus));
        assert!(venus.contains(&G::Tessera) && venus.contains(&G::Corona));
        assert!(!venus.contains(&G::ImpactBasin) && !venus.contains(&G::CrateredHighlands), "young surface");

        // Lithology and rock ages come along, from the built-in catalog.
        let lith = mars.lithology().unwrap();
        lith.validate().unwrap();
        assert!(lith.owner.values.iter().all(|&o| o != 0));
        assert!(mars.rock_age().unwrap().values.iter().all(|a| (0.0..=4400.0).contains(a)));
    }
}
