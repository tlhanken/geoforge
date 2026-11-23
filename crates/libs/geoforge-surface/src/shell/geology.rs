//! Geologic provinces of an ice shell.
//!
//! The base is the shell's split into old plains and renewed ice: on a thin shell ridged
//! plains and chaos terrain, on a thick one dark cratered terrain and bright grooved terrain.
//! On top, in order of age: impact basins (palimpsests on old thick shells), then on a thin
//! shell its bands (dilational, as wide as 60 km, hundreds of kilometres long) and its ridge
//! lineae (quasi-linear double ridges running for thousands of kilometres, the youngest
//! features; Leonard, Patthoff & Senske 2024). Each feature is a province row with its start
//! (`anchor`), width and strength, so the elevation layer can draw it and an editor can move
//! it.

use crate::geology::lithology::{self, RockDef};
use crate::lid::geology::Painter;
use crate::lid::landforms::{craters, BASIN_MIN_KM};
use crate::preconditions;
use crate::shell::tectonics::{IceShellTectonicsParams, ShellStyle};
use geoforge_core::grid::{cells_within, components, Cell, Grid};
use geoforge_core::product::{Generator, LayerId, Precondition, ProductData, ProductId, Products, RunContext};
use geoforge_core::shapes::{EventKind, Field, Partition, Timeline};
use geoforge_core::surface::{CrustType, GeologicProvince as G, Province};
use geoforge_core::{BodyContext, Error, UnitVec3};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde::{Deserialize, Serialize};

/// `None` = chosen from the shell's style.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct IceShellGeologyParams {
    /// Share of the surface in bands.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub band_share: Option<f64>,
    /// Ridge lineae drawn.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub lineae: Option<u32>,
    /// Rocks added to the built-in catalog (see [`RockDef`]).
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub extra_rocks: Vec<RockDef>,
}

pub struct IceShellGeology;

impl Generator for IceShellGeology {
    type Params = IceShellGeologyParams;
    const NAME: &'static str = "ice_shell_geology";
    const VERSION: u32 = 1;
    const LAYER: LayerId = LayerId::Geology;

    fn needs() -> &'static [ProductId] {
        &[ProductId::CrustType, ProductId::CrustAge]
    }
    fn applies(body: &BodyContext) -> Result<(), Precondition> {
        preconditions::ice_shell(body)
    }

    fn run(ctx: &RunContext<'_>, inputs: &Products, p: &IceShellGeologyParams) -> Result<Products, Error> {
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

/// A point anywhere on the sphere, evenly.
fn anywhere(rng: &mut StdRng) -> UnitVec3 {
    UnitVec3::from_lat_lon(rng.gen::<f64>().mul_add(2.0, -1.0).asin().to_degrees(), rng.gen_range(-180.0..180.0))
}

/// Cells within `half_km` of a great-circle arc from `start` along `bearing` for `length_km`.
fn arc(grid: &dyn Grid, radius_km: f64, start: UnitVec3, bearing: f64, length_km: f64, half_km: f64) -> Vec<Cell> {
    let length = length_km / radius_km;
    let steps = (length / grid.spacing_rad()).ceil().max(1.0) as usize;
    let spine: Vec<Cell> = (0..=steps).map(|k| grid.locate(start.travel(bearing, length * k as f64 / steps as f64))).collect();
    cells_within(grid, &spine, (half_km / radius_km).max(0.5 * grid.spacing_rad()), |_| true)
}

/// Paint the provinces and record their events.
fn build(ctx: &RunContext<'_>, crust: &[CrustType], age: &[f32], p: &IceShellGeologyParams) -> (Partition<Province>, Timeline) {
    let grid = ctx.grid;
    let r = ctx.radius_km();
    let style = ShellStyle::of(ctx.body, &IceShellTectonicsParams::default());
    let thin = style.thick < 0.5;
    let mut rng = StdRng::seed_from_u64(ctx.seed(ProductId::Provinces).value());
    let mut paint = Painter::new(grid, r);
    let mean_age = |cells: &[Cell]| cells.iter().map(|&c| f64::from(age[c as usize])).sum::<f64>() / cells.len().max(1) as f64;

    // 1. The shell's own regions.
    let all: Vec<Cell> = grid.cells().collect();
    for old in [true, false] {
        let member = |c: Cell| (crust[c as usize] == CrustType::Continental) == old;
        for region in components(grid, &all, member) {
            let kind = match (old, thin) {
                (true, true) => G::RidgedPlains,
                (false, true) => G::ChaosTerrain,
                (true, false) => G::CrateredIce,
                (false, false) => G::GroovedTerrain,
            };
            let width = (region.iter().map(|&c| grid.area_sr(c)).sum::<f64>()).sqrt() * r;
            paint.add(Painter::province(kind, width, 0.5, mean_age(&region), None), &region, None);
        }
    }

    // 2. Impact basins, oldest (largest) first: only old ground keeps them.
    let basins = craters(grid, age, r, BASIN_MIN_KM, 1.5 * r, ctx.seed(ProductId::Provinces).child(b"basins").value());
    for b in &basins {
        let a = f64::from(age[grid.locate(b.at) as usize]) * rng.gen_range(0.9..1.0);
        paint.circle(G::ImpactBasin, b.at, b.diameter_km / 2.0, 1.0, a, Some(EventKind::ImpactBasin));
    }

    // 3. Bands, until they cover their share of the surface.
    let band_share = p.band_share.unwrap_or(if thin { 0.05 * (1.0 - style.thick) } else { 0.0 });
    let (mut covered, total) = (0.0, 4.0 * std::f64::consts::PI);
    while covered < band_share * total {
        let (start, bearing) = (anywhere(&mut rng), rng.gen_range(0.0..360.0));
        let (length_km, half_km) = (rng.gen_range(300.0..1500.0), rng.gen_range(5.0..30.0));
        let cells = arc(grid, r, start, bearing, length_km, half_km);
        covered += cells.iter().map(|&c| grid.area_sr(c)).sum::<f64>();
        let a = style.plains_age_ma * rng.gen_range(0.3..0.9);
        paint.add(Painter::province(G::DilationalBand, 2.0 * half_km, rng.gen_range(0.5..1.0), a, Some(start)), &cells, None);
    }

    // 4. Ridge lineae: the youngest, crossing everything.
    let lineae = p.lineae.unwrap_or((6.0 + 30.0 * (1.0 - style.thick)).round() as u32);
    for _ in 0..lineae {
        let (start, bearing) = (anywhere(&mut rng), rng.gen_range(0.0..360.0));
        let half_km = rng.gen_range(2.0..6.0);
        let cells = arc(grid, r, start, bearing, rng.gen_range(800.0..3000.0), half_km);
        let a = style.plains_age_ma.min(style.renewed_age_ma * 2.0) * rng.gen_range(0.05..0.5);
        paint.add(Painter::province(G::RidgeLinea, 2.0 * half_km, rng.gen_range(0.5..1.0), a, Some(start)), &cells, None);
    }
    paint.finish()
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::shell::tectonics::IceShellTectonics;
    use geoforge_core::grid::GridSpec;
    use geoforge_core::key::Key;
    use geoforge_core::{Archetype, Seed};

    pub(crate) fn run(a: Archetype) -> Products {
        let body = a.context(Key(1), a.slug());
        let grid = GridSpec::Equirect { width: 720, height: 360 }.build().unwrap();
        let ctx = RunContext { root_seed: Seed(5), body: &body, grid: &grid };
        let mut all = IceShellTectonics::run(&ctx, &Products::new(), &IceShellTectonicsParams::default()).unwrap();
        all.extend(IceShellGeology::run(&ctx, &all, &IceShellGeologyParams::default()).unwrap());
        all
    }

    #[test]
    fn europa_and_titan_differ() {
        let europa = run(Archetype::Europa);
        let pv = europa.provinces().unwrap();
        let count = |k: G| pv.table.rows.values().filter(|p| p.kind == k).count();
        assert!(count(G::RidgedPlains) > 0 && count(G::ChaosTerrain) > 10, "chaos in many patches: {}", count(G::ChaosTerrain));
        assert!(count(G::DilationalBand) > 3 && count(G::RidgeLinea) > 20);
        assert_eq!(count(G::ImpactBasin), 0, "a young shell keeps no basins");
        assert!(count(G::GroovedTerrain) == 0 && count(G::CrateredIce) == 0);
        let lith = europa.lithology().unwrap();
        lith.validate().unwrap();
        let titan = run(Archetype::Titan);
        let kinds: Vec<G> = titan.provinces().unwrap().table.rows.values().map(|p| p.kind).collect();
        assert!(kinds.contains(&G::GroovedTerrain) && kinds.contains(&G::CrateredIce) && !kinds.contains(&G::ChaosTerrain));
    }
}
