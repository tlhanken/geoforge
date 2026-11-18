//! Tectonics layer: plates, crust type, crust age and the inferred tectonic timeline.
//!
//! Plate geometry is a faithful port of geoforge v1 Stages 1.1–1.3 (same seed, same plates on
//! the same grid; see `tests/v1_parity.rs`). Motion, boundaries, crust and areas are
//! corrected and extended as described in each module.

pub mod boundaries;
pub mod crust;
pub mod electrostatic;
pub mod motion;
pub mod partition;

use crate::tectonics::boundaries::BoundaryParams;
use crate::tectonics::crust::CrustMode;
use crate::tectonics::partition::{IslandStats, RefineParams};
use geoforge_core::grid::Grid;
use geoforge_core::key::Key;
use geoforge_core::product::{Generator, LayerId, Precondition, ProductData, ProductId, Products, RunContext};
use geoforge_core::shapes::{Field, Partition, Table};
use geoforge_core::surface::{CrustType, Plate};
use geoforge_core::{BodyContext, Error, UnitVec3};
use rand::rngs::StdRng;
use rand::SeedableRng;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct TectonicsParams {
    pub plate_count: usize,
    pub refine: RefineParams,
    pub min_speed_cm_yr: f64,
    pub max_speed_cm_yr: f64,
    pub boundaries: BoundaryParams,
    pub crust: CrustMode,
}

impl Default for TectonicsParams {
    fn default() -> Self {
        Self {
            plate_count: 20,
            refine: RefineParams::default(),
            min_speed_cm_yr: 1.0,
            max_speed_cm_yr: 10.0,
            boundaries: BoundaryParams::default(),
            crust: CrustMode::default(),
        }
    }
}

/// Plate ids per cell and seed points, exactly as v1 Stages 1.1–1.3 produce them for `seed`.
#[must_use]
pub fn v1_plate_partition(grid: &dyn Grid, plate_count: usize, seed: u64, refine: &RefineParams) -> (Vec<u32>, Vec<UnitVec3>, IslandStats) {
    let cfg = electrostatic::ElectrostaticConfig::default();
    let mut rng = StdRng::seed_from_u64(seed);
    let mut charges = electrostatic::random_charges(plate_count, &mut rng, &cfg);
    electrostatic::simulate_equilibrium(&mut charges, &cfg);
    let seeds = electrostatic::seed_points(&charges);
    let owner = partition::voronoi(grid, &seeds);
    let mut owner = partition::domain_warp(grid, &owner, seed, refine);
    if refine.smoothing_iterations > 0 {
        partition::smooth(grid, &mut owner, refine.smoothing_iterations);
    }
    let stats = partition::remove_islands(grid, &mut owner);
    (owner, seeds, stats)
}

/// Generate the plates partition (geometry + Euler poles) for a seed.
#[must_use]
pub fn generate_plates(grid: &dyn Grid, seed: u64, radius_km: f64, p: &TectonicsParams) -> Partition<Plate> {
    let (owner, seeds, _) = v1_plate_partition(grid, p.plate_count, seed, &p.refine);
    let motion = motion::v1_random_motion(seeds.len(), seed, p.min_speed_cm_yr, p.max_speed_cm_yr);
    let mut table = Table::new();
    for (i, (pt, (dir, speed))) in seeds.iter().zip(motion).enumerate() {
        let (axis, rate) = motion::euler_pole_for(*pt, dir, speed, radius_km);
        table.insert(
            Key(i as u64 + 1),
            Plate { name: format!("Plate {}", i + 1), euler_pole: axis, rate_deg_per_myr: rate, seed_point: *pt },
        );
    }
    let mut plates = Partition { owner: Field::from_values(grid.spec(), owner), table };
    plates.prune();
    plates
}

pub struct Tectonics;

impl Generator for Tectonics {
    type Params = TectonicsParams;
    const LAYER: LayerId = LayerId::Tectonics;

    fn needs() -> &'static [ProductId] {
        &[]
    }

    fn applies(body: &BodyContext) -> Result<(), Precondition> {
        if body.interior.has_lithosphere() {
            Ok(())
        } else {
            Err(Precondition("tectonics needs a rocky crust (has_lithosphere)".into()))
        }
    }

    fn run(ctx: &RunContext<'_>, inputs: &Products, p: &TectonicsParams) -> Result<Products, Error> {
        let grid = ctx.grid;
        let r = ctx.radius_km();
        let mut out = Products::new();

        let plates = match inputs.plates() {
            Some(pl) => pl.clone(),
            None => {
                let pl = generate_plates(grid, ctx.seed(ProductId::Plates).value(), r, p);
                out.insert(ProductData::Plates(pl.clone()));
                pl
            }
        };
        plates.validate()?;
        let owner = &plates.owner.values;

        let crust: Vec<CrustType> = match inputs.crust_type() {
            Some(f) => f.values.clone(),
            None => {
                let areas = partition::areas_km2(grid, owner, r);
                let kinds = crust::plate_kinds(&areas, p.crust, ctx.seed(ProductId::CrustType));
                let v: Vec<CrustType> = owner.iter().map(|o| kinds.get(o).copied().unwrap_or(CrustType::Oceanic)).collect();
                out.insert(ProductData::CrustType(Field::from_values(grid.spec(), v.clone())));
                v
            }
        };

        let segments = boundaries::analyze(grid, &plates, r, &p.boundaries);
        let (age, age_src) = crust::crust_age(grid, owner, &crust, &segments, r);
        let age = match inputs.crust_age() {
            Some(f) => f.values.clone(),
            None => {
                out.insert(ProductData::CrustAge(Field::from_values(grid.spec(), age.clone())));
                age
            }
        };
        if inputs.tectonic_events().is_none() {
            let tl = crust::infer_events(grid, owner, &crust, &age, &age_src, &segments);
            out.insert(ProductData::TectonicEvents(tl));
        }
        Ok(out)
    }
}
