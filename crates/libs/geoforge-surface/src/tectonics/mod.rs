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
    /// How strongly a seed's electrostatic charge sets its plate's size: cells go to the seed
    /// minimising `distance / charge^size_weighting`. 0 reproduces v1 exactly (charge only
    /// affects where seeds settle); the default is tuned so plate sizes follow Earth's.
    pub size_weighting: f64,
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
            size_weighting: SIZE_WEIGHTING,
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

/// Default charge weighting (see [`TectonicsParams::size_weighting`]).
pub const SIZE_WEIGHTING: f64 = 0.22;

/// v1's plate pipeline with the partition weighted by charge; `alpha = 0` is
/// [`v1_plate_partition`].
#[must_use]
pub fn weighted_plate_partition(grid: &dyn Grid, plate_count: usize, seed: u64, refine: &RefineParams, alpha: f64) -> (Vec<u32>, Vec<UnitVec3>, IslandStats) {
    if alpha == 0.0 {
        return v1_plate_partition(grid, plate_count, seed, refine);
    }
    let cfg = electrostatic::ElectrostaticConfig::default();
    let mut rng = StdRng::seed_from_u64(seed);
    let mut charges = electrostatic::random_charges(plate_count, &mut rng, &cfg);
    electrostatic::simulate_equilibrium(&mut charges, &cfg);
    let seeds = electrostatic::seed_points(&charges);
    let q: Vec<f64> = charges.iter().map(|c| c.charge).collect();
    let owner = partition::weighted_voronoi(grid, &seeds, &q, alpha);
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
    let (owner, seeds, _) = weighted_plate_partition(grid, p.plate_count, seed, &p.refine, p.size_weighting);
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
    const NAME: &'static str = "tectonics";
    const VERSION: u32 = 2;
    const LAYER: LayerId = LayerId::Tectonics;

    fn needs() -> &'static [ProductId] {
        &[]
    }

    fn applies(body: &BodyContext) -> Result<(), Precondition> {
        crate::preconditions::plate_tectonics(body)
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
                let v = crust::assign(grid, owner, p.crust, ctx.seed(ProductId::CrustType), r);
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
