//! The ported plate partition reproduces geoforge v1 (Stages 1.1–1.3) for the same seed on
//! the same 1800 × 900 grid. Differences are limited to v1's seam-column and pole-row
//! smoothing skips, which the grid-generic port does not replicate.

use geoforge::WorldMap;
use geoforge_core::grid::EquirectGrid;
use geoforge_surface::tectonics::{partition::RefineParams, v1_plate_partition};

fn v1_plates(seed: u64) -> Vec<u16> {
    let mut world = WorldMap::new(1800, 900, seed).unwrap();
    world.tectonics().generate_plates(20).unwrap();
    world.tectonics().roughen_boundaries(None).unwrap();
    world.tectonics().deisland(None).unwrap();
    world.tectonics.as_ref().unwrap().data.clone()
}

#[test]
fn plate_partition_matches_v1() {
    let grid = EquirectGrid::new(1800, 900).unwrap();
    for seed in [2837u64, 97_243_067, 42] {
        let v1 = v1_plates(seed);
        let (v2, _, _) = v1_plate_partition(&grid, 20, seed, &RefineParams::default());
        let same = v1.iter().zip(&v2).filter(|(a, b)| u32::from(**a) == **b).count();
        let frac = same as f64 / v1.len() as f64;
        eprintln!("seed {seed}: {:.5}% of cells identical ({} differ)", frac * 100.0, v1.len() - same);
        assert!(frac > 0.999, "seed {seed}: only {:.4}% identical", frac * 100.0);
    }
}
