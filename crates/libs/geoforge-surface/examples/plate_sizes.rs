//! Plate size distribution for a range of charge weightings.
//! `cargo run --profile iter -p geoforge-surface --example plate_sizes -- 0 0.3 0.5 0.7`
use geoforge_core::grid::{EquirectGrid, Grid};
use geoforge_surface::tectonics::{partition::RefineParams, weighted_plate_partition};
use std::collections::BTreeMap;

fn main() {
    let alphas: Vec<f64> = std::env::args().skip(1).map(|a| a.parse().unwrap()).collect();
    let grid = EquirectGrid::new(512, 256).unwrap();
    let area: Vec<f64> = (0..grid.len() as u32).map(|c| grid.area_sr(c)).collect();
    let total: f64 = area.iter().sum();
    for a in alphas {
        let (mut l, mut t7, mut n1, mut tiny) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        for seed in [1u64, 2, 3, 7, 42, 99, 1234, 2837, 5555, 97243067] {
            let (owner, _, _) = weighted_plate_partition(&grid, 20, seed, &RefineParams::default(), a);
            let mut m: BTreeMap<u32, f64> = BTreeMap::new();
            for (c, o) in owner.iter().enumerate() {
                *m.entry(*o).or_default() += area[c] / total;
            }
            let mut f: Vec<f64> = m.into_values().collect();
            f.sort_by(|a, b| b.total_cmp(a));
            l.push(f[0]);
            t7.push(f.iter().take(7).sum::<f64>());
            n1.push(f.iter().filter(|x| **x >= 0.01).count() as f64);
            tiny.push(f.iter().filter(|x| **x < 0.002).count() as f64);
        }
        let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len() as f64;
        let mx = |v: &[f64]| v.iter().copied().fold(0.0, f64::max);
        let mn = |v: &[f64]| v.iter().copied().fold(1.0, f64::min);
        println!("alpha {a:4.2}: largest {:.3} ({:.3}-{:.3})  top7 {:.3}  plates>=1% {:.1}  plates<0.2% {:.1}", mean(&l), mn(&l), mx(&l), mean(&t7), mean(&n1), mean(&tiny));
    }
}
