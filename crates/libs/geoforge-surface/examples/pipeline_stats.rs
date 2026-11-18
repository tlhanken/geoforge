//! Run tectonics → geology → elevation for a few seeds; print realism stats and dump fields.
use geoforge_core::grid::{EquirectGrid, Grid};
use geoforge_core::product::{Generator, Products, RunContext};
use geoforge_core::surface::{Category, CrustType};
use geoforge_core::{BodyContext, Key, Seed};
use geoforge_surface::{Elevation, Geology, Tectonics};
use std::io::Write;

fn main() {
    let w: u32 = std::env::var("W").ok().and_then(|v| v.parse().ok()).unwrap_or(1800);
    let grid = EquirectGrid::new(w, w / 2).unwrap();
    let body = BodyContext::earth_like(Key(1), "Test");
    let seeds: Vec<u64> = std::env::args().skip(1).filter_map(|a| a.parse().ok()).collect();
    for seed in seeds {
        let ctx = RunContext { root_seed: Seed(seed), body: &body, grid: &grid };
        let mut all = Products::new();
        let t = std::time::Instant::now();
        all.extend(Tectonics::run(&ctx, &all, &Default::default()).unwrap());
        let t1 = t.elapsed();
        all.extend(Geology::run(&ctx, &all, &Default::default()).unwrap());
        let t2 = t.elapsed();
        all.extend(Elevation::run(&ctx, &all, &Default::default()).unwrap());
        let t3 = t.elapsed();
        let e = &all.elevation().unwrap().values;
        let crust = &all.crust_type().unwrap().values;
        let prov = all.provinces().unwrap();
        let area: Vec<f64> = grid.cells().map(|c| grid.area_sr(c)).collect();
        let tot: f64 = area.iter().sum();
        let land: f64 = grid.cells().filter(|&c| e[c as usize] >= 0.0).map(|c| area[c as usize]).sum::<f64>() / tot;
        let cont: f64 = grid.cells().filter(|&c| crust[c as usize] == CrustType::Continental).map(|c| area[c as usize]).sum::<f64>() / tot;
        let (mut od, mut ow, mut le, mut lw) = (0.0, 0.0, 0.0, 0.0);
        for c in grid.cells() { let i = c as usize; let z = f64::from(e[i]); if z < 0.0 { od += -z * area[i]; ow += area[i]; } else { le += z * area[i]; lw += area[i]; } }
        let shelf: f64 = grid.cells().filter(|&c| { let z = e[c as usize]; z < 0.0 && z > -200.0 }).map(|c| area[c as usize]).sum::<f64>() / tot;
        let mx = e.iter().cloned().fold(f32::MIN, f32::max); let mn = e.iter().cloned().fold(f32::MAX, f32::min);
        println!("seed {seed}: [tect {t1:.1?}, geo {:.1?}, elev {:.1?}] provinces {}", t2 - t1, t3 - t2, prov.table.len());
        println!("   land {:.1}% (cont crust {:.1}%), shelf(0..-200m) {:.1}%, mean ocean depth {:.0} m, mean land elev {:.0} m, max {:.0} m, min {:.0} m", land * 100.0, cont * 100.0, shelf * 100.0, od / ow, le / lw, mx, mn);
        // Hypsometry, 1000 m bins.
        let mut bins = [0.0f64; 21];
        for c in grid.cells() { let b = ((f64::from(e[c as usize]) + 11000.0) / 1000.0).floor().clamp(0.0, 20.0) as usize; bins[b] += area[c as usize] / tot * 100.0; }
        print!("   hypsometry %/km from -11 km: ");
        for b in bins { print!("{b:.1} "); }
        println!();
        let mut counts = std::collections::BTreeMap::new();
        for (_, p) in prov.table.iter() { *counts.entry(p.kind.name()).or_insert(0) += 1; }
        println!("   {counts:?}");
        println!("   events: tectonic {}, geology {}", all.tectonic_events().unwrap().events.len(), all.geology_events().unwrap().events.len());
        if let Ok(d) = std::env::var("DUMP") {
            let mut f = std::fs::File::create(format!("{d}/p{seed}.bin")).unwrap();
            for c in grid.cells() { let i = c as usize; let k = prov.owner_of(c).and_then(|k| prov.table.get(k)).map_or(255u8, |p| p.kind.to_u8()); f.write_all(&e[i].to_le_bytes()).unwrap(); f.write_all(&[k, crust[i] as u8]).unwrap(); }
        }
    }
}
