use geoforge_core::grid::{EquirectGrid, Grid};
use geoforge_core::surface::{BoundaryKind, CrustType};
use geoforge_surface::tectonics::{boundaries, crust, generate_plates, partition, TectonicsParams};
use std::io::Write;

fn main() {
    let grid = EquirectGrid::new(1800, 900).unwrap();
    let r = 6371.0;
    let total_area = 4.0 * std::f64::consts::PI * r * r;
    let dump = std::env::var("DUMP").ok();
    for seed in [2837u64, 97_243_067, 42, 7, 1234] {
        let p = TectonicsParams::default();
        let plates = generate_plates(&grid, seed, r, &p);
        let owner = &plates.owner.values;
        let areas = partition::areas_km2(&grid, owner, r);
        let mut sorted: Vec<f64> = areas.values().map(|a| a / total_area * 100.0).collect();
        sorted.sort_by(|a, b| b.total_cmp(a));
        let v1k = crust::plate_kinds(&areas, crust::CrustMode::PlateSize { continental_percentile: 0.3 }, geoforge_core::Seed(seed));
        let v1_cont: f64 = areas.iter().filter(|(id, _)| v1k[id] == CrustType::Continental).map(|(_, a)| a).sum::<f64>() / total_area;
        let tk = crust::plate_kinds(&areas, crust::CrustMode::default(), geoforge_core::Seed(seed));
        let t_cont: f64 = areas.iter().filter(|(id, _)| tk[id] == CrustType::Continental).map(|(_, a)| a).sum::<f64>() / total_area;
        let t0 = std::time::Instant::now();
        let segs = boundaries::analyze(&grid, &plates, r, &p.boundaries);
        let dt = t0.elapsed();
        let len = |k| segs.iter().filter(|s| s.kind == k).map(|s| s.length_km).sum::<f64>();
        let (lc, ld, lt) = (len(BoundaryKind::Convergent), len(BoundaryKind::Divergent), len(BoundaryKind::Transform));
        let crust_v: Vec<CrustType> = owner.iter().map(|o| tk[o]).collect();
        let (age, _) = crust::crust_age(&grid, owner, &crust_v, &segs, r);
        let mut oa = 0.0;
        let mut ow = 0.0;
        let mut max_o: f32 = 0.0;
        for c in grid.cells() {
            if crust_v[c as usize] == CrustType::Oceanic {
                let w = grid.area_sr(c);
                oa += f64::from(age[c as usize]) * w;
                ow += w;
                max_o = max_o.max(age[c as usize]);
            }
        }
        println!("seed {seed}: plates {} | largest {:.1}% 2nd {:.1}% 3rd {:.1}% smallest {:.2}%", sorted.len(), sorted[0], sorted[1], sorted[2], sorted.last().unwrap());
        println!("   continental area: v1 rule {:.1}%  target-fraction {:.1}%", v1_cont * 100.0, t_cont * 100.0);
        println!("   boundaries: {} segments, total {:.0} km (conv {:.0}, div {:.0}, transform {:.0}) [{:.0?}]", segs.len(), lc + ld + lt, lc, ld, lt, dt);
        let rates: Vec<f64> = segs.iter().filter(|s| s.kind != BoundaryKind::Transform).map(|s| s.rate_cm_yr).collect();
        println!(
            "   conv/div rates cm/yr: mean {:.1} max {:.1}; ocean crust mean age {:.0} Ma, max {:.0} Ma",
            rates.iter().sum::<f64>() / rates.len() as f64,
            rates.iter().cloned().fold(0.0, f64::max),
            oa / ow,
            max_o
        );
        if let Some(d) = &dump {
            let bf = boundaries::boundary_kind_field(grid.len(), &segs);
            let mut f = std::fs::File::create(format!("{d}/s{seed}.bin")).unwrap();
            for c in grid.cells() {
                let i = c as usize;
                f.write_all(&owner[i].to_le_bytes()).unwrap();
                f.write_all(&[bf[i], crust_v[i] as u8]).unwrap();
                f.write_all(&age[i].to_le_bytes()).unwrap();
            }
        }
    }
}
