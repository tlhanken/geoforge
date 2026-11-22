//! End to end: create, generate, export, paint, import as a patch, regenerate downstream;
//! then author a layer outright. Small grid so it runs in seconds.

use geoforge_core::grid::GridSpec;
use geoforge_core::product::{LayerId, ProductData, ProductId};
use geoforge_core::Seed;
use geoforge_pipeline::{self as pipe, GenerateOptions, LayerAction, NewWorld, State};
use geoforge_store::World;

fn tmp(name: &str) -> std::path::PathBuf {
    let d = std::env::temp_dir().join(format!("geoforge-pipe-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    d
}

fn ran(reports: &[pipe::LayerReport]) -> Vec<LayerId> {
    reports.iter().filter(|r| matches!(r.action, LayerAction::Ran { .. })).map(|r| r.layer).collect()
}

#[test]
fn generate_edit_regenerate() {
    let dir = tmp("cycle");
    let mut w =
        pipe::create_world(&dir, &NewWorld { recipe: geoforge_cosmos::Recipe::Habitable, seed: Seed(7), grid: GridSpec::Equirect { width: 256, height: 128 }, body: None, focus: None }).unwrap();

    // First run builds everything; a second run reuses everything.
    let (p1, r1) = pipe::generate(&mut w, &GenerateOptions::default(), |_| {}).unwrap();
    assert_eq!(ran(&r1), LayerId::ALL.into_iter().filter(|l| *l != LayerId::Bands).collect::<Vec<_>>());
    let (_, r2) = pipe::generate(&mut w, &GenerateOptions::default(), |_| {}).unwrap();
    assert!(ran(&r2).is_empty(), "nothing should rerun: {:?}", ran(&r2));
    // Reopening from disk gives the same products.
    let w2 = World::open(&dir).unwrap();
    let loaded = w2.load_products().unwrap();
    let differ: Vec<ProductId> = ProductId::ALL.into_iter().filter(|&p| loaded.get(p) != p1.get(p)).collect();
    assert!(differ.is_empty(), "differ after reload: {differ:?}");

    // Export the rough heightmap, raise a block by ~2 km, import as a patch.
    let edit = pipe::edit_folder(&w);
    pipe::export(&w, &[ProductId::RoughElevation], &edit, &Default::default()).unwrap();
    let png = edit.join("rough_elevation.png");
    let mut img = image::open(&png).unwrap().into_luma16();
    let step = 24000.0 / 65535.0;
    for y in 60..68 {
        for x in 100..108 {
            let px = img.get_pixel_mut(x, y);
            px.0[0] = (f64::from(px.0[0]) + 2000.0 / step).round().min(65535.0) as u16;
        }
    }
    img.save(&png).unwrap();
    let rep = pipe::import(&mut w, ProductId::RoughElevation, &edit, false).unwrap();
    assert_eq!(rep.mode, "patch");
    assert_eq!(rep.changed_cells, 64);

    // Only elevation and what depends on it rerun; the edit is in the result.
    let plan = pipe::plan(&w).unwrap();
    let stale: Vec<LayerId> = plan.layers.iter().filter(|l| l.needs_run()).map(|l| l.layer).collect();
    assert_eq!(stale, vec![LayerId::Elevation, LayerId::Circulation, LayerId::Climate, LayerId::Erosion, LayerId::Hydrology, LayerId::Biomes, LayerId::Resources]);
    let (p3, r3) = pipe::generate(&mut w, &GenerateOptions::default(), |_| {}).unwrap();
    assert_eq!(ran(&r3), stale);
    let (Some(e1), Some(e3)) = (p1.rough_elevation(), p3.rough_elevation()) else { panic!() };
    let c = 64 * 256 + 104;
    assert!((e3.values[c] - e1.values[c] - 2000.0).abs() < 1.0, "{} {}", e1.values[c], e3.values[c]);
    assert_eq!(e3.values[0], e1.values[0]);
    // Erosion carries the raised block into the final heightmap.
    assert!(p3.elevation().unwrap().values[c] > p1.elevation().unwrap().values[c] + 1000.0);
    // Higher ground is colder.
    assert!(p3.temperature().unwrap().values[c] < p1.temperature().unwrap().values[c] - 8.0);

    // Changing an upstream parameter regenerates under the edit: the patch rides on top.
    w.manifest.params.get_mut("elevation").unwrap()["continental_base_m"] = serde_json::json!(300.0);
    w.save_manifest().unwrap();
    let (p4, _) = pipe::generate(&mut w, &GenerateOptions::default(), |_| {}).unwrap();
    assert!(pipe::edit_status(&w).get("rough_elevation").is_some_and(|s| s.rebased));
    let e4 = p4.rough_elevation().unwrap();
    assert!(e4.values[c] > e1.values[c] + 1500.0);

    // Author the crust type outright (all continental): tectonics keeps it, the rest follows.
    pipe::export(&w, &[ProductId::CrustType], &edit, &Default::default()).unwrap();
    let cpng = edit.join("crust_type.png");
    let mut ci = image::open(&cpng).unwrap().into_rgb8();
    for px in ci.pixels_mut() {
        px.0 = [196, 160, 0];
    }
    ci.save(&cpng).unwrap();
    let rep = pipe::import(&mut w, ProductId::CrustType, &edit, true).unwrap();
    assert_eq!(rep.mode, "authored");
    let (p5, _) = pipe::generate(&mut w, &GenerateOptions::default(), |_| {}).unwrap();
    let ProductData::CrustType(ct) = p5.get(ProductId::CrustType).unwrap() else { panic!() };
    assert!(ct.values.iter().all(|c| *c == geoforge_core::surface::CrustType::Continental));
    let plan = pipe::plan(&w).unwrap();
    assert!(plan.layers.iter().filter(|l| l.not_applicable.is_none()).flat_map(|l| &l.products).all(|p| p.state == State::Fresh));

    // Reset drops the edit; it is generated again.
    pipe::reset(&mut w, ProductId::CrustType).unwrap();
    assert!(pipe::plan(&w).unwrap().layers[0].needs_run());
}
