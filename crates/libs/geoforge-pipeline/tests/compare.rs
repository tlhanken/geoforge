//! The A/B page end to end: an "observed" world whose temperature is authored, set against
//! the model world it was copied from, and a second model that only copied the input.

use geoforge_core::grid::GridSpec;
use geoforge_core::product::ProductId;
use geoforge_core::Seed;
use geoforge_pipeline::{self as pipe, compare, GenerateOptions, NewWorld};
use geoforge_store::World;
use std::path::{Path, PathBuf};

fn tmp(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("geoforge-compare-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    d
}

fn copy(from: &Path, to: &Path) {
    std::fs::create_dir_all(to).unwrap();
    for e in std::fs::read_dir(from).unwrap() {
        let e = e.unwrap();
        let dst = to.join(e.file_name());
        if e.file_type().unwrap().is_dir() {
            copy(&e.path(), &dst);
        } else {
            std::fs::copy(e.path(), dst).unwrap();
        }
    }
}

#[test]
fn observed_products_are_compared_and_copied_inputs_are_not() {
    let model_dir = tmp("model");
    let mut model = pipe::create_world(
        &model_dir,
        &NewWorld { recipe: geoforge_cosmos::Recipe::Habitable, seed: Seed(3), grid: GridSpec::Equirect { width: 128, height: 64 }, body: None, system: None, focus: None },
    )
    .unwrap();
    pipe::generate(&mut model, &GenerateOptions::default(), |_| {}).unwrap();

    // The "observed" world: the same temperature map, warmer, authored.
    let real_dir = tmp("real");
    copy(&model_dir, &real_dir);
    let mut real = World::open(&real_dir).unwrap();
    let edit = pipe::edit_folder(&real);
    pipe::export(&real, &[ProductId::Temperature], &edit, &Default::default()).unwrap();
    let png = edit.join("temperature.png");
    let mut img = image::open(&png).unwrap().into_luma16();
    for px in img.pixels_mut() {
        px.0[0] = px.0[0].saturating_add(3000);
    }
    img.save(&png).unwrap();
    assert_eq!(pipe::import(&mut real, ProductId::Temperature, &edit, true).unwrap().mode, "authored");
    pipe::generate(&mut real, &GenerateOptions::default(), |_| {}).unwrap();

    // A model that took the observed temperature as its input.
    let copied_dir = tmp("copied");
    copy(&real_dir, &copied_dir);
    let copied = World::open(&copied_dir).unwrap();

    let model = World::open(&model_dir).unwrap();
    let models = [compare::Model { label: "model".into(), world: &model }, compare::Model { label: "copied".into(), world: &copied }];
    // A place on land (temperature is compared over land).
    let grid = model.manifest.grid.build().unwrap();
    let elev = model.load_products().unwrap();
    let Some(geoforge_core::product::ProductData::Elevation(e)) = elev.get(ProductId::Elevation) else { panic!("no elevation") };
    let land = (0..e.values.len()).find(|&c| e.values[c] > 200.0).expect("some land") as u32;
    let ll = geoforge_core::grid::Grid::center(&grid, land).to_lat_lon();
    let probes = vec![compare::Probe { name: "inland".into(), lat: (ll.lat_deg, ll.lat_deg), lon: (ll.lon_deg, ll.lon_deg) }];
    let c = compare::build(&real, &models, &compare::Options { width: 256, probes, sources: vec![("temperature".into(), "a warmer copy".into())] }).unwrap();

    assert_eq!(c.sections.iter().map(|s| s.name.as_str()).collect::<Vec<_>>(), ["temperature"]);
    let s = &c.sections[0];
    assert_eq!(s.models.len(), 1, "the copied input is not compared");
    let compare::Metrics::Scalar { bias, corr, over, .. } = &s.models[0].1 else { panic!("{:?}", s.models) };
    assert_eq!(*over, "land");
    assert!(*bias < -1.0 && *corr > 0.95, "bias {bias}, corr {corr}");
    // The probe: observed, then the model, colder.
    let (_, v) = &s.probes[0];
    assert!(v[1].unwrap() < v[0].unwrap());
    // Maps: observed and the model, with the observed key; the numbers in the page.
    assert_eq!(c.images.iter().map(|(n, _)| n.as_str()).collect::<Vec<_>>(), ["real_temperature", "model_temperature"]);
    assert!(c.html.contains(r#"id="metrics""#) && c.html.contains("Annual mean temperature") && c.html.contains("a warmer copy"));
}
