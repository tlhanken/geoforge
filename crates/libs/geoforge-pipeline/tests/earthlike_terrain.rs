//! An Earth-like world's terrain and drainage, generated end to end, against Earth's: how much
//! land stands high, how tall its peaks are, how its land parts into drainage basins, and how
//! much of it drains to no sea.

use geoforge_core::grid::GridSpec;
use geoforge_core::{Archetype, Key, Seed};
use geoforge_pipeline::{self as pipe, GenerateOptions, NewWorld};

#[test]
fn an_earthlike_world_stands_and_drains_like_earth() {
    let dir = std::env::temp_dir().join(format!("geoforge-earthlike-terrain-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    let body = Archetype::Earth.context(Key(1), "Earth");
    let mut w = pipe::create_world(
        &dir,
        &NewWorld { recipe: geoforge_cosmos::Recipe::Habitable, seed: Seed(42), grid: GridSpec::Equirect { width: 1024, height: 512 }, body: Some(body.clone()), system: None, focus: None },
    )
    .unwrap();
    let (products, _) = pipe::generate(&mut w, &GenerateOptions::default(), |_| {}).unwrap();
    let grid = w.manifest.grid.build().unwrap();
    let s = pipe::stats::summarize(&products, &grid, body.exterior.radius_km, true);
    let v = |sec: &str, k: &str| s[sec][k].value;
    // Earth: 27% of land above 1000 m (its ice sheets included), its highest 8.85 km.
    assert!((0.15..0.40).contains(&v("elevation", "land_above_1000m_fraction")), "{}", v("elevation", "land_above_1000m_fraction"));
    assert!((7000.0..9000.0).contains(&v("elevation", "max_elevation_m")), "{}", v("elevation", "max_elevation_m"));
    // The largest basin no more than twice the Amazon's share of the land (7 Mkm² of 149).
    let land_km2 = v("elevation", "land_fraction") * 4.0 * std::f64::consts::PI * body.exterior.radius_km.powi(2);
    let largest = v("hydrology", "largest_basin_km2") / land_km2;
    assert!(largest < 2.0 * 7.0 / 149.0, "the largest basin drains {:.1}% of the land", 100.0 * largest);
    // Earth: 18% of land drains to no sea.
    assert!((0.08..0.3).contains(&v("hydrology", "endorheic_fraction_of_land")), "{}", v("hydrology", "endorheic_fraction_of_land"));
    let _ = std::fs::remove_dir_all(&dir);
}
