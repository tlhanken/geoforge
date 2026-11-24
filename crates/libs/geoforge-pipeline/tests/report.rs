//! The report page of a generated world: it opens on the world whole in true colour, every
//! image it shows is made, and maps carry their keys (the deposits' kinds, the chance ramps).

use geoforge_core::grid::GridSpec;
use geoforge_core::Seed;
use geoforge_pipeline::{self as pipe, report, GenerateOptions, NewWorld};

#[test]
fn a_report_shows_every_map_with_its_key() {
    let dir = std::env::temp_dir().join(format!("geoforge-report-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    let new = NewWorld { recipe: geoforge_cosmos::Recipe::Habitable, seed: Seed(12345), grid: GridSpec::Equirect { width: 128, height: 64 }, body: None, system: None, focus: None };
    let mut w = pipe::create_world(&dir, &new).unwrap();
    pipe::generate(&mut w, &GenerateOptions::default(), |_| {}).unwrap();
    let r = report::build(&w, &report::Options { width: 256, full: true, ..report::Options::default() }).unwrap();

    // Every image the page shows is made, at the width asked for (keys aside).
    let made: Vec<&str> = r.images.iter().map(|(n, _)| n.as_str()).collect();
    for src in r.html.split("src=\"img/").skip(1) {
        let name = src.split(".jpg").next().unwrap();
        assert!(made.contains(&name), "{name} is shown but not made");
    }
    for (name, img) in &r.images {
        assert_eq!(img.width(), 256, "{name}");
    }
    for name in ["globe", "appearance", "elevation", "provinces", "koppen", "deposits", "iron", "earthquakes"] {
        assert!(made.contains(&name), "no {name}");
    }

    // The deposits map names its kinds; chance maps carry their ramp.
    let fig = |name: &str| r.html.split(&format!("src=\"img/{name}.jpg\"")).nth(1).and_then(|s| s.split("</figure>").next()).unwrap_or_default().to_string();
    assert!(fig("deposits").contains("class=\"key\""), "{}", fig("deposits"));
    assert!(fig("iron").contains("class=\"ramp\"") && fig("iron").contains("100%"));
    assert!(fig("elevation").contains("class=\"ramp\""));
    assert!(fig("provinces").contains("Craton") || fig("provinces").contains("Platform"));

    // The globe comes first, before any map; places, weather and sky are left to their commands.
    assert!(r.html.find("img/globe.jpg") < r.html.find("img/appearance.jpg"));
    assert!(!r.html.contains("id=\"places\"") && !r.html.contains("id=\"sky\""));

    // The lite report keeps only the main maps, and no system printout.
    let lite = report::build(&w, &report::Options { width: 256, ..report::Options::default() }).unwrap();
    let shown: Vec<&str> = lite.images.iter().map(|(n, _)| n.as_str()).filter(|n| *n != "globe").collect();
    assert!(shown.iter().all(|n| report::LITE.contains(n)) && shown.len() >= 10, "{shown:?}");
    assert!(!lite.html.contains("id=\"system\"") && r.html.contains("id=\"system\""));
    let _ = std::fs::remove_dir_all(&dir);
}
