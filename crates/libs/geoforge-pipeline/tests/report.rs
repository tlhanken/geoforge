//! The report page of a generated world: every image it shows is made, maps carry their
//! keys (the deposits' kinds, the chance ramps), places are picked by climate, and the sky
//! has its chart and key.

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
    let r = report::build(&w, &report::Options { width: 256, days: 3, ..report::Options::default() }).unwrap();

    // Every image the page shows is made, at the width asked for (keys aside).
    let made: Vec<&str> = r.images.iter().map(|(n, _)| n.as_str()).collect();
    for src in r.html.split("src=\"img/").skip(1) {
        let name = src.split(".jpg").next().unwrap();
        assert!(made.contains(&name), "{name} is shown but not made");
    }
    for (name, img) in &r.images {
        if !name.starts_with("sky_key") {
            assert_eq!(img.width(), 256, "{name}");
        }
    }
    for name in ["elevation", "provinces", "koppen", "deposits", "iron", "earthquakes", "weather_map", "sky", "sky_key_magnitude", "sky_key_colour"] {
        assert!(made.contains(&name), "no {name}");
    }

    // The deposits map names its kinds; chance maps carry their ramp.
    let fig = |name: &str| r.html.split(&format!("src=\"img/{name}.jpg\"")).nth(1).and_then(|s| s.split("</figure>").next()).unwrap_or_default().to_string();
    assert!(fig("deposits").contains("class=\"key\""), "{}", fig("deposits"));
    assert!(fig("iron").contains("class=\"ramp\"") && fig("iron").contains("100%"));
    assert!(fig("elevation").contains("class=\"ramp\""));
    assert!(fig("provinces").contains("Craton") || fig("provinces").contains("Platform"));

    // Places picked by climate, each with its weather.
    assert!(r.html.matches("class=\"card\"><h3>").count() >= 2, "too few places");
    assert!(r.html.contains("Low / high"));
    // The sky's key explains sizes, colours and rings.
    assert!(r.html.contains("visual magnitude") && r.html.contains("spectral class") && r.html.contains("green ring"));
    let _ = std::fs::remove_dir_all(&dir);
}
