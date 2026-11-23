//! A place's almanac on a generated world: a whole year of days with sun times and moons, its
//! months read against the drought record, the same almanac every time, and the same days
//! in every calendar.

use geoforge_core::grid::GridSpec;
use geoforge_core::Seed;
use geoforge_pipeline::{self as pipe, GenerateOptions, NewWorld};

#[test]
fn a_year_at_a_place_day_by_day() {
    let dir = std::env::temp_dir().join(format!("geoforge-almanac-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    let mut w =
        pipe::create_world(&dir, &NewWorld { recipe: geoforge_cosmos::Recipe::Habitable, seed: Seed(3), grid: GridSpec::Equirect { width: 128, height: 64 }, body: None, system: None, focus: None })
            .unwrap();
    pipe::generate(&mut w, &GenerateOptions::default(), |_| {}).unwrap();
    // A land cell at mid-latitude.
    let grid = w.manifest.grid.build().unwrap();
    let p = w.load_products().unwrap();
    let e = &p.elevation().unwrap().values;
    let c = (0..e.len()).find(|&c| e[c] > 100.0 && (30.0..55.0).contains(&geoforge_core::grid::Grid::center(&grid, c as u32).to_lat_lon().lat_deg.abs())).expect("mid-latitude land") as u32;
    let ll = geoforge_core::grid::Grid::center(&grid, c).to_lat_lon();

    use geoforge_almanac::calendar::{self, Choice};
    let twelfths = calendar::of_world(&w, &Choice::World).unwrap();
    assert_eq!(twelfths.name(), "twelfths", "a world without its own calendar");
    let a = geoforge_almanac::year::year(&w, ll.lat_deg, ll.lon_deg, 1, twelfths.as_ref()).unwrap();
    let body = w.load_body().unwrap();
    let cal = geoforge_almanac::ephemeris::Calendar::new(&body, &[]);
    let per_year = cal.days_per_year().unwrap();
    assert!((a.days.len() as f64 - per_year).abs() <= 1.0, "{} days, {per_year} a year", a.days.len());
    assert!(a.days.iter().all(|d| d.weather.year == 1 && d.date.year == 1));
    // Dates run day by day through the months.
    for p in a.days.windows(2) {
        let (x, y) = (p[0].date, p[1].date);
        assert!((y.month == x.month && y.day == x.day + 1) || (y.month == x.month + 1 && y.day == 0), "{x:?} then {y:?}");
    }
    // Days lengthen and shorten through the year.
    let lengths: Vec<f64> = a.days.iter().filter_map(|d| d.sun.map(|(r, s)| s - r)).collect();
    let (short, long) = lengths.iter().fold((f64::MAX, f64::MIN), |(lo, hi), &l| (lo.min(l), hi.max(l)));
    assert!(long - short > 0.1 * a.day_h.unwrap(), "{short} to {long} h");
    assert_eq!(a.months.len(), 12);
    // The months' rain is the days' rain.
    let days_rain: f64 = a.days.iter().map(|d| d.weather.precipitation_mm).sum();
    let months_rain: f64 = a.months.iter().map(|m| m.rain_mm).sum();
    assert!((days_rain - months_rain).abs() < 0.15 * days_rain.max(1.0), "{days_rain} vs {months_rain}");
    assert_eq!(a.droughts.years, geoforge_almanac::year::DROUGHT_YEARS);
    assert_eq!(a, geoforge_almanac::year::year(&w, ll.lat_deg, ll.lon_deg, 1, twelfths.as_ref()).unwrap());

    // The world's own calendar names the same days.
    let json = format!(
        r#"{{"name": "Reckoning", "months": [{{"name": "Thaw", "days": {}}}, {{"name": "Frost"}}], "week": ["One", "Two", "Three"], "present_year": 400, "era": "AR"}}"#,
        (per_year / 2.0).floor()
    );
    std::fs::write(w.body_dir().join("calendar.json"), json).unwrap();
    let own = calendar::of_world(&w, &Choice::World).unwrap();
    let b = geoforge_almanac::year::year(&w, ll.lat_deg, ll.lon_deg, 1, own.as_ref()).unwrap();
    assert_eq!((b.calendar.as_str(), b.year_name.as_str(), b.calendar_months.len()), ("Reckoning", "401 AR", 2));
    assert!(a.days.iter().zip(&b.days).all(|(x, y)| x.weather == y.weather), "a calendar changes only the names");
    assert_eq!(b.days.iter().filter(|d| d.date.month == 0).count() as f64, (per_year / 2.0).floor());
    let _ = std::fs::remove_dir_all(&dir);
}
