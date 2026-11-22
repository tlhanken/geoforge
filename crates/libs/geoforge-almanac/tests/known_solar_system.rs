//! The real Solar System (`geoforge_cosmos::sol`) against the published almanac of 2026–27:
//! the eclipses (NASA, Espenak) and the planets' conjunctions and oppositions (NAOJ, Celestial
//! Phenomena in 2026, in Japan time, UTC + 9 h). Days count from the present, the 2026 March
//! equinox (14:46 UTC).

use geoforge_almanac::ephemeris::{eclipses, position_of, Coverage, EclipseKind, Equatorial};
use geoforge_core::body::SolarSystem;
use geoforge_core::{Archetype, Key, Seed};
use geoforge_cosmos::sol;

fn earth() -> (SolarSystem, Key) {
    let (s, body) = sol::archetype(Seed(1), Archetype::Earth);
    (s, body.key)
}

/// Days from the present to a Julian date (TT).
fn day(jd_tt: f64) -> f64 {
    jd_tt - sol::EQUINOX_2026_JD_TT
}

#[test]
fn the_present_is_the_march_equinox() {
    let (s, e) = earth();
    let sun = s.by_name("Sun").unwrap().key;
    let frame = Equatorial::of(&s, e).unwrap();
    let dec = |t: f64| frame.star_declination(&s, e, sun, t);
    assert!(dec(0.0).abs() < 0.01, "{}", dec(0.0));
    assert!(dec(1.0) > 0.3 && dec(-1.0) < -0.3);
    // The June solstice, 2026-06-21 08:24 UTC, is when the Sun stands highest.
    let june = day(2_461_212.850_8);
    assert!(dec(june) > 23.4 && dec(june) > dec(june - 2.0) && dec(june) > dec(june + 2.0), "{}", dec(june));
}

#[test]
fn eclipses_of_2026_and_2027() {
    let (s, e) = earth();
    let found = eclipses(&s, e, -40.0, 400.0);
    for (when, kind, coverage) in [
        (2_461_089.009_1, EclipseKind::Solar, None),                    // 2026-02-17 annular, grazing Antarctica
        (2_461_102.982_6, EclipseKind::Lunar, Some(Coverage::Total)),   // 2026-03-03
        (2_461_265.241_8, EclipseKind::Solar, Some(Coverage::Total)),   // 2026-08-12, Iceland and Spain
        (2_461_280.675_6, EclipseKind::Lunar, Some(Coverage::Partial)), // 2026-08-28
        (2_461_443.167_0, EclipseKind::Solar, Some(Coverage::Annular)), // 2027-02-06
    ] {
        // The Moon here keeps its mean orbit; its periodic terms (evection 1.3°, variation
        // 0.7°) move the real eclipses by up to a few hours.
        let hit = found.iter().find(|x| x.kind == kind && (x.day - day(when)).abs() * 24.0 < 4.0);
        let Some(hit) = hit else { panic!("no {kind:?} eclipse near day {:.2}: {found:#?}", day(when)) };
        if let Some(c) = coverage {
            assert_eq!(hit.coverage, c, "day {:.2}", hit.day);
        }
    }
    assert_eq!(found.len(), 5, "{found:#?}");
}

/// Geocentric ecliptic longitude of `name` minus the Sun's, degrees in (−180, 180].
fn from_sun(s: &SolarSystem, earth: Key, name: &str, t: f64) -> f64 {
    let at = |n: &str| position_of(s, s.by_name(n).unwrap().key, t);
    let (me, sun, p) = (position_of(s, earth, t), at("Sun"), at(name));
    let lon = |v: [f64; 3]| (v[1] - me[1]).atan2(v[0] - me[0]).to_degrees();
    (lon(p) - lon(sun) + 540.0).rem_euclid(360.0) - 180.0
}

#[test]
fn conjunctions_and_oppositions_of_2026() {
    let (s, e) = earth();
    // NAOJ's table is in Japan time.
    let jst = 9.0 / 24.0;
    for (name, when, target) in [
        ("Venus", 2_461_047.567_2, 0.0),     // superior conjunction 2026-01-07 01:36 JST
        ("Mars", 2_461_050.362_7, 0.0),      // 2026-01-09 20:41
        ("Jupiter", 2_461_051.238_3, 180.0), // opposition 2026-01-10 17:42
        ("Neptune", 2_461_122.347_0, 0.0),   // 2026-03-22 20:19
        ("Saturn", 2_461_125.247_2, 0.0),    // 2026-03-25 17:55
        ("Uranus", 2_461_183.477_4, 0.0),    // 2026-05-22 23:26
    ] {
        let off = |t: f64| (from_sun(&s, e, name, t) - target + 540.0).rem_euclid(360.0) - 180.0;
        let when = when - jst;
        let (mut lo, mut hi) = (day(when) - 10.0, day(when) + 10.0);
        assert!(off(lo) * off(hi) < 0.0, "{name}: no crossing near day {:.1}", day(when));
        for _ in 0..50 {
            let mid = 0.5 * (lo + hi);
            if off(lo) * off(mid) <= 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        let hours = (lo - day(when)) * 24.0;
        // Light time and the table's 15–600″ (Saturn's the worst) leave an hour or two.
        assert!(hours.abs() < 3.0, "{name}: {hours:.1} h off");
    }
}
