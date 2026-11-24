//! How a planet turns once its system has settled: whether its tilt holds or wanders, and
//! whether its air holds it off the one face to its star that tides would give it.
//!
//! **Tilt.** A planet's tilt is set at birth by its last giant impacts; what happens after
//! depends on how fast its axis precesses. A big moon speeds the precession past the
//! frequencies at which the other planets tug on the orbit, and the tilt stays: Earth's varies
//! by ±1.3° (Laskar, Joutel & Robutel 1993). Without one it can meet those frequencies and
//! wander: Mars between 0 and 60° (Laskar et al. 2004), a moonless Earth over 0–85° by Laskar &
//! Robutel (1993), but only 10–20° either way over billions of years by Lissauer, Barnes &
//! Chambers (2012). Between those views, a moonless world's tilt today lies anywhere within
//! 10–45° of its birth tilt. A moon counts as big when its pull on the planet's bulge, G·m/a³,
//! is at least the star's, G·M/a_p³ (the Moon's is 2.2 times the Sun's on Earth).
//!
//! **Thermal tides.** The star's heat raises a bulge in the air that lags noon, and the star's
//! pull on it spins the planet away from one face (Venus's retrograde day). It wins against
//! the solid body's tide beyond a_c ≈ (0.4 − 0.12·M/M☉)·(1 + 0.8e²) AU for a 1-bar air
//! (Leconte et al. 2015; Valente et al. 2024, their fit and their eccentric runs), where the
//! solar day comes in at ~30 days (the thermal tide's peak, Leconte's σ₀ = 145 per year) and
//! shortens to one day by 0.62 AU around a 0.82 M☉ star: 30 d·(a/a_c)^−4.7. That is worked
//! out for air near a bar, and is used from 0.5 bar up; with less the planet locks. Thicker
//! air holds its heat longer (its radiative time grows with its mass), so the bulge lags less
//! and the day it settles to is longer: past 5 bar the day grows as (p/5 bar)², which brings
//! Venus's 92 bar at 0.72 AU to its 117-day solar day. Up to 5 bar the planet keeps the
//! prograde of the two stable states, the one it falls into while being slowed; thick air
//! drives it toward either (Correia & Laskar 2001, 2003: which one turns on the planet's
//! early spin, tilt and the friction between its core and mantle), so a thick-aired world
//! turns backward half the time, as Venus does (tilt 177°).

use crate::context::body_context;
use crate::draw::Draw;
use geoforge_core::body::{BodyExterior, NodeKind, OrbitalNode, SolarSystem};
use geoforge_core::key::Key;

/// Birth-to-now reach of a moonless world's tilt, degrees: from Lissauer et al. to Laskar.
const WANDER_DEG: (f64, f64) = (10.0, 45.0);
/// Air (bar) for which the thermal tide is worked out, and past which it is thick.
const THERMAL_TIDE_BAR: (f64, f64) = (0.5, 5.0);
/// How fast the day lengthens with thick air's pressure, (p/5 bar)^this.
const THICK_AIR_SLOPE: f64 = 2.0;
/// Solar day where the thermal tide first wins, days.
const THERMAL_TIDE_DAY_D: f64 = 30.0;
/// How fast the solar day shortens outward of a_c.
const THERMAL_TIDE_SLOPE: f64 = 4.7;

/// Whether a moon of `moon_kg` at `moon_a_km` steadies the tilt of a planet at `planet_a_km`
/// from a star of `star_kg`.
#[must_use]
pub fn moon_steadies(moon_kg: f64, moon_a_km: f64, star_kg: f64, planet_a_km: f64) -> bool {
    moon_kg / moon_a_km.powi(3) >= star_kg / planet_a_km.powi(3)
}

/// Today's tilt (degrees) of a world born with `birth_deg` and no big moon.
pub(crate) fn wandered_tilt(d: &mut Draw, birth_deg: f64) -> f64 {
    let reach = d.range(WANDER_DEG.0, WANDER_DEG.1);
    d.range((birth_deg - reach).max(0.0), (birth_deg + reach).min(90.0))
}

/// Inner edge (AU) of the thermal-tide zone around a star of `star_mass` M☉, for an orbit of
/// eccentricity `e`.
#[must_use]
pub fn thermal_tide_edge_au(star_mass: f64, e: f64) -> f64 {
    (0.4 - 0.12 * star_mass).max(0.05) * (1.0 + 0.8 * e * e)
}

/// Solar day (hours) that the thermal tide holds a planet to, at `a_au` with eccentricity `e`
/// from a star of `star_mass` M☉, under `pressure_bar` of air: `None` if it locks after all.
#[must_use]
pub fn thermal_tide_solar_day_h(a_au: f64, e: f64, star_mass: f64, pressure_bar: f64) -> Option<f64> {
    let edge = thermal_tide_edge_au(star_mass, e);
    (pressure_bar >= THERMAL_TIDE_BAR.0 && a_au > edge)
        .then(|| 24.0 * THERMAL_TIDE_DAY_D * (a_au / edge).powf(-THERMAL_TIDE_SLOPE) * (pressure_bar / THERMAL_TIDE_BAR.1).max(1.0).powf(THICK_AIR_SLOPE))
}

/// Sidereal rotation (hours) of such a planet turning forward, with a year of `year_h` hours.
#[must_use]
pub fn thermal_tide_rotation_h(a_au: f64, e: f64, star_mass: f64, pressure_bar: f64, year_h: f64) -> Option<f64> {
    thermal_tide_solar_day_h(a_au, e, star_mass, pressure_bar).map(|solar_h| 1.0 / (1.0 / solar_h + 1.0 / year_h))
}

/// Spin of such a planet turning backward: its sidereal rotation (hours) and whether that
/// spin, measured along the orbit's pole, is backward (a tilt past 90°). A solar day longer
/// than the year leaves the planet turning slowly forward.
#[must_use]
pub fn backward_rotation_h(solar_h: f64, year_h: f64) -> (f64, bool) {
    let spin = 1.0 / year_h - 1.0 / solar_h;
    (1.0 / spin.abs().max(1e-12), spin < 0.0)
}

/// Rotation (hours) of a world tidally locked on an orbit of `year_h` hours and eccentricity
/// `e`: once per orbit, or Mercury's 3:2 when the orbit is eccentric.
#[must_use]
pub fn locked_rotation_h(year_h: f64, e: f64) -> f64 {
    if e > 0.1 {
        year_h * 2.0 / 3.0
    } else {
        year_h
    }
}

/// Change the body at `key`.
pub(crate) fn set_body(nodes: &mut [OrbitalNode], key: Key, f: impl FnOnce(&mut BodyExterior)) {
    if let Some(NodeKind::Body(b)) = nodes.iter_mut().find(|n| n.key == key).map(|n| &mut n.node) {
        f(b);
    }
}

/// Free the locked rocky planets of `system` (keys, with their hosts' masses in M☉) whose air's
/// thermal tide holds them off; a thick-aired one turns backward or forward by a draw of its own.
pub(crate) fn thermal_tides(system: &mut SolarSystem, locked: &[(Key, f64)]) {
    for &(key, host_mass) in locked {
        let Some(ctx) = body_context(system, key) else { continue };
        let year_h = ctx.irradiance.year_days * 24.0;
        let pressure = ctx.interior.surface_pressure_bar();
        let Some(solar_h) = thermal_tide_solar_day_h(ctx.orbit.a_au, ctx.orbit.e, host_mass, pressure) else { continue };
        let backward = pressure > THERMAL_TIDE_BAR.1 && Draw::new(system.root.seed.child_u64(b"thermal-tide", key.0)).u() < 0.5;
        let (rotation_h, flip) = if backward { backward_rotation_h(solar_h, year_h) } else { (1.0 / (1.0 / solar_h + 1.0 / year_h), false) };
        set_body(&mut system.nodes, key, |b| {
            b.rotation_h = rotation_h;
            if flip {
                b.obliquity_deg = 180.0 - b.obliquity_deg;
            }
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use geoforge_core::body::{AU_KM, EARTH_MASS_KG, SOLAR_MASS_KG};
    use geoforge_core::seed::Seed;

    #[test]
    fn the_moon_steadies_earth() {
        assert!(moon_steadies(0.0123 * EARTH_MASS_KG, 384_400.0, SOLAR_MASS_KG, AU_KM));
        // A Moon a third as heavy and half again as far would not.
        assert!(!moon_steadies(0.004 * EARTH_MASS_KG, 576_600.0, SOLAR_MASS_KG, AU_KM));
    }

    #[test]
    fn moonless_tilts_wander_but_stay_prograde() {
        let mut d = Draw::new(Seed(4));
        let tilts: Vec<f64> = (0..2000).map(|_| wandered_tilt(&mut d, 23.4)).collect();
        assert!(tilts.iter().all(|t| (0.0..=68.4).contains(t)));
        let far = tilts.iter().filter(|t| (*t - 23.4).abs() > 10.0).count() as f64 / 2000.0;
        assert!((0.3..0.8).contains(&far), "{far} moved more than 10°");
    }

    #[test]
    fn thermal_tides_free_k_dwarf_worlds() {
        // An Earth twin around a 0.82 M☉ star at 0.62 AU turns once a day (Valente et al. 2024).
        let year_h = 0.62f64.powf(1.5) / 0.82f64.sqrt() * 8766.0;
        let rot = thermal_tide_rotation_h(0.62, 0.0, 0.82, 1.0, year_h).unwrap();
        let solar = 1.0 / (1.0 / rot - 1.0 / year_h);
        assert!((18.0..32.0).contains(&solar), "solar day {solar} h");
        // Just past the edge the day is about a month; inside it, or with little air, locked.
        let edge = thermal_tide_edge_au(0.7, 0.0);
        assert!((0.3..0.34).contains(&edge), "{edge}");
        let near = thermal_tide_rotation_h(1.02 * edge, 0.0, 0.7, 1.0, 2000.0).unwrap();
        assert!((20.0..30.0).contains(&(1.0 / (1.0 / near - 1.0 / 2000.0) / 24.0)), "{near} h");
        assert!(thermal_tide_rotation_h(0.25, 0.0, 0.7, 1.0, 2000.0).is_none());
        assert!(thermal_tide_rotation_h(0.5, 0.0, 0.7, 0.1, 2000.0).is_none());
        // Eccentric orbits need to be farther out (Valente et al.: 0.35 AU needs e under 0.33).
        assert!(thermal_tide_edge_au(0.7, 0.33) > 1.08 * edge);
    }

    /// Venus: 92 bar at 0.72 AU gives its 117-day solar day, and turning backward, its
    /// 243-day sidereal day; air past 5 bar lengthens the day.
    #[test]
    fn thick_air_gives_venus_its_day() {
        let solar = thermal_tide_solar_day_h(0.7233, 0.0068, 1.0, 92.0).unwrap() / 24.0;
        assert!((100.0..135.0).contains(&solar), "solar day {solar} d");
        let (sidereal, backward) = backward_rotation_h(solar * 24.0, 224.7 * 24.0);
        assert!(backward && (215.0..275.0).contains(&(sidereal / 24.0)), "{} d", sidereal / 24.0);
        let at = |p: f64| thermal_tide_solar_day_h(0.5, 0.0, 0.7, p).unwrap();
        assert!((at(5.0) - at(1.0)).abs() < 1e-9 && at(20.0) > 10.0 * at(5.0));
        // A day longer than the year leaves the planet turning slowly forward.
        assert!(!backward_rotation_h(500.0 * 24.0, 224.7 * 24.0).1);
    }
}
