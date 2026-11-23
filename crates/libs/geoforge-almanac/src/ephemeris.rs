//! Where everything in a system is at a moment, and the calendar a world keeps.
//!
//! Each node moves on a fixed Kepler orbit around its parent from its elements, which hold
//! at the system's present (t = 0); time is in Earth days. Around a barycentre a star's
//! elements describe its own orbit about it, so its period is the pair's relative orbit's.
//! Planets do not perturb each other here.
//!
//! Eclipses are found by stepping the geometry as seen from the world's centre: a moon (or,
//! for a moon, its planet) crossing a star's disc is a solar eclipse, seen from part of the
//! world when the separation is within the two discs plus the parallax of the world's own
//! radius; a moon passing through the world's umbra is a lunar eclipse. Orbits are fixed
//! unless a node says how it turns (the real Moon's node and perigee do, so its eclipse
//! seasons drift over the 18.6-year nodal cycle; drawn moons' orbits stay fixed).
//!
//! A world's solar day follows from its spin and its year, 1/day = 1/spin − 1/year (plus for
//! a retrograde spinner such as Venus, whose tilt is over 90°); a world locked to its star has
//! none. A month is a moon's synodic period, new moon to new moon: 1/month = 1/orbit − 1/year.

use geoforge_core::body::{NodeKind, SolarSystem};
use geoforge_core::derived::orbital_period_days;
use geoforge_core::key::Key;
use geoforge_core::time::{Clock, Span};
use geoforge_core::{BodyContext, Elements};
use std::collections::BTreeMap;

/// A node's place at a moment.
#[derive(Debug, Clone, PartialEq)]
pub struct Placed {
    pub key: Key,
    pub name: String,
    pub parent: Option<Key>,
    /// Position relative to the system's root (its star or barycentre), AU, in the plane the
    /// orbits are measured from.
    pub position_au: [f64; 3],
    /// Distance from the parent, AU.
    pub distance_au: f64,
    /// Orbital period around the parent, Earth days.
    pub period_days: f64,
}

/// Orbital period of `key` around its parent, Earth days; `None` for a node with no orbit.
#[must_use]
pub fn period_days(system: &SolarSystem, key: Key) -> Option<f64> {
    let node = system.node(key)?;
    let (e, parent) = (node.elements?, node.parent?);
    let m = system.mass_of(key);
    let is_star = |n: &geoforge_core::body::OrbitalNode| matches!(n.node, NodeKind::Star(_));
    let barycentre = matches!(system.node(parent).map(|p| &p.node), Some(NodeKind::Barycenter));
    // What is orbited: a barycentre's stars, or the parent.
    let central = if barycentre { system.children(parent).into_iter().filter(|c| is_star(c)).map(|c| system.mass_of(c.key)).sum() } else { system.mass_of(parent) };
    // A member star's elements are about the barycentre: the pair's relative orbit is wider
    // by M/(M − m) and goes round the pair's whole mass. Otherwise the two-body period.
    let (a_rel, total) = if barycentre && is_star(node) { (e.a_au * central / (central - m).max(central * 1e-9), central) } else { (e.a_au, central + m) };
    Some(orbital_period_days(a_rel, total))
}

/// Position relative to the parent `t` days after the present, AU.
#[must_use]
pub fn relative_position(e: &Elements, period_days: f64, t_days: f64) -> [f64; 3] {
    e.position_at(e.m0_deg + 360.0 * t_days / period_days)
}

/// Every star and body `t` days after the present (belts and rings have no single place).
#[must_use]
pub fn positions(system: &SolarSystem, t_days: f64) -> Vec<Placed> {
    let mut at: BTreeMap<Key, [f64; 3]> = BTreeMap::new();
    let mut out = Vec::new();
    for n in &system.nodes {
        if matches!(n.node, NodeKind::Belt(_) | NodeKind::Ring(_)) {
            continue;
        }
        let pos = place(system, n.key, t_days, &mut at);
        let (rel, period) = match (orbit_at(n, t_days), period_days(system, n.key)) {
            (Some(e), Some(p)) => (relative_position(&e, p, t_days), p),
            _ => ([0.0; 3], 0.0),
        };
        out.push(Placed { key: n.key, name: n.name.clone(), parent: n.parent, position_au: pos, distance_au: norm(rel), period_days: period });
    }
    out
}

fn place(system: &SolarSystem, key: Key, t: f64, at: &mut BTreeMap<Key, [f64; 3]>) -> [f64; 3] {
    if let Some(&p) = at.get(&key) {
        return p;
    }
    let node = system.node(key).expect("a node of the system");
    let pos = match (node.parent, orbit_at(node, t), period_days(system, key)) {
        (Some(parent), Some(e), Some(period)) => {
            let base = place(system, parent, t, at);
            let rel = relative_position(&e, period, t);
            [base[0] + rel[0], base[1] + rel[1], base[2] + rel[2]]
        }
        (Some(parent), _, _) => place(system, parent, t, at),
        _ => [0.0; 3],
    };
    at.insert(key, pos);
    pos
}

/// A node's elements `t` days after the present, turned by its precession.
fn orbit_at(node: &geoforge_core::body::OrbitalNode, t: f64) -> Option<Elements> {
    node.elements.map(|e| node.precession.map_or(e, |p| e.precessed(p, t)))
}

fn norm(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// A body's centre `t` days after the present, AU from the system's root.
#[must_use]
pub fn position_of(system: &SolarSystem, key: Key, t_days: f64) -> [f64; 3] {
    place(system, key, t_days, &mut BTreeMap::new())
}

/// TT − UTC since 2017 (32.184 s plus 37 leap seconds), days.
const TT_MINUS_UTC_DAYS: f64 = 69.184 / 86_400.0;

/// A Julian date (TT) as a UTC calendar date and time, "2026-03-20 14:46 UTC" (Meeus,
/// *Astronomical Algorithms*, ch. 7).
#[must_use]
pub fn utc_date(jd_tt: f64) -> String {
    let jd = jd_tt - TT_MINUS_UTC_DAYS + 0.5;
    let (z, f) = (jd.floor(), jd.fract());
    let alpha = ((z - 1_867_216.25) / 36_524.25).floor();
    let a = if z < 2_299_161.0 { z } else { z + 1.0 + alpha - (alpha / 4.0).floor() };
    let b = a + 1524.0;
    let c = ((b - 122.1) / 365.25).floor();
    let d = (365.25 * c).floor();
    let e = ((b - d) / 30.6001).floor();
    let day = b - d - (30.6001 * e).floor();
    let month = if e < 14.0 { e - 1.0 } else { e - 13.0 };
    let year = if month > 2.0 { c - 4716.0 } else { c - 4715.0 };
    let minutes = (f * 1440.0).round() as i64;
    // Rounding to the minute can carry into the next day; days are whole here, so carry it.
    let (day, minutes) = if minutes == 1440 { (day + 1.0, 0) } else { (day, minutes) };
    format!("{year:.0}-{month:02.0}-{day:02.0} {:02}:{:02} UTC", minutes / 60, minutes % 60)
}

/// The real calendar date `t_days` after the present, for a system that keeps one (the Solar
/// System).
#[must_use]
pub fn real_date(system: &SolarSystem, t_days: f64) -> Option<String> {
    system.root.present_jd_tt.map(|jd| utc_date(jd + t_days))
}

/// What a world's seasons follow: the nearest star or barycentre up the tree.
#[must_use]
pub fn light_source(system: &SolarSystem, world: Key) -> Option<Key> {
    let mut up = system.node(world)?.parent;
    while let Some(k) = up {
        let n = system.node(k)?;
        if matches!(n.node, NodeKind::Star(_) | NodeKind::Barycenter) {
            return Some(k);
        }
        up = n.parent;
    }
    None
}

/// A world's equatorial frame in the system frame: `z` its north pole, `x` toward its
/// northern spring equinox, `y = z × x`. The present is that equinox: the star (or the pair
/// its orbit circles) lies along `x`, crossing the equator northward.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Equatorial {
    pub x: [f64; 3],
    pub y: [f64; 3],
    pub z: [f64; 3],
}

impl Equatorial {
    /// The frame of `world`, if it orbits something. Its pole leans by its obliquity from the
    /// normal of its own orbit, away from its direction of motion at the equinox.
    #[must_use]
    pub fn of(system: &SolarSystem, world: Key) -> Option<Self> {
        let node = system.node(world)?;
        let (NodeKind::Body(b), Some(e)) = (&node.node, node.elements) else {
            return None;
        };
        let (lan, inc) = (e.lan_deg.to_radians(), e.i_deg.to_radians());
        let normal = [lan.sin() * inc.sin(), -lan.cos() * inc.sin(), inc.cos()];
        let star = position_of(system, light_source(system, world)?, 0.0);
        let to_star = sub(star, position_of(system, world, 0.0));
        // The equinox lies in the orbit's plane.
        let along = sub(to_star, scale(normal, dot(to_star, normal)));
        let x = scale(along, 1.0 / norm(along).max(1e-300));
        let obliquity = b.obliquity_deg.to_radians();
        let z = add(scale(normal, obliquity.cos()), scale(cross(normal, x), obliquity.sin()));
        Some(Self { x, y: cross(z, x), z })
    }

    /// Right ascension and declination of direction `v`, degrees (RA 0–360).
    #[must_use]
    pub fn ra_dec(&self, v: [f64; 3]) -> (f64, f64) {
        let (a, b, c) = (dot(v, self.x), dot(v, self.y), dot(v, self.z));
        (b.atan2(a).to_degrees().rem_euclid(360.0), c.atan2(a.hypot(b)).to_degrees())
    }

    /// The star's declination `t` days after the present, degrees: what drives the seasons
    /// (positive in northern summer).
    #[must_use]
    pub fn star_declination(&self, system: &SolarSystem, world: Key, star: Key, t_days: f64) -> f64 {
        self.ra_dec(sub(position_of(system, star, t_days), position_of(system, world, t_days))).1
    }
}

fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn scale(a: [f64; 3], k: f64) -> [f64; 3] {
    [a[0] * k, a[1] * k, a[2] * k]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

/// Solar (a star covered) or lunar (a moon in the world's shadow).
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum EclipseKind {
    Solar,
    Lunar,
}

/// How much is covered at greatest eclipse.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Coverage {
    Total,
    Annular,
    Partial,
}

#[derive(Debug, Clone, PartialEq, serde::Serialize)]
pub struct Eclipse {
    pub kind: EclipseKind,
    pub coverage: Coverage,
    /// Earth days from the present, at greatest eclipse.
    pub day: f64,
    /// What covers (solar) or is covered (lunar).
    pub body: String,
    /// The star involved.
    pub star: String,
    /// Share of the star's diameter covered where the eclipse is greatest on the world
    /// (solar), or of the moon's diameter inside the umbra (lunar); over 1 when total.
    pub magnitude: f64,
}

/// Eclipses seen from `world` over `days` from `from_day`, by greatest eclipse.
#[must_use]
pub fn eclipses(system: &SolarSystem, world: Key, from_day: f64, days: f64) -> Vec<Eclipse> {
    let Some(node) = system.node(world) else {
        return Vec::new();
    };
    let NodeKind::Body(me) = &node.node else {
        return Vec::new();
    };
    let radius_au = |k: Key| match system.node(k).map(|n| &n.node) {
        Some(NodeKind::Body(b)) => b.radius_km / geoforge_core::body::AU_KM,
        Some(NodeKind::Star(s)) => s.radius_solar * SOLAR_RADIUS_AU,
        _ => 0.0,
    };
    let stars: Vec<&geoforge_core::body::OrbitalNode> = system.nodes.iter().filter(|n| matches!(n.node, NodeKind::Star(_))).collect();
    let moons: Vec<Key> = system.children(world).iter().filter(|n| matches!(n.node, NodeKind::Body(_))).map(|n| n.key).collect();
    // A moon's planet covers the star too.
    let planet = node.parent.filter(|&p| matches!(system.node(p).map(|n| &n.node), Some(NodeKind::Body(_))));
    let occulters: Vec<Key> = moons.iter().copied().chain(planet).collect();
    let my_radius = me.radius_km / geoforge_core::body::AU_KM;
    let mut out = Vec::new();
    for star in &stars {
        for &o in &occulters {
            let is_moon = moons.contains(&o);
            // Step finely against the faster of the occulter's and the world's orbits.
            let period = [period_days(system, o), period_days(system, world)].into_iter().flatten().fold(f64::MAX, f64::min);
            if !period.is_finite() || period <= 0.0 {
                continue;
            }
            let dt = period / 2000.0;
            let geometry = |t: f64, lunar: bool| {
                let (w, b, s) = (position_of(system, world, t), position_of(system, o, t), position_of(system, star.key, t));
                let (to_b, to_s) = (sub(b, w), sub(s, w));
                let (db, ds) = (norm(to_b), norm(to_s));
                let cos = dot(to_b, to_s) / (db * ds).max(1e-300);
                let sep = if lunar { (-cos).clamp(-1.0, 1.0).acos() } else { cos.clamp(-1.0, 1.0).acos() };
                (sep, db, ds)
            };
            for lunar in [false, true] {
                if lunar && !is_moon {
                    continue;
                }
                let steps = (days / dt).ceil() as usize;
                let mut prev = (f64::MAX, f64::MAX);
                for i in 0..=steps + 1 {
                    let t = from_day + i as f64 * dt;
                    let sep = geometry(t, lunar).0;
                    // A minimum one step back.
                    if prev.1 < prev.0 && prev.1 <= sep {
                        let (tm, (sep_m, db, ds)) = refine(|x| geometry(x, lunar), t - 2.0 * dt, t);
                        let (rb, rs) = ((radius_au(o) / db).asin(), (radius_au(star.key) / ds).asin());
                        if let Some(e) = classify(lunar, sep_m, rb, rs, my_radius, db, ds, radius_au(star.key)) {
                            if (from_day..=from_day + days).contains(&tm) {
                                let name = system.node(o).map_or_else(String::new, |n| n.name.clone());
                                out.push(Eclipse { kind: e.0, coverage: e.1, day: tm, body: name, star: star.name.clone(), magnitude: e.2 });
                            }
                        }
                    }
                    prev = (prev.1, sep);
                }
            }
        }
    }
    out.sort_by(|a, b| a.day.total_cmp(&b.day));
    out
}

/// Kind, coverage and magnitude of an eclipse at separation `sep` (radians), if any: `rb`
/// and `rs` the angular radii of the body and the star, `my_radius`, `db` and `ds` the
/// world's radius and the distances to the body and the star (AU), `star_r` the star's
/// radius (AU).
#[allow(clippy::too_many_arguments)]
fn classify(lunar: bool, sep: f64, rb: f64, rs: f64, my_radius: f64, db: f64, ds: f64, star_r: f64) -> Option<(EclipseKind, Coverage, f64)> {
    if lunar {
        // The world's umbra at the moon's distance, as an angle from the world's centre.
        let umbra = ((my_radius - (star_r - my_radius) * db / ds) / db).atan();
        if umbra <= 0.0 || sep >= umbra + rb {
            return None;
        }
        let coverage = if sep + rb <= umbra { Coverage::Total } else { Coverage::Partial };
        return Some((EclipseKind::Lunar, coverage, (umbra + rb - sep) / (2.0 * rb)));
    }
    // Seen from somewhere on the world: the world's own radius shifts the body by up to this.
    let parallax = (my_radius / db).min(1.0).asin();
    if sep >= rb + rs + parallax {
        return None;
    }
    let coverage = if sep <= (rb - rs).abs() + parallax {
        if rb >= rs {
            Coverage::Total
        } else {
            Coverage::Annular
        }
    } else {
        Coverage::Partial
    };
    // Magnitude where the eclipse is greatest on the world: the share of the star's diameter
    // the body covers along the line through both centres (over 1 when total).
    let s = (sep - parallax).max(0.0);
    let magnitude = if rb >= rs { (rs + rb - s) / (2.0 * rs) } else { (rs.min(s + rb) - (-rs).max(s - rb)) / (2.0 * rs) };
    Some((EclipseKind::Solar, coverage, magnitude))
}

/// Golden-section search for the minimum of `f`'s first value on [a, b].
fn refine(f: impl Fn(f64) -> (f64, f64, f64), mut a: f64, mut b: f64) -> (f64, (f64, f64, f64)) {
    let g = 0.5 * (5f64.sqrt() - 1.0);
    for _ in 0..60 {
        let (c, d) = (b - g * (b - a), a + g * (b - a));
        if f(c).0 < f(d).0 {
            b = d;
        } else {
            a = c;
        }
    }
    let t = 0.5 * (a + b);
    (t, f(t))
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// The Sun's radius, AU.
const SOLAR_RADIUS_AU: f64 = 0.004_650_47;
/// Share of a moon's face lit as seen from its planet, with the star at `star`.
#[must_use]
pub fn lit_fraction(moon: [f64; 3], planet: [f64; 3], star: [f64; 3]) -> f64 {
    let to_star = [star[0] - moon[0], star[1] - moon[1], star[2] - moon[2]];
    let to_planet = [planet[0] - moon[0], planet[1] - moon[1], planet[2] - moon[2]];
    let cos_phase = (to_star[0] * to_planet[0] + to_star[1] * to_planet[1] + to_star[2] * to_planet[2]) / (norm(to_star) * norm(to_planet)).max(1e-30);
    0.5 * (1.0 + cos_phase)
}

/// The days, year and months a world keeps.
#[derive(Debug, Clone, PartialEq)]
pub struct Calendar {
    /// Mean solar day, hours; `None` for a world locked to its star.
    pub solar_day_h: Option<f64>,
    /// The year, Earth days.
    pub year_days: f64,
    /// Synodic months (new moon to new moon) by moon, Earth days.
    pub months: Vec<(String, f64)>,
}

impl Calendar {
    /// The calendar of `body`, its moons the tidal partners lighter than it.
    #[must_use]
    pub fn of_body(body: &BodyContext) -> Self {
        let moons: Vec<(String, f64)> = body
            .tides
            .iter()
            .filter(|t| t.mass_kg < body.exterior.mass_kg)
            .map(|t| (t.name.clone(), orbital_period_days(t.distance_km / geoforge_core::body::AU_KM, body.exterior.mass_kg + t.mass_kg)))
            .collect();
        Self::new(body, &moons)
    }

    /// The calendar of `body`, with its moons' orbital periods (Earth days).
    #[must_use]
    pub fn new(body: &BodyContext, moons: &[(String, f64)]) -> Self {
        let clock = Clock::of_body(body);
        let year = clock.year.earth_days();
        let months = moons.iter().map(|(n, p)| (n.clone(), 1.0 / (1.0 / p - 1.0 / year))).collect();
        Self { solar_day_h: clock.day.map(Span::hours), year_days: year, months }
    }

    /// The world's clock: its day and year as spans.
    #[must_use]
    pub fn clock(&self) -> Clock {
        Clock::new(self.solar_day_h.map(Span::from_hours), Span::from_earth_days(self.year_days))
    }

    /// Local days in a year, if the world has days.
    #[must_use]
    pub fn days_per_year(&self) -> Option<f64> {
        self.solar_day_h.map(|d| self.year_days * 24.0 / d)
    }

    /// The date `t` Earth days after the present: the year (0 the present one) and the day of
    /// that year, in local days and their fraction (Earth days if the world has no days).
    #[must_use]
    pub fn date(&self, t_days: f64) -> (i64, f64) {
        let years = (t_days / self.year_days).floor();
        let into = t_days - years * self.year_days;
        (years as i64, self.solar_day_h.map_or(into, |d| into * 24.0 / d))
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn julian_dates_read_as_calendar_dates() {
        assert_eq!(utc_date(2_451_545.0 + TT_MINUS_UTC_DAYS), "2000-01-01 12:00 UTC");
        assert_eq!(utc_date(2_461_120.116_08), "2026-03-20 14:46 UTC");
    }

    use super::*;
    use geoforge_core::{Archetype, Key};

    fn days(a: Archetype) -> Calendar {
        Calendar::new(&a.context(Key(1), "x"), &[("Moon".into(), 27.3217)])
    }

    #[test]
    fn solar_days_and_months() {
        let earth = days(Archetype::Earth);
        assert!((earth.solar_day_h.unwrap() - 24.0).abs() < 0.01 && (earth.months[0].1 - 29.53).abs() < 0.01, "{earth:?}");
        assert!((earth.days_per_year().unwrap() - 365.25).abs() < 0.02);
        assert!((days(Archetype::Mars).solar_day_h.unwrap() - 24.66).abs() < 0.01);
        assert!((days(Archetype::Venus).solar_day_h.unwrap() / 24.0 - 116.75).abs() < 0.1, "retrograde");
        assert!((days(Archetype::Mercury).solar_day_h.unwrap() / 24.0 - 175.9).abs() < 0.3);
        let (y, d) = earth.date(400.0);
        assert!(y == 1 && (d - 34.7).abs() < 0.1, "{y} {d}");
    }

    #[test]
    fn kepler_orbits() {
        let e = Elements { a_au: 2.0, e: 0.5, i_deg: 30.0, lan_deg: 40.0, argp_deg: 50.0, m0_deg: 0.0 };
        // At periapsis now, apoapsis half a period on, back after one.
        let r = |t: f64| norm(relative_position(&e, 100.0, t));
        assert!((r(0.0) - 1.0).abs() < 1e-9 && (r(50.0) - 3.0).abs() < 1e-9 && (r(100.0) - 1.0).abs() < 1e-9);
        // Equal areas: the planet is near apoapsis far longer than near periapsis.
        let near = (0..100).filter(|&t| r(f64::from(t)) < 2.0).count();
        assert!(near < 40, "{near} of 100 days inside the mean distance");
        // Full moon: the moon beyond its planet from the star is fully lit.
        assert!((lit_fraction([2.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0; 3]) - 1.0).abs() < 1e-12);
        assert!(lit_fraction([0.5, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0; 3]) < 1e-12);
    }

    /// Sun, Earth and Moon on their real orbits: about 2.4 solar and 1.5 umbral lunar
    /// eclipses a year (NASA's five-millennium catalogues: 2.38 and 1.45), both total and
    /// annular solar eclipses, and none outside the eclipse seasons.
    #[test]
    fn earths_eclipses() {
        let earth = Archetype::Earth.context(Key(1), "Earth").exterior;
        let mut moon = earth;
        moon.mass_kg = 7.342e22;
        moon.radius_km = 1737.4;
        let sys: SolarSystem = serde_json::from_value(serde_json::json!({
            "root": { "key": 1, "seed": 1, "name": "Sol", "epoch": { "system_age_gyr": 4.57 } },
            "nodes": [
                { "key": 2, "name": "Sun", "node": { "kind": "star", "mass_solar": 1.0, "luminosity_solar": 1.0, "radius_solar": 1.0, "temperature_k": 5772.0, "age_gyr": 4.57, "metallicity_dex": 0.0, "spectral_class": "G2V" } },
                { "key": 3, "name": "Earth", "parent": 2, "node": serde_json::to_value(NodeKind::Body(earth)).unwrap(),
                  "elements": { "a_au": 1.000001, "e": 0.0167, "i_deg": 0.0, "lan_deg": 0.0, "argp_deg": 102.9, "m0_deg": 0.0 } },
                { "key": 4, "name": "Moon", "parent": 3, "node": serde_json::to_value(NodeKind::Body(moon)).unwrap(),
                  "elements": { "a_au": 0.002_569_6, "e": 0.0549, "i_deg": 5.145, "lan_deg": 125.0, "argp_deg": 318.0, "m0_deg": 115.0 } }
            ]
        }))
        .expect("a system");
        let years = 30.0;
        let found = eclipses(&sys, Key(3), 0.0, years * 365.25);
        let count = |k: EclipseKind| found.iter().filter(|e| e.kind == k).count() as f64 / years;
        let (solar, lunar) = (count(EclipseKind::Solar), count(EclipseKind::Lunar));
        assert!((1.8..3.2).contains(&solar), "{solar} solar eclipses a year");
        assert!((0.9..2.2).contains(&lunar), "{lunar} lunar eclipses a year");
        for c in [Coverage::Total, Coverage::Annular] {
            assert!(found.iter().any(|e| e.kind == EclipseKind::Solar && e.coverage == c), "no {c:?} solar eclipse");
        }
        assert!(found.iter().any(|e| e.kind == EclipseKind::Lunar && e.coverage == Coverage::Total));
        // The present is the March equinox: the Sun at RA 0, declination 0, then +23.4° at
        // the June solstice (93 days on) and −23.4° in December.
        let eq = Equatorial::of(&sys, Key(3)).expect("a frame");
        let dec = |t: f64| eq.star_declination(&sys, Key(3), Key(2), t);
        assert!(dec(0.0).abs() < 1e-9 && dec(5.0) > 1.0, "{}", dec(5.0));
        let june = (60..130).map(|d| (d, dec(f64::from(d)))).max_by(|a, b| a.1.total_cmp(&b.1)).unwrap();
        assert!((88..98).contains(&june.0) && (june.1 - 23.44).abs() < 0.05, "{june:?}");
        assert!((dec(276.0) + 23.44).abs() < 0.2, "{}", dec(276.0));
        for e in found.iter().filter(|e| e.kind == EclipseKind::Solar) {
            let ok = match e.coverage {
                Coverage::Total => e.magnitude >= 1.0,
                Coverage::Annular => (0.9..1.0).contains(&e.magnitude),
                Coverage::Partial => (0.0..1.0).contains(&e.magnitude),
            };
            assert!(ok, "{e:?}");
        }
    }
    #[test]
    fn a_system_moves() {
        let sys: SolarSystem = serde_json::from_str(include_str!("../tests/fixtures/system_binary.json")).expect("fixture");
        let now = positions(&sys, 0.0);
        for p in &now {
            if let Some(e) = sys.node(p.key).and_then(|n| n.elements) {
                assert!(p.distance_au >= e.a_au * (1.0 - e.e) - 1e-9 && p.distance_au <= e.a_au * (1.0 + e.e) + 1e-9, "{}", p.name);
                let later = positions(&sys, p.period_days);
                let back = later.iter().find(|q| q.key == p.key).unwrap();
                assert!((back.distance_au - p.distance_au).abs() < 1e-6, "{} returns after a period", p.name);
            }
        }
        // The pair turns about its barycentre: the same period, opposite at every moment,
        // their mass-weighted positions cancelling.
        let stars: Vec<&geoforge_core::body::OrbitalNode> = sys.nodes.iter().filter(|n| matches!(n.node, NodeKind::Star(_))).collect();
        assert_eq!(stars.len(), 2);
        let (pa, pb) = (period_days(&sys, stars[0].key).unwrap(), period_days(&sys, stars[1].key).unwrap());
        assert!((pa / pb - 1.0).abs() < 1e-9, "{pa} vs {pb}");
        for t in [0.0, 3.3, 17.0, 40.1] {
            let at = positions(&sys, t);
            let pos = |k: Key| at.iter().find(|p| p.key == k).unwrap().position_au;
            let (ma, mb) = (sys.mass_of(stars[0].key), sys.mass_of(stars[1].key));
            let (ra, rb) = (pos(stars[0].key), pos(stars[1].key));
            for i in 0..3 {
                assert!((ma * ra[i] + mb * rb[i]).abs() / (ma + mb) < 1e-9, "t {t}: {ra:?} {rb:?}");
            }
        }
    }
}
