//! Where everything in a system is at a moment, and the calendar a world keeps.
//!
//! Each node moves on a fixed Kepler orbit around its parent from its elements, which hold
//! at the system's present (t = 0); time is in Earth days. Around a barycentre a star's
//! elements describe its own orbit about it, so its period is the pair's relative orbit's.
//! Planets do not perturb each other here.
//!
//! A world's solar day follows from its spin and its year, 1/day = 1/spin − 1/year (plus for
//! a retrograde spinner such as Venus, whose tilt is over 90°); a world locked to its star has
//! none. A month is a moon's synodic period, new moon to new moon: 1/month = 1/orbit − 1/year.

use crate::body::{NodeKind, SolarSystem};
use crate::derived::orbital_period_days;
use crate::key::Key;
use crate::{BodyContext, Elements};
use std::collections::BTreeMap;
use std::f64::consts::TAU;

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
    let is_star = |n: &crate::body::OrbitalNode| matches!(n.node, NodeKind::Star(_));
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
    let mean = (e.m0_deg.to_radians() + TAU * t_days / period_days).rem_euclid(TAU);
    let ecc = e.e.clamp(0.0, 0.99);
    // Kepler's equation, E − e·sin E = M, by Newton's method.
    let mut big_e = if ecc > 0.8 { std::f64::consts::PI } else { mean };
    for _ in 0..50 {
        let step = (big_e - ecc * big_e.sin() - mean) / (1.0 - ecc * big_e.cos());
        big_e -= step;
        if step.abs() < 1e-13 {
            break;
        }
    }
    let nu = 2.0 * ((1.0 + ecc).sqrt() * (big_e / 2.0).sin()).atan2((1.0 - ecc).sqrt() * (big_e / 2.0).cos());
    let r = e.a_au * (1.0 - ecc * big_e.cos());
    let u = e.argp_deg.to_radians() + nu;
    let (node, inc) = (e.lan_deg.to_radians(), e.i_deg.to_radians());
    [
        r * (node.cos() * u.cos() - node.sin() * u.sin() * inc.cos()),
        r * (node.sin() * u.cos() + node.cos() * u.sin() * inc.cos()),
        r * u.sin() * inc.sin(),
    ]
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
        let (rel, period) = match (n.elements, period_days(system, n.key)) {
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
    let pos = match (node.parent, node.elements, period_days(system, key)) {
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

fn norm(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

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
    /// The calendar of `body`, with its moons' orbital periods (Earth days).
    #[must_use]
    pub fn new(body: &BodyContext, moons: &[(String, f64)]) -> Self {
        let (spin_d, year) = (body.exterior.rotation_h.abs() / 24.0, body.irradiance.year_days);
        let retrograde = body.exterior.obliquity_deg > 90.0;
        let solar_day_h = (!body.is_tidally_locked()).then(|| 24.0 / (1.0 / spin_d + if retrograde { 1.0 } else { -1.0 } / year));
        let months = moons.iter().map(|(n, p)| (n.clone(), 1.0 / (1.0 / p - 1.0 / year))).collect();
        Self { solar_day_h, year_days: year, months }
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
    use super::*;
    use crate::{Archetype, Key};

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
        let stars: Vec<&crate::body::OrbitalNode> = sys.nodes.iter().filter(|n| matches!(n.node, NodeKind::Star(_))).collect();
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
