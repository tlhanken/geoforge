//! A planet's moons: regular ones born in its disk, small inner ones near its rings, big ones
//! from giant impacts, and captured ones.
//!
//! - **Regular moons** of giants, born in the planet's own disk, hold a share of its mass that
//!   the disk's inflow sets: ~2×10⁻⁴ for a gas giant, ~1×10⁻⁴ for an ice giant (Canup & Ward
//!   2006). Gas giants have 3–6 (Jupiter's four, Saturn's seven round ones), ice giants 3–5
//!   (Uranus's five). How the mass is shared is drawn over two and a half decades, so one moon
//!   may hold nearly all of it (Titan, 96% of Saturn's) or several share it (the Galileans).
//!   They lie on a widening ladder from ~6 planet radii, each 1.5–2 times the last.
//! - **Small inner moons** inside them, 1–3 shepherds and fragments from 1.8 to 3.5 radii
//!   (Amalthea, Saturn's ring moons), 10⁻⁹–10⁻⁸ of the planet's mass.
//! - **A captured giant moon**: an ice giant takes one a third of the time (Neptune's
//!   Triton, from a binary that passed too close: Agnor & Hamilton 2006). It orbits backward,
//!   tilted, ~14 planet radii out, and its capture scatters the regular moons: one small one
//!   is left inside it.
//! - **Irregular moons** of every giant, captured from passing small bodies onto far, tilted,
//!   often backward orbits between a tenth and half of its Hill sphere: a population (tens to
//!   hundreds over 1 km, up to ~100–200 km: Himalia, Phoebe), not listed one by one.
//! - **A rocky planet's big moon**, from a giant impact (Earth's): `big_moon_chance` of
//!   worlds of 0.3 Earth masses or more that are not tidally locked; 0.4–2% of the planet's
//!   mass, 30–80 radii out.
//! - **Captured specks** of a rocky planet beside an asteroid belt: half of those whose orbit
//!   lies within 1.5 times of the belt's inner edge, unlocked, take one or two (Mars's Phobos
//!   and Deimos, ~10¹⁶ kg), 2.5–7 radii out.
//!
//! - **Rings** ([`rings`]): faint dust on every giant, bright ice on some cold ones.
//!
//! Every moon lies well inside the planet's Hill sphere (the caller keeps those inside 0.4 of
//! it). Moons of dwarf planets are not drawn.

use crate::draw::Draw;
use crate::system::{big_moon_orbit_km, moon, roche_km, MoonSpec};
use geoforge_core::body::{BodyClass, BodyExterior, Composition, NodeKind, OrbitalNode, Population, PowerLaw, Species, AU_KM, EARTH_MASS_KG};
use geoforge_core::key::{Key, KeyAllocator};

/// The planet moons are drawn for.
pub(crate) struct Host<'a> {
    pub key: Key,
    pub name: &'a str,
    pub ext: &'a BodyExterior,
    /// Its Hill radius, km.
    pub hill_km: f64,
    /// Tidally locked to its star (a rocky planet then keeps no big moon).
    pub locked: bool,
    /// Its orbit against the inner edge of the system's asteroid belt (a/inner), if any.
    pub belt_ratio: Option<f64>,
}

/// Share of a rocky planet's big moons, by the generator calling.
pub(crate) const BIG_MOON_CHANCE: f64 = 0.2;
/// Share of ice giants with a captured giant moon.
const CAPTURE_CHANCE: f64 = 1.0 / 3.0;
const LETTERS: &str = "abcdefghijklmnopqrstuvwxyz";

/// The moons of `host`, nearest first and lettered so, and its irregular moons as a
/// population (giants only).
pub(crate) fn moons(d: &mut Draw, keys: &mut KeyAllocator, host: &Host<'_>, big_moon_chance: f64) -> (Vec<OrbitalNode>, Option<OrbitalNode>) {
    let p = host.ext;
    let mut specs: Vec<(MoonSpec, Orbit)> = Vec::new();
    let regular = |d: &mut Draw, count: usize, share: f64, specs: &mut Vec<(MoonSpec, Orbit)>| {
        let weights: Vec<f64> = (0..count).map(|_| 10f64.powf(d.range(-2.5, 0.0))).collect();
        let total: f64 = weights.iter().sum();
        let budget = p.mass_kg * share * d.range(0.5, 1.5);
        let mut a_km = p.radius_km * d.range(5.5, 7.0);
        for w in weights {
            specs.push((MoonSpec { name: String::new(), mass_kg: budget * w / total, a_km, density_g_cc: 1.9 }, Orbit::Regular));
            a_km *= d.range(1.5, 2.0);
        }
    };
    let inner = |d: &mut Draw, count: usize, specs: &mut Vec<(MoonSpec, Orbit)>| {
        let floor = roche_km(p, 1.0).max(1.8 * p.radius_km);
        for _ in 0..count {
            let a_km = (p.radius_km * d.range(1.8, 3.5)).max(floor * 1.05);
            specs.push((MoonSpec { name: String::new(), mass_kg: p.mass_kg * 10f64.powf(d.range(-9.0, -8.0)), a_km, density_g_cc: 1.0 }, Orbit::Regular));
        }
    };
    match p.class {
        BodyClass::GasGiant => {
            let count = d.int(3, 6) as usize;
            regular(d, count, 2e-4, &mut specs);
            let n = d.int(1, 3) as usize;
            inner(d, n, &mut specs);
        }
        BodyClass::IceGiant => {
            if d.u() < CAPTURE_CHANCE {
                let a_km = p.radius_km * d.range(12.0, 16.0);
                specs.push((MoonSpec { name: String::new(), mass_kg: p.mass_kg * d.range(1e-4, 3e-4), a_km, density_g_cc: 2.06 }, Orbit::Captured));
                inner(d, 1, &mut specs);
            } else {
                let count = d.int(3, 5) as usize;
                regular(d, count, 1e-4, &mut specs);
                let n = d.int(1, 2) as usize;
                inner(d, n, &mut specs);
            }
        }
        BodyClass::Terrestrial => {
            if !host.locked && p.mass_kg > 0.3 * EARTH_MASS_KG && d.u() < big_moon_chance {
                let mass_kg = p.mass_kg * d.range(0.004, 0.02);
                let drawn = d.range(30.0, 80.0) * p.radius_km;
                if let Some(a_km) = big_moon_orbit_km(drawn, p, host.hill_km) {
                    specs.push((MoonSpec { name: String::new(), mass_kg, a_km, density_g_cc: 3.3 }, Orbit::Regular));
                }
            }
            if !host.locked && host.belt_ratio.is_some_and(|r| (1.0 / 1.5..1.0).contains(&r)) && d.u() < 0.5 {
                for _ in 0..d.int(1, 2) {
                    let a_km = (p.radius_km * d.range(2.5, 7.0)).max(1.1 * roche_km(p, 1.9));
                    specs.push((MoonSpec { name: String::new(), mass_kg: 10f64.powf(d.range(15.0, 16.5)), a_km, density_g_cc: 1.9 }, Orbit::Captured));
                }
            }
        }
        BodyClass::SmallBody => {}
    }
    specs.sort_by(|a, b| a.0.a_km.total_cmp(&b.0.a_km));
    let nodes = specs
        .into_iter()
        .enumerate()
        .map(|(i, (mut spec, orbit))| {
            spec.name = format!("{} {}", host.name, LETTERS.get(i..=i).unwrap_or("+"));
            let mut node = moon(d, keys.fresh(), host.key, p.mass_kg, spec);
            if orbit == Orbit::Captured {
                if let Some(e) = node.elements.as_mut() {
                    // Captured: a tilted orbit, backward for a giant's big catch.
                    e.i_deg = if p.class == BodyClass::Terrestrial { d.range(0.0, 3.0) } else { d.range(130.0, 170.0) };
                }
            }
            node
        })
        .collect();
    let irregular = matches!(p.class, BodyClass::GasGiant | BodyClass::IceGiant).then(|| {
        let count = if p.class == BodyClass::GasGiant { 10f64.powf(d.range(1.5, 2.5)) } else { 10f64.powf(d.range(1.0, 2.0)) };
        OrbitalNode {
            key: keys.fresh(),
            name: format!("{} irregular moons", host.name),
            parent: Some(host.key),
            node: NodeKind::Belt(Population {
                inner_au: 0.1 * host.hill_km / AU_KM,
                outer_au: 0.5 * host.hill_km / AU_KM,
                count_estimate: count,
                size_law: PowerLaw { min_km: 1.0, max_km: d.range(60.0, 200.0), exponent: 3.5, knee: None },
                composition: Composition::from_weights(&[(Species::Silicate, 0.5), (Species::Ice, 0.5)]),
            }),
            elements: None,
            precession: None,
        }
    });
    (nodes, irregular)
}

/// Coldest a giant can be (equilibrium, K) and still keep bright ice rings: water ice
/// sublimes away nearer its star (Saturn's 81 K and Jupiter's 110 K both keep ice; Jupiter's
/// rings are dust for want of a ring-making moon, not of cold).
const ICE_RING_MAX_K: f64 = 150.0;
/// Share of cold giants whose rings are bright: Saturn's, from a moon torn apart or a disk
/// of debris, may be young and passing (Iess et al. 2019), so one giant in four.
const BRIGHT_RING_CHANCE: f64 = 0.25;

/// A giant's rings, inside its Roche limit: faint dust on every giant (Jupiter's, Uranus's
/// and Neptune's, optical depth 10⁻⁶–10⁻⁵, 1–100 µm grains, from 1.3 to 2.5 radii); on a
/// quarter of those cold enough for ice (`equilibrium_k` under 150 K), bright ice like
/// Saturn's instead (depth 0.2–1.6, 1 cm–10 m blocks, from 1.2 radii to the Roche limit for
/// ice). `None` for other bodies. Drawn from its own stream: rings change no other draw.
pub(crate) fn rings(d: &mut Draw, key: Key, host: &Host<'_>, equilibrium_k: f64) -> Option<OrbitalNode> {
    let p = host.ext;
    if !matches!(p.class, BodyClass::GasGiant | BodyClass::IceGiant) {
        return None;
    }
    let bright = equilibrium_k < ICE_RING_MAX_K && d.u() < BRIGHT_RING_CHANCE;
    let (inner_km, outer_km, depth, size_law, composition) = if bright {
        let outer = roche_km(p, 0.9).max(1.6 * p.radius_km) * d.range(0.85, 1.0);
        (
            p.radius_km * d.range(1.15, 1.3),
            outer,
            10f64.powf(d.range(-0.7, 0.2)),
            PowerLaw { min_km: 1e-5, max_km: 0.01, exponent: 3.0, knee: None },
            Composition::from_weights(&[(Species::Ice, 0.95), (Species::Silicate, 0.05)]),
        )
    } else {
        (
            p.radius_km * d.range(1.25, 1.4),
            p.radius_km * d.range(1.8, 3.2),
            10f64.powf(d.range(-6.5, -5.0)),
            PowerLaw { min_km: 1e-9, max_km: 1e-7, exponent: 3.5, knee: None },
            Composition::from_weights(&[(Species::Silicate, 0.8), (Species::Ice, 0.2)]),
        )
    };
    let mut pop = Population { inner_au: inner_km / AU_KM, outer_au: outer_km / AU_KM, count_estimate: 1.0, size_law, composition };
    pop.count_estimate = depth / pop.optical_depth();
    Some(OrbitalNode { key, name: format!("{}'s rings", host.name), parent: Some(host.key), node: NodeKind::Ring(pop), elements: None, precession: None })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Orbit {
    Regular,
    Captured,
}

#[cfg(test)]
mod tests {
    use super::*;
    const JUPITER_MASS_KG: f64 = 1.898e27;
    use geoforge_core::seed::Seed;

    fn giant(class: BodyClass, mass_kg: f64, radius_km: f64) -> BodyExterior {
        BodyExterior { class, mass_kg, radius_km, rotation_h: 10.0, obliquity_deg: 3.0 }
    }

    fn draw(ext: &BodyExterior, s: u64, belt_ratio: Option<f64>) -> (Vec<OrbitalNode>, Option<OrbitalNode>) {
        let mut keys = KeyAllocator::starting_at(10);
        let host = Host { key: Key(1), name: "P", ext, hill_km: 5e7, locked: false, belt_ratio };
        moons(&mut Draw::new(Seed(s)), &mut keys, &host, BIG_MOON_CHANCE)
    }

    fn masses(nodes: &[OrbitalNode]) -> Vec<f64> {
        nodes.iter().filter_map(|n| if let NodeKind::Body(b) = &n.node { Some(b.mass_kg) } else { None }).collect()
    }

    /// Jupiter-like giants: 4–9 moons in all, their regular ones holding ~2×10⁻⁴ of the
    /// planet (the Galileans: 2.1×10⁻⁴), sometimes one dominant (Titan's 96%), sometimes
    /// shared; lettered outward; irregulars as a population.
    #[test]
    fn gas_giants_have_moon_systems() {
        let jupiter = giant(BodyClass::GasGiant, JUPITER_MASS_KG, 69_911.0);
        let (mut dominant, mut shared) = (0, 0);
        for s in 0..200 {
            let (nodes, irregular) = draw(&jupiter, s, None);
            assert!((4..=9).contains(&nodes.len()), "{}", nodes.len());
            let m = masses(&nodes);
            let total: f64 = m.iter().sum();
            assert!((0.9e-4..3.1e-4).contains(&(total / JUPITER_MASS_KG)), "{}", total / JUPITER_MASS_KG);
            let top = m.iter().copied().fold(0.0, f64::max) / total;
            dominant += usize::from(top > 0.9);
            shared += usize::from(top < 0.6);
            let a: Vec<f64> = nodes.iter().map(|n| n.elements.unwrap().a_au).collect();
            assert!(a.windows(2).all(|w| w[0] <= w[1]) && nodes[0].name == "P a");
            assert!(nodes.iter().all(|n| n.elements.unwrap().a_au * AU_KM > roche_km(&jupiter, 1.0)));
            let Some(OrbitalNode { node: NodeKind::Belt(pop), .. }) = irregular else { panic!("irregulars") };
            assert!(pop.count_estimate >= 30.0);
        }
        assert!(dominant > 10 && shared > 10, "{dominant} dominant, {shared} shared");
    }

    /// A third of ice giants hold a captured moon on a backward orbit with one small moon
    /// inside it (Neptune); the rest 3–5 regular moons and 1–2 small ones (Uranus).
    #[test]
    fn ice_giants_regular_or_captured() {
        let neptune = giant(BodyClass::IceGiant, 1.02e26, 24_622.0);
        let mut captured = 0;
        for s in 0..300 {
            let (nodes, _) = draw(&neptune, s, None);
            if nodes.iter().any(|n| n.elements.unwrap().i_deg > 90.0) {
                captured += 1;
                assert_eq!(nodes.len(), 2);
            } else {
                assert!((4..=7).contains(&nodes.len()));
            }
        }
        assert!((70..130).contains(&captured), "{captured} of 300");
    }

    /// Every giant has rings: faint dust, or bright ice on a quarter of cold ones, never on
    /// warm ones; bright rings lie inside the Roche limit for ice, as dense as Saturn's.
    #[test]
    fn giants_have_rings() {
        let saturn = giant(BodyClass::GasGiant, 5.68e26, 58_232.0);
        let host = Host { key: Key(1), name: "S", ext: &saturn, hill_km: 6e7, locked: false, belt_ratio: None };
        let depth = |n: &OrbitalNode| if let NodeKind::Ring(p) = &n.node { (p.optical_depth(), p.inner_au * AU_KM, p.outer_au * AU_KM) } else { panic!("a ring") };
        let mut bright = 0;
        for s in 0..400 {
            let r = rings(&mut Draw::new(Seed(s)), Key(9), &host, 81.0).expect("rings");
            let (tau, inner, outer) = depth(&r);
            assert!(inner > saturn.radius_km && outer > inner, "{inner}–{outer}");
            if tau > 0.1 {
                bright += 1;
                assert!(tau < 1.7 && outer <= roche_km(&saturn, 0.9) * 1.0001, "{tau} to {outer}");
            } else {
                assert!(tau < 1e-4);
            }
            let warm = rings(&mut Draw::new(Seed(s)), Key(9), &host, 300.0).unwrap();
            assert!(depth(&warm).0 < 1e-4, "no ice rings near the star");
        }
        assert!((70..130).contains(&bright), "{bright} of 400");
        let earth = BodyExterior { class: BodyClass::Terrestrial, mass_kg: EARTH_MASS_KG, radius_km: 6371.0, rotation_h: 24.0, obliquity_deg: 23.0 };
        assert!(rings(&mut Draw::new(Seed(1)), Key(9), &Host { ext: &earth, ..host }, 255.0).is_none());
    }

    /// A rocky planet just inside an asteroid belt may hold captured specks (Mars); far from
    /// one, only a big moon now and then (Earth).
    #[test]
    fn rocky_planets_by_a_belt_catch_specks() {
        let mars = BodyExterior { class: BodyClass::Terrestrial, mass_kg: 6.4e23, radius_km: 3390.0, rotation_h: 24.6, obliquity_deg: 25.0 };
        let near: usize = (0..200).map(|s| draw(&mars, s, Some(0.7)).0.len()).sum();
        let far: usize = (0..200).map(|s| draw(&mars, s, Some(0.3)).0.len()).sum();
        assert!(near > 100 && far == 0, "{near} near, {far} far");
        let earth = BodyExterior { mass_kg: EARTH_MASS_KG, radius_km: 6371.0, ..mars };
        let big: usize = (0..500).map(|s| draw(&earth, s, None).0.len()).sum();
        assert!((60..140).contains(&big), "{big} of 500");
    }
}
