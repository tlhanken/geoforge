//! Systems drawn from a galaxy: any star, living or dead, alone or in a multiple, with the
//! planets it formed and kept. The focus is the most Earth-like solid body, if there is one.
//!
//! The tree: a single star at the root; a wide binary as the primary with its companion in
//! orbit (planets around the primary, inside the stable zone); a close binary (under 1 AU at
//! birth) as a barycentre with both stars and the planets around both; a triple's outer star
//! orbits whichever of these is the root. A supernova that blows off more than half a pair's
//! mass unbinds it, and the companion is gone.
//!
//! Planets close to their star are tidally locked: within 0.35 AU·M^(1/3)·(age/4.5 Gyr)^(1/6)
//! (Mercury's distance for the Sun; the whole habitable zone of a red dwarf). An eccentric one
//! settles into Mercury's 3:2 resonance instead.

use crate::context::{belt_member, body_context};
use crate::draw::Draw;
use crate::galaxy::Galaxy;
use crate::planets::{self, binary_survival, p_type_limit, s_type_limit, Architecture, Belt, Disk};
use crate::population::{draw as draw_phenotype, Companion, Phenotype};
use crate::stars::evolve;
use crate::system::{big_moon_orbit_km, hill_km, moon, moon_system, planet_radius_earth, Generated, MoonSpec, MOON_HILL_SHARE, ROMAN};
use geoforge_core::body::{
    BodyClass, BodyContext, BodyExterior, Composition, Elements, Epoch, NodeKind, OrbitalNode, Population, PowerLaw, SolarSystem, Species, Star, StellarStage, SystemRoot, TectonicRegime, AU_KM,
    EARTH_MASS_KG, EARTH_RADIUS_KM, SOLAR_MASS_KG,
};
use geoforge_core::derived::{main_sequence_luminosity, orbital_period_days};
use geoforge_core::key::{Key, KeyAllocator};
use geoforge_core::seed::Seed;
use std::fmt;

/// A system with nothing solid in it to make a world of.
#[derive(Debug, Clone, PartialEq)]
pub struct Barren {
    pub system: Box<SolarSystem>,
    /// What it is, in a sentence.
    pub reason: String,
}

impl fmt::Display for Barren {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.reason)
    }
}

impl std::error::Error for Barren {}

/// Binaries closer than this at birth (AU) are orbited as a pair.
const CLOSE_PAIR_AU: f64 = 1.0;

/// A system drawn from `galaxy`, or [`Barren`] when it has no solid body (a black hole, say).
///
/// # Errors
/// When the system has no solid body to focus on.
pub fn galactic_system(seed: Seed, galaxy: &Galaxy) -> Result<Generated, Barren> {
    build(seed, &draw_phenotype(galaxy, seed))
}

/// The system of phenotype `ph`.
///
/// # Errors
/// When the system has no solid body to focus on.
pub fn build(seed: Seed, ph: &Phenotype) -> Result<Generated, Barren> {
    let mut d = Draw::new(seed.child(b"galactic-system"));
    let mut keys = KeyAllocator::starting_at(0);
    let root_key = keys.fresh();
    let sys_name = geoforge_names::name(seed.0);
    let (age, feh) = (ph.age_gyr, ph.metallicity_dex);
    let primary = evolve(ph.primary_mass, age, feh);
    let companion = ph.companions.first().map(|c| (evolve(c.initial_mass, age, feh), *c));
    let outer = ph.companions.get(1).map(|c| (evolve(c.initial_mass, age, feh), *c));
    let lettered = |l: &str| if ph.companions.is_empty() { sys_name.clone() } else { format!("{sys_name} {l}") };

    // The inner pair, if it is still bound.
    let pair_birth = ph.primary_mass + ph.companions.first().map_or(0.0, |c| c.initial_mass);
    let pair_now = primary.mass_solar + companion.as_ref().map_or(0.0, |c| c.0.mass_solar);
    let exploded = |s: &Star| matches!(s.stage, StellarStage::NeutronStar | StellarStage::BlackHole);
    let bound = |stars: &[&Star], birth: f64, now: f64| !(stars.iter().any(|s| exploded(s)) && birth - now > 0.5 * birth);
    let companion = companion.filter(|(s, _)| bound(&[&primary, s], pair_birth, pair_now));
    let close = companion.as_ref().is_some_and(|(_, c)| c.a_au < CLOSE_PAIR_AU);

    let mut nodes = Vec::new();
    let a_key = keys.fresh();
    let host_key;
    if let (true, Some((b, c))) = (close, &companion) {
        // A barycentre with both stars around it.
        host_key = keys.fresh();
        nodes.push(OrbitalNode { key: host_key, name: sys_name.clone(), parent: None, node: NodeKind::Barycenter, elements: None });
        let a = c.a_au * pair_birth / pair_now;
        let share_b = b.mass_solar / pair_now;
        // One orbit seen from the barycentre: B mirrors A, opposite it at every moment.
        let orbit_a = orbit(&mut d, a * share_b, c.e);
        let orbit_b = Elements { a_au: a * (1.0 - share_b), argp_deg: (orbit_a.argp_deg + 180.0) % 360.0, ..orbit_a };
        nodes.push(star_node(a_key, lettered("A"), Some(host_key), primary.clone(), Some(orbit_a)));
        nodes.push(star_node(keys.fresh(), lettered("B"), Some(host_key), b.clone(), Some(orbit_b)));
    } else {
        host_key = a_key;
        nodes.push(star_node(a_key, lettered("A"), None, primary.clone(), None));
        if let Some((b, c)) = &companion {
            nodes.push(star_node(keys.fresh(), lettered("B"), Some(a_key), b.clone(), Some(orbit(&mut d, c.a_au * pair_birth / pair_now, c.e))));
        }
    }
    let total_birth = pair_birth + ph.companions.get(1).map_or(0.0, |c| c.initial_mass);
    let total_now = pair_now + outer.as_ref().map_or(0.0, |o| o.0.mass_solar);
    let outer = outer.filter(|(s, _)| companion.is_some() && bound(&[&primary, s], total_birth, total_now));
    if let Some((c_star, c)) = &outer {
        nodes.push(star_node(keys.fresh(), lettered("C"), Some(host_key), c_star.clone(), Some(orbit(&mut d, c.a_au * total_birth / total_now, c.e))));
    }

    // Planets form around the host as it was born, then live through what its star became.
    let disk = birth_disk(ph, close);
    let host_now = if close { pair_now } else { primary.mass_solar };
    let host_birth = if close { pair_birth } else { ph.primary_mass };
    let formed = planets::form(&mut d, &disk);
    let arch = planets::fate(&mut d, formed, &primary, host_now, host_birth);

    let mut planets = Vec::new();
    for (i, p) in arch.planets.iter().enumerate() {
        let key = keys.fresh();
        let nm = format!("{sys_name} {}", ROMAN.get(i).copied().unwrap_or("+"));
        let ext = planet(&mut d, p, host_now, age);
        let el = Elements { a_au: p.a_au, e: p.e, i_deg: d.range(0.0, 3.5), lan_deg: d.range(0.0, 360.0), argp_deg: d.range(0.0, 360.0), m0_deg: d.range(0.0, 360.0) };
        nodes.push(OrbitalNode { key, name: nm.clone(), parent: Some(host_key), node: NodeKind::Body(ext), elements: Some(el) });
        planets.push((key, nm, ext, p.a_au));
    }
    nodes.extend(belts(&mut d, &mut keys, host_key, &sys_name, &arch));
    for (key, nm, ext, a_au) in &planets {
        let hill = hill_km(*a_au, ext.mass_kg, host_now);
        let locked = *a_au < tidal_lock_au(host_now, age);
        let moons = match ext.class {
            BodyClass::GasGiant => {
                let count = 2 + (d.u() * 4.0) as usize;
                moon_system(&mut d, &mut keys, *key, nm, ext, count, 2e-4)
            }
            BodyClass::IceGiant => {
                let count = (d.u() * 3.0) as usize;
                moon_system(&mut d, &mut keys, *key, nm, ext, count, 1e-4)
            }
            // A big moon from a giant impact, as Earth's: uncommon.
            BodyClass::Terrestrial if !locked && ext.mass_kg > 0.3 * EARTH_MASS_KG && d.u() < 0.2 => {
                let mass_kg = ext.mass_kg * d.range(0.004, 0.02);
                let drawn = d.range(30.0, 80.0) * ext.radius_km;
                big_moon_orbit_km(drawn, ext, hill)
                    .map(|a_km| moon(&mut d, keys.fresh(), *key, ext.mass_kg, MoonSpec { name: format!("{nm} a"), mass_kg, a_km, density_g_cc: 3.3 }))
                    .into_iter()
                    .collect()
            }
            _ => Vec::new(),
        };
        // Moons survive only well inside the planet's Hill sphere.
        nodes.extend(moons.into_iter().filter(|m| m.elements.is_some_and(|e| e.a_au * AU_KM < MOON_HILL_SHARE * hill)));
    }

    let root = SystemRoot { key: root_key, seed, name: sys_name.clone(), epoch: Epoch { system_age_gyr: age }, placement: None, site: Some(ph.site.clone()) };
    let system = SolarSystem { root, nodes };
    match focus(&system) {
        Some(focus) => Ok(Generated { system, focus }),
        None => Err(Barren { reason: format!("{sys_name} is {} with nothing solid in orbit", describe_star(&primary)), system: Box::new(system) }),
    }
}

/// "a red dwarf (M3V, 0.31 M☉)", "a black hole (12.3 M☉)".
#[must_use]
pub fn describe_star(s: &Star) -> String {
    let what = match s.stage {
        StellarStage::MainSequence if s.temperature_k < 3_700.0 => "red dwarf",
        StellarStage::MainSequence if s.temperature_k > 10_000.0 => "blue star",
        stage => stage.label(),
    };
    match s.stage {
        StellarStage::NeutronStar | StellarStage::BlackHole => format!("a {what} ({:.1} M☉)", s.mass_solar),
        _ => format!("a {what} ({}, {:.2} M☉)", s.spectral_class, s.mass_solar),
    }
}

fn star_node(key: Key, name: String, parent: Option<Key>, star: Star, elements: Option<Elements>) -> OrbitalNode {
    OrbitalNode { key, name, parent, node: NodeKind::Star(star), elements }
}

fn orbit(d: &mut Draw, a_au: f64, e: f64) -> Elements {
    Elements { a_au, e, i_deg: d.range(0.0, 10.0), lan_deg: 0.0, argp_deg: d.range(0.0, 360.0), m0_deg: d.range(0.0, 360.0) }
}

/// The disk the planets formed in: around the primary inside a wide companion's reach, or
/// around a close pair outside its stars' orbits, and inside a triple's outer star.
fn birth_disk(ph: &Phenotype, close: bool) -> Disk {
    let m1 = ph.primary_mass;
    let reach = |c: &Companion, inner: f64| s_type_limit(c.initial_mass / (inner + c.initial_mass), c.e) * c.a_au;
    let outer_limit = ph.companions.get(1).map_or(f64::INFINITY, |c| reach(c, m1 + ph.companions[0].initial_mass));
    match ph.companions.first() {
        Some(c) if close => {
            let mu = c.initial_mass / (m1 + c.initial_mass);
            Disk {
                mass: m1 + c.initial_mass,
                luminosity: main_sequence_luminosity(m1) + main_sequence_luminosity(c.initial_mass),
                metallicity: ph.metallicity_dex,
                inner_au: p_type_limit(mu, c.e) * c.a_au,
                outer_au: outer_limit,
                survival: 1.0,
            }
        }
        Some(c) => {
            Disk { mass: m1, luminosity: main_sequence_luminosity(m1), metallicity: ph.metallicity_dex, inner_au: 0.0, outer_au: reach(c, m1).min(outer_limit), survival: binary_survival(c.a_au) }
        }
        None => Disk { mass: m1, luminosity: main_sequence_luminosity(m1), metallicity: ph.metallicity_dex, inner_au: 0.0, outer_au: f64::INFINITY, survival: 1.0 },
    }
}

/// Orbits closer than this to a host of `mass` M☉ (AU) are tidally locked after `age_gyr`.
#[must_use]
pub fn tidal_lock_au(mass: f64, age_gyr: f64) -> f64 {
    0.35 * mass.cbrt() * (age_gyr.max(0.01) / 4.5).powf(1.0 / 6.0)
}

fn planet(d: &mut Draw, p: &planets::Formed, host_mass: f64, age: f64) -> BodyExterior {
    let period_h = orbital_period_days(p.a_au, host_mass * SOLAR_MASS_KG) * 24.0;
    let locked = p.a_au < tidal_lock_au(host_mass, age);
    let giant = matches!(p.class, BodyClass::GasGiant | BodyClass::IceGiant);
    let spin = if giant { d.range(9.0, 18.0) } else { d.range(14.0, 40.0) };
    let tilt = d.range(0.0, 35.0);
    let rotation_h = match locked {
        true if p.e > 0.1 => period_h * 2.0 / 3.0,
        true => period_h,
        false => spin,
    };
    BodyExterior {
        class: p.class,
        mass_kg: p.mass_earth * EARTH_MASS_KG,
        radius_km: planet_radius_earth(p.class, p.mass_earth) * EARTH_RADIUS_KM,
        rotation_h,
        obliquity_deg: if locked { tilt / 35.0 } else { tilt },
    }
}

fn belts(d: &mut Draw, keys: &mut KeyAllocator, host: Key, sys_name: &str, arch: &Architecture) -> Vec<OrbitalNode> {
    let make = |keys: &mut KeyAllocator, name: String, (inner, outer): Belt, pop: Population| OrbitalNode {
        key: keys.fresh(),
        name,
        parent: Some(host),
        node: NodeKind::Belt(Population { inner_au: inner, outer_au: outer, ..pop }),
        elements: None,
    };
    let mut out = Vec::new();
    if let Some(b) = arch.asteroid_belt {
        let pop = Population {
            inner_au: 0.0,
            outer_au: 0.0,
            count_estimate: 10f64.powf(d.range(5.5, 6.5)),
            size_law: PowerLaw { min_km: 1.0, max_km: d.range(400.0, 950.0), exponent: 3.5 },
            composition: Composition::from_weights(&[(Species::Silicate, 0.75), (Species::Iron, 0.15), (Species::Ice, 0.10)]),
        };
        out.push(make(keys, format!("{sys_name} Belt"), b, pop));
    }
    if let Some(b) = arch.outer_belt {
        let pop = Population {
            inner_au: 0.0,
            outer_au: 0.0,
            count_estimate: 10f64.powf(d.range(8.0, 9.0)),
            size_law: PowerLaw { min_km: 1.0, max_km: d.range(1000.0, 2500.0), exponent: 3.0 },
            composition: Composition::from_weights(&[(Species::Ice, 0.6), (Species::Silicate, 0.4)]),
        };
        out.push(make(keys, format!("{sys_name} Outer Belt"), b, pop));
    }
    out
}

/// How Earth-like a body is: liquid water at the surface above all, then plates and air, then
/// starlight and mass near Earth's. `None` for bodies without a solid surface.
#[must_use]
pub fn earthlike_score(ctx: &BodyContext) -> Option<f64> {
    let i = &ctx.interior;
    if !i.has_solid_surface() {
        return None;
    }
    let s = (ctx.irradiance.flux_w_m2 / 1361.0).max(1e-9);
    let m = ctx.exterior.mass_kg / EARTH_MASS_KG;
    let mut score = -s.ln().abs() - 0.5 * m.ln().abs();
    if i.has_surface_liquid_water() {
        score += 5.0;
    }
    if i.regime == TectonicRegime::Plates {
        score += 1.0;
    }
    if i.regime == TectonicRegime::IceShell {
        score += 0.5;
    }
    if i.has_atmosphere() {
        score += 1.0;
    }
    Some(score)
}

/// The most Earth-like round body; else a belt's largest member; else `None`.
fn focus(system: &SolarSystem) -> Option<Key> {
    let best = system
        .nodes
        .iter()
        .filter(|n| matches!(n.node, NodeKind::Body(b) if b.class == BodyClass::Terrestrial))
        .filter_map(|n| body_context(system, n.key).and_then(|c| earthlike_score(&c)).map(|s| (n.key, s)))
        .max_by(|a, b| a.1.total_cmp(&b.1));
    if let Some((key, _)) = best {
        return Some(key);
    }
    system.nodes.iter().filter(|n| matches!(n.node, NodeKind::Belt(_))).find_map(|b| belt_member(system, b.key, 0)).map(|c| c.key)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::context::any_body;

    #[test]
    fn deterministic_and_complete() {
        let g = Galaxy::milky_way();
        let mut barren = 0;
        for s in 0..300u64 {
            let a = galactic_system(Seed(s), &g);
            assert_eq!(a, galactic_system(Seed(s), &g));
            match a {
                Ok(gen) => {
                    assert!(any_body(&gen.system, gen.focus).is_some(), "seed {s}: focus has a context");
                    for n in &gen.system.nodes {
                        if let NodeKind::Body(_) = n.node {
                            let c = body_context(&gen.system, n.key).unwrap_or_else(|| panic!("seed {s}: {}", n.name));
                            assert!(c.irradiance.flux_w_m2.is_finite() && c.irradiance.year_days.is_finite());
                        }
                    }
                }
                Err(_) => barren += 1,
            }
        }
        assert!(barren < 150, "most systems have something solid: {barren} of 300 barren");
    }

    #[test]
    fn red_dwarf_worlds_are_locked() {
        let g = Galaxy::milky_way();
        let mut checked = 0;
        for s in 0..400u64 {
            let Ok(gen) = galactic_system(Seed(s), &g) else { continue };
            let ctx = any_body(&gen.system, gen.focus).unwrap();
            let host = gen.system.nodes.iter().find_map(|n| if let NodeKind::Star(st) = &n.node { Some(st.clone()) } else { None }).unwrap();
            if host.stage.is_main_sequence() && host.mass_solar < 0.4 && ctx.interior.has_surface_liquid_water() && gen.system.nodes.iter().filter(|n| matches!(n.node, NodeKind::Star(_))).count() == 1
            {
                let period_h = orbital_period_days(ctx.orbit.a_au, host.mass_solar * SOLAR_MASS_KG) * 24.0;
                assert!((ctx.exterior.rotation_h / period_h - 1.0).abs() < 0.34, "seed {s}: day {} h, year {period_h} h", ctx.exterior.rotation_h);
                checked += 1;
            }
        }
        assert!(checked > 0, "some red dwarf has a wet world");
    }
}
