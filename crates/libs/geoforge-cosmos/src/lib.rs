//! Thin cosmology: one solar system as the world's root.
//!
//! The system is a pure function of its own key and seed (never of ancestors), so it can
//! later be attached to a galaxy sector without changing anything inside it. Its job is to
//! produce a [`BodyContext`](geoforge_core::BodyContext) for any of its bodies; the surface
//! pipeline never sees this crate.
//!
//! Two recipes make a system ([`Recipe`]): the habitable one, a Sun-like star with an
//! Earth-like world ([`generate_system`]), and a draw from a galaxy's stellar population, any
//! star alive or dead with whatever planets it kept ([`galactic_system`]).
//!
//! | Module | Concern |
//! |---|---|
//! | [`galaxy`] | galaxy models: components, shapes, arms, dust, morphologies |
//! | [`sectors`] | a galaxy's stars generated where they are looked at, by nested cubes |
//! | [`starfield`] | the stars with a world's edits over them |
//! | [`population`] | a system's phenotype drawn from a galaxy (or taken from a star): site, age, metals, stars at birth |
//! | [`stars`] | stellar evolution: main sequence, giants, white dwarfs, neutron stars, black holes |
//! | `planets` | planets for any star at birth, and what survives the star's death |
//! | [`galactic`] | the orbital tree of a drawn system, and its focus |
//! | [`census`] | statistics of many drawn systems, with example seeds |
//! | [`sky`] | the naked-eye stars around a system, dimmed by dust, and its galaxy's glow |
//! | [`sol`] | the real Solar System from published elements, at a moment of its calendar |
//! | [`system`] | the habitable recipe's orbital tree: star, planets, belts, moons |
//! | [`worlds`] | what each body is made of: shells, air, water, regime; the CO₂ thermostat |
//! | [`context`] | the handoff for a body, belt members on demand, lookup by name |

mod draw;
mod dynamics;

pub mod census;
pub mod context;
pub mod galactic;
pub mod galaxy;
mod planets;
pub mod population;
pub mod sectors;
pub mod sky;
pub mod sol;
pub mod starfield;
pub mod stars;
pub mod system;
pub mod worlds;

pub use context::{any_body, belt_member, body_context, find_body};
pub use galactic::{galactic_system, Barren};
pub use galaxy::{Galaxy, GalaxyPreset};
pub use starfield::Starfield;
pub use system::{generate_system, Generated};
pub use worlds::{describe, interior_for, regulate_co2, Setting};

/// How a system is made.
#[derive(Debug, Clone, PartialEq)]
pub enum Recipe {
    /// A Sun-like star with an Earth-like world in its habitable zone ([`generate_system`]).
    Habitable,
    /// A system drawn from a galaxy's stellar population ([`galactic_system`]).
    Galaxy(Galaxy),
}

impl Default for Recipe {
    fn default() -> Self {
        Self::Galaxy(Galaxy::milky_way())
    }
}

impl Recipe {
    /// The system of `seed`.
    ///
    /// # Errors
    /// When a drawn system has no solid body to focus on.
    pub fn generate(&self, seed: geoforge_core::seed::Seed) -> Result<Generated, Barren> {
        match self {
            Self::Habitable => Ok(generate_system(seed)),
            Self::Galaxy(g) => galactic_system(seed, g),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use geoforge_core::body::{BodyClass, NodeKind, TectonicRegime};
    use geoforge_core::radiation::balance_temperature_c;
    use geoforge_core::seed::Seed;

    #[test]
    fn systems_are_plausible_and_deterministic() {
        for s in 0..40u64 {
            let g = generate_system(Seed(s));
            assert_eq!(g.system, generate_system(Seed(s)).system);
            let ctx = body_context(&g.system, g.focus).unwrap();
            // Focus world: Earth-ish starlight, an ocean, plate tectonics, a mild climate.
            assert!((0.8..1.25).contains(&(ctx.irradiance.flux_w_m2 / 1361.0)), "flux {}", ctx.irradiance.flux_w_m2);
            assert_eq!(ctx.interior.regime, TectonicRegime::Plates);
            assert!(ctx.interior.has_lithosphere() && ctx.interior.surface_ocean_gel_km() > 1.0);
            let t = balance_temperature_c(&ctx);
            // The thermostat holds weathering, not the mean: fast spinners run colder overall.
            assert!((5.0..26.0).contains(&t), "seed {s}: balance temperature {t}");
            assert!((5.0..15.0).contains(&ctx.gravity_m_s2()), "gravity {}", ctx.gravity_m_s2());
            assert!((50.0..1000.0).contains(&ctx.irradiance.year_days), "year {}", ctx.irradiance.year_days);
            // Every body of the tree has a context.
            for n in &g.system.nodes {
                if let NodeKind::Body(_) = n.node {
                    assert!(body_context(&g.system, n.key).is_some(), "{}", n.name);
                }
            }
        }
    }

    #[test]
    fn giants_have_moons_and_moons_share_the_year() {
        let mut found = false;
        for s in 0..20u64 {
            let g = generate_system(Seed(s));
            for giant in g.system.nodes.iter().filter(|n| matches!(n.node, NodeKind::Body(b) if b.class == BodyClass::GasGiant)) {
                let planet = body_context(&g.system, giant.key).unwrap();
                for m in g.system.children(giant.key) {
                    found = true;
                    let moon = body_context(&g.system, m.key).unwrap();
                    assert_eq!(moon.irradiance.flux_w_m2, planet.irradiance.flux_w_m2);
                    assert_eq!(moon.irradiance.year_days, planet.irradiance.year_days);
                    assert_eq!(moon.irradiance.obliquity_deg, planet.exterior.obliquity_deg);
                    assert_eq!(moon.tides[0].name, giant.name);
                    // Tidally locked.
                    let period_h = geoforge_core::derived::orbital_period_days(moon.orbit.a_au, planet.exterior.mass_kg) * 24.0;
                    assert!((moon.exterior.rotation_h - period_h).abs() < 1e-6);
                }
            }
        }
        assert!(found, "some gas giant has moons");
    }

    /// Planets of one host at least 8 mutual Hill radii apart (eccentric survivors of
    /// scattering aside), moons between the Roche limit and 0.4 of the Hill radius, and an
    /// asteroid belt inside its giant's 2:1 resonance, in both recipes.
    #[test]
    fn orbits_are_spaced_and_moons_held() {
        use geoforge_core::body::{EARTH_MASS_KG, SOLAR_MASS_KG};
        let galaxy = galaxy::Galaxy::milky_way();
        let systems = (0..30u64).map(|s| generate_system(Seed(s))).chain((0..200u64).filter_map(|s| galactic_system(Seed(s), &galaxy).ok()));
        let (mut planets_seen, mut moons_seen) = (0, 0);
        for g in systems {
            let sys = &g.system;
            let star_mass = |key| match &sys.node(key).map(|n| &n.node) {
                Some(NodeKind::Star(s)) => s.mass_solar,
                _ => sys.children(key).iter().filter_map(|c| if let NodeKind::Star(s) = &c.node { Some(s.mass_solar) } else { None }).sum(),
            };
            for host in sys.nodes.iter().filter(|n| !matches!(n.node, NodeKind::Body(_) | NodeKind::Belt(_) | NodeKind::Ring(_))) {
                let mut kids: Vec<(f64, f64)> = sys
                    .children(host.key)
                    .iter()
                    .filter_map(|c| match (&c.node, c.elements) {
                        (NodeKind::Body(b), Some(e)) => Some((e.a_au, b.mass_kg / EARTH_MASS_KG)),
                        _ => None,
                    })
                    .collect();
                kids.sort_by(|a, b| a.0.total_cmp(&b.0));
                let m = star_mass(host.key);
                for w in kids.windows(2) {
                    let h = ((w[0].1 + w[1].1) * EARTH_MASS_KG / (3.0 * m * SOLAR_MASS_KG)).cbrt() * 0.5 * (w[0].0 + w[1].0);
                    assert!((w[1].0 - w[0].0) / h >= 7.99, "{}: {:?} only {:.1} mutual Hill radii apart", sys.root.name, w, (w[1].0 - w[0].0) / h);
                    planets_seen += 1;
                }
            }
            for n in &sys.nodes {
                let (Some(parent), Some(e), NodeKind::Body(moon)) = (n.parent, n.elements, &n.node) else { continue };
                let Some(NodeKind::Body(planet)) = sys.node(parent).map(|p| &p.node) else { continue };
                let p_el = sys.node(parent).and_then(|p| p.elements).unwrap();
                let moon_density = moon.mass_kg / (4.0 / 3.0 * std::f64::consts::PI * (moon.radius_km * 1e5).powi(3)) * 1000.0;
                let a_km = e.a_au * geoforge_core::body::AU_KM;
                let hill = system::hill_km(p_el.a_au, planet.mass_kg, star_mass(sys.node(parent).unwrap().parent.unwrap()));
                assert!(a_km > system::roche_km(planet, moon_density), "{} inside the Roche limit", n.name);
                assert!(a_km < system::MOON_HILL_SHARE * hill, "{} outside the Hill sphere", n.name);
                moons_seen += 1;
            }
        }
        assert!(planets_seen > 300 && moons_seen > 100, "{planets_seen} neighbours, {moons_seen} moons");
    }

    #[test]
    fn belt_members_on_demand() {
        let g = generate_system(Seed(3));
        let belt = g.system.nodes.iter().find(|n| n.name.ends_with("Outer Belt")).unwrap();
        let first = belt_member(&g.system, belt.key, 0).unwrap();
        let tenth = belt_member(&g.system, belt.key, 9).unwrap();
        assert!(first.exterior.radius_km >= tenth.exterior.radius_km);
        assert_eq!(belt_member(&g.system, belt.key, 9).unwrap(), tenth, "deterministic");
        assert_ne!(first.key, tenth.key);
        // By name, 1-based.
        let by_name = find_body(&g.system, &format!("{} #10", belt.name)).unwrap();
        assert_eq!(by_name, tenth);
        assert!(find_body(&g.system, "nowhere #1").is_none());
        // A small member is an inactive body with a crust: the lid generators crater it.
        let small = belt_member(&g.system, belt.key, 5000).unwrap();
        assert_eq!(small.interior.regime, TectonicRegime::Inactive);
        assert!(small.interior.has_lithosphere());
    }
}
