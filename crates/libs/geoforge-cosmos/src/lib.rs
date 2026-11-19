//! Thin cosmology: one solar system as the world's root.
//!
//! The system is a pure function of its own key and seed (never of ancestors), so it can
//! later be attached to a galaxy sector without changing anything inside it. Its job is to
//! produce a [`BodyContext`](geoforge_core::BodyContext) for any of its bodies; the surface
//! pipeline never sees this crate.
//!
//! | Module | Concern |
//! |---|---|
//! | [`system`] | the orbital tree: star, planets, belts, moons |
//! | [`worlds`] | what each body is made of: shells, air, water, regime; the CO₂ thermostat |
//! | [`context`] | the handoff for a body, belt members on demand, lookup by name |

mod draw;

pub mod context;
pub mod system;
pub mod worlds;

pub use context::{belt_member, body_context, find_body};
pub use draw::name;
pub use system::{generate_system, Generated};
pub use worlds::{describe, interior_for, regulate_co2, Setting};

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
            assert!((10.0..26.0).contains(&t), "seed {s}: balance temperature {t}");
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
