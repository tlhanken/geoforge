//! What a body must have for a generator to apply. Each generator's `applies` combines these,
//! so the rules that pick between implementations of a layer live in one place.

use geoforge_core::body::TectonicRegime;
use geoforge_core::product::Precondition;
use geoforge_core::BodyContext;

type Check = Result<(), Precondition>;

fn require(ok: bool, why: &str) -> Check {
    if ok { Ok(()) } else { Err(Precondition(why.into())) }
}

/// Rock crust broken into moving plates.
pub fn plate_tectonics(body: &BodyContext) -> Check {
    require(body.interior.has_lithosphere(), "needs a rocky crust")?;
    require(body.interior.regime == TectonicRegime::Plates, "needs plate tectonics")
}

/// Rock crust under one rigid lid: stagnant, episodic, or geologically dead.
pub fn rigid_lid(body: &BodyContext) -> Check {
    require(body.interior.has_lithosphere(), "needs a rocky crust")?;
    require(
        matches!(body.interior.regime, TectonicRegime::StagnantLid | TectonicRegime::Episodic | TectonicRegime::Inactive),
        "needs a rigid lid (stagnant, episodic or inactive)",
    )
}

/// Ground under air.
pub fn surface_air(body: &BodyContext) -> Check {
    require(body.interior.has_solid_surface(), "needs a solid surface")?;
    require(body.interior.has_atmosphere(), "needs an atmosphere")
}

/// Liquid water on the surface, so it rains, flows and pools.
pub fn surface_water(body: &BodyContext) -> Check {
    surface_air(body)?;
    require(body.interior.has_surface_liquid_water(), "needs liquid water on the surface")
}

#[cfg(test)]
mod tests {
    use super::*;
    use geoforge_core::{Archetype, Key};

    #[test]
    fn regimes_pick_one_tectonic_model() {
        for a in Archetype::ALL {
            let b = a.context(Key(1), a.slug());
            let both = plate_tectonics(&b).is_ok() && rigid_lid(&b).is_ok();
            assert!(!both, "{a}: plate and lid models must exclude each other");
        }
        let ok = |a: Archetype, f: fn(&BodyContext) -> Check| f(&a.context(Key(1), "x")).is_ok();
        assert!(ok(Archetype::Earth, plate_tectonics));
        assert!(ok(Archetype::Mars, rigid_lid) && ok(Archetype::Venus, rigid_lid) && ok(Archetype::Moon, rigid_lid));
        assert!(!ok(Archetype::Europa, rigid_lid) && !ok(Archetype::GasGiant, rigid_lid));
        assert!(ok(Archetype::Earth, surface_water) && !ok(Archetype::Mars, surface_water));
        assert!(ok(Archetype::Titan, surface_air) && !ok(Archetype::Moon, surface_air));
    }
}
