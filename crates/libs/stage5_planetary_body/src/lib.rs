//! Pipeline stage 5 — planetary body parameters from solar-system zoom.

#![deny(clippy::todo, clippy::unimplemented)]

use geoforge_types::cosmology::{PlanetMassClass, PlanetSlot, SolarSystem};
use geoforge_types::planetary::{
    PlanetBodyRecord, PlanetaryBodies, PlanetaryParams, DEFAULT_MAP_HEIGHT, DEFAULT_MAP_WIDTH,
};
use geoforge_types::seed::{draw, range, unit};

/// Build planetary parameter records for every planet slot (no surface meshes yet).
#[must_use]
pub fn generate(solar: &SolarSystem) -> PlanetaryBodies {
    let primary = solar
        .stars
        .first()
        .expect("solar system has at least one star");

    let bodies: Vec<PlanetBodyRecord> = solar
        .planet_slots
        .iter()
        .map(|slot| record_for_slot(solar.system_seed, slot, primary.luminosity_solar))
        .collect();

    let focus_slot_index = pick_focus_slot(&bodies);

    PlanetaryBodies {
        focus_slot_index,
        bodies,
    }
}

fn record_for_slot(
    system: geoforge_types::Seed,
    slot: &PlanetSlot,
    stellar_luminosity: f64,
) -> PlanetBodyRecord {
    let slot_seed = system.child(b"planet").child_bytes(&[slot.index]);
    let mut params = base_params(slot_seed, slot.mass_class);
    apply_orbit_and_star(&mut params, slot, stellar_luminosity);
    apply_variation(&mut params, slot_seed, slot.mass_class);

    PlanetBodyRecord {
        slot_index: slot.index,
        mass_class: slot.mass_class,
        orbit: slot.orbit,
        params,
        map_width: DEFAULT_MAP_WIDTH,
        map_height: DEFAULT_MAP_HEIGHT,
    }
}

fn base_params(seed: geoforge_types::Seed, class: PlanetMassClass) -> PlanetaryParams {
    match class {
        PlanetMassClass::Rocky => PlanetaryParams::from_radius_and_density(
            range(draw(seed, b"radius"), 2_400.0, 7_200.0),
            range(draw(seed, b"density"), 3_500.0, 5_800.0),
        ),
        PlanetMassClass::SuperEarth => PlanetaryParams::from_radius_and_density(
            range(draw(seed, b"radius"), 10_000.0, 20_000.0),
            range(draw(seed, b"density"), 4_500.0, 6_500.0),
        ),
        PlanetMassClass::GasGiant => PlanetaryParams::from_radius_and_density(
            range(draw(seed, b"radius"), 40_000.0, 80_000.0),
            range(draw(seed, b"density"), 600.0, 1_400.0),
        ),
        PlanetMassClass::IceGiant => PlanetaryParams::from_radius_and_density(
            range(draw(seed, b"radius"), 20_000.0, 35_000.0),
            range(draw(seed, b"density"), 1_200.0, 1_800.0),
        ),
        PlanetMassClass::Dwarf => PlanetaryParams::from_radius_and_density(
            range(draw(seed, b"radius"), 500.0, 2_400.0),
            range(draw(seed, b"density"), 2_000.0, 3_500.0),
        ),
    }
}

fn apply_orbit_and_star(params: &mut PlanetaryParams, slot: &PlanetSlot, stellar_luminosity: f64) {
    params.semi_major_axis_au = slot.orbit.semi_major_axis_au;
    params.orbital_period_days = slot.orbit.period_days;
    params.orbital_eccentricity = slot.orbit.eccentricity;
    params.orbital_inclination_degrees = slot.orbit.inclination_deg;
    params.stellar_luminosity = stellar_luminosity;
    let ecc = slot.orbit.eccentricity;
    params.perihelion_au = params.semi_major_axis_au * (1.0 - ecc);
    params.aphelion_au = params.semi_major_axis_au * (1.0 + ecc);
}

fn apply_variation(params: &mut PlanetaryParams, seed: geoforge_types::Seed, class: PlanetMassClass) {
    params.axial_tilt_degrees = range(draw(seed, b"tilt"), 0.0, 45.0);
    params.rotation_period_hours =
        range(draw(seed, b"rotation"), 10.0, 240.0) * (0.5 + unit(draw(seed, b"rot2")));
    let pressure_scale = range(draw(seed, b"atmo"), 0.0, 1.5);
    let atmo_mult = match class {
        PlanetMassClass::GasGiant | PlanetMassClass::IceGiant => 5.0,
        PlanetMassClass::Dwarf => 0.01,
        PlanetMassClass::Rocky | PlanetMassClass::SuperEarth => 1.0,
    };
    params.atmospheric_pressure_kpa =
        PlanetaryParams::earth().atmospheric_pressure_kpa * pressure_scale * atmo_mult;
    params.greenhouse_factor =
        range(draw(seed, b"gh"), 0.1, 1.5) * PlanetaryParams::earth().greenhouse_factor;
}

fn pick_focus_slot(bodies: &[PlanetBodyRecord]) -> u8 {
    let habitable: Vec<_> = bodies
        .iter()
        .filter(|b| {
            matches!(
                b.mass_class,
                PlanetMassClass::Rocky | PlanetMassClass::SuperEarth
            ) && (0.4..=2.5).contains(&b.orbit.semi_major_axis_au)
        })
        .collect();

    if let Some(best) = habitable
        .iter()
        .min_by(|a, b| {
            a.orbit
                .semi_major_axis_au
                .partial_cmp(&b.orbit.semi_major_axis_au)
                .unwrap_or(std::cmp::Ordering::Equal)
        }) {
        return best.slot_index;
    }

    bodies
        .first()
        .map(|b| b.slot_index)
        .unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use geoforge_types::cosmology::{PlanetMassClass, PlanetSlot};
    use geoforge_types::coordinates::OrbitalElements;
    use geoforge_types::cosmology::{Multiplicity, SolarSystem, SystemPhenotype};
    use geoforge_types::cosmology::SpectralClass;
    use geoforge_types::Seed;

    fn minimal_solar() -> SolarSystem {
        SolarSystem {
            system_seed: Seed::new(1),
            galaxy_index: 0,
            region: geoforge_types::cosmology::RegionId::Cylindrical {
                ring: 0,
                wedge: 0,
                layer: 0,
            },
            index_in_region: 0,
            barycenter_ly: geoforge_types::coordinates::GalacticPoint::origin(),
            phenotype: SystemPhenotype {
                multiplicity: Multiplicity::Single,
                primary_class: SpectralClass::G,
                secondary_class: None,
            },
            stars: vec![geoforge_types::cosmology::Star {
                luminosity_solar: 1.0,
                mass_solar: 1.0,
                temperature_k: 5778.0,
                age_gyr: 4.6,
                metallicity_dex: 0.0,
                spectral_class: SpectralClass::G,
                position_au: geoforge_types::cosmology::BarycentricOffset::origin(),
            }],
            planet_slots: vec![PlanetSlot {
                index: 0,
                mass_class: PlanetMassClass::Rocky,
                orbit: OrbitalElements {
                    period_days: 365.0,
                    eccentricity: 0.0,
                    semi_major_axis_au: 1.0,
                    inclination_deg: 0.0,
                    phase: 0.0,
                },
            }],
            belts: vec![],
        }
    }

    #[test]
    fn one_record_per_slot() {
        let solar = minimal_solar();
        let n = solar.planet_slots.len();
        let bodies = generate(&solar);
        assert_eq!(bodies.bodies.len(), n);
        assert!(bodies.focus().is_some());
    }
}
