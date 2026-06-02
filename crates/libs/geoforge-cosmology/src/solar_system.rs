//! Full solar-system generation at zoom LOD.

use geoforge_types::cosmology::{
    BeltSlot, Multiplicity, PlanetMassClass, PlanetSlot, SolarSystem, SpectralClass, Star,
    StellarSystemRef, SystemPhenotype,
};
use geoforge_types::coordinates::OrbitalElements;
use geoforge_types::cosmology::BarycentricOffset;
use geoforge_types::Seed;

use crate::marker::marker_position;
use crate::rng::{draw, range, unit};
use crate::seeds::{galaxy_seed, system_seed};

/// Generate a solar system from a stable reference (regenerable).
#[must_use]
pub fn generate_solar_system(
    root: Seed,
    galaxy_index: u32,
    profile: &geoforge_types::cosmology::GalaxyProfile,
    reference: StellarSystemRef,
    phenotype: SystemPhenotype,
) -> SolarSystem {
    let gseed = galaxy_seed(root, galaxy_index);
    let sys = system_seed(
        crate::seeds::region_seed(gseed, reference.region),
        reference.index_in_region,
    );

    let barycenter_ly = marker_position(
        profile,
        reference.region,
        reference.index_in_region,
        sys,
    );

    let stars = generate_stars(sys, &phenotype);
    let planet_slots = generate_planet_slots(sys, stars.len());
    let belts = generate_belts(sys);

    SolarSystem {
        galaxy_index,
        region: reference.region,
        index_in_region: reference.index_in_region,
        barycenter_ly,
        phenotype,
        stars,
        planet_slots,
        belts,
    }
}

fn generate_stars(sys: Seed, phenotype: &SystemPhenotype) -> Vec<Star> {
    let mut stars = Vec::new();
    let n = match phenotype.multiplicity {
        Multiplicity::Single => 1,
        Multiplicity::Binary => 2,
        Multiplicity::Trinary => 3,
    };

    for i in 0..n {
        let class = if i == 0 {
            phenotype.primary_class
        } else {
            phenotype.secondary_class.unwrap_or(SpectralClass::M)
        };
        let offset = if i == 0 {
            BarycentricOffset::origin()
        } else {
            let sep = range(draw(sys, b"sep"), 5.0, 40.0);
            BarycentricOffset {
                x_au: sep,
                y_au: 0.0,
                z_au: 0.0,
            }
        };
        stars.push(star_from_class(class, offset));
    }
    stars
}

fn star_from_class(class: SpectralClass, position_au: BarycentricOffset) -> Star {
    let (lum, mass, temp) = match class {
        SpectralClass::O => (30_000.0, 16.0, 40_000.0),
        SpectralClass::B => (1_000.0, 7.0, 20_000.0),
        SpectralClass::A => (80.0, 2.5, 8_500.0),
        SpectralClass::F => (6.0, 1.5, 6_500.0),
        SpectralClass::G => (1.0, 1.0, 5_778.0),
        SpectralClass::K => (0.4, 0.8, 4_500.0),
        SpectralClass::M => (0.05, 0.3, 3_000.0),
    };
    Star {
        luminosity_solar: lum,
        mass_solar: mass,
        temperature_k: temp,
        age_gyr: 4.6,
        metallicity_dex: 0.0,
        spectral_class: class,
        position_au,
    }
}

fn generate_planet_slots(sys: Seed, star_count: usize) -> Vec<PlanetSlot> {
    let n = (unit(draw(sys, b"planet_count")) * 10.0) as u8;
    let mut slots = Vec::new();
    let mut au = 0.3 + 0.2 * star_count as f64;
    for i in 0..n {
        au *= 1.4 + unit(draw(sys, b"orbit").child_bytes(&[i]));
        let mass_class = match (unit(draw(sys, b"mass").child_bytes(&[i])) * 5.0) as u8 {
            0 => PlanetMassClass::Rocky,
            1 => PlanetMassClass::SuperEarth,
            2 => PlanetMassClass::GasGiant,
            3 => PlanetMassClass::IceGiant,
            _ => PlanetMassClass::Dwarf,
        };
        slots.push(PlanetSlot {
            index: i,
            orbit: OrbitalElements {
                period_days: 365.25 * au.powf(1.5),
                eccentricity: range(draw(sys, b"ecc").child_bytes(&[i]), 0.0, 0.2),
                semi_major_axis_au: au,
                inclination_deg: range(draw(sys, b"inc").child_bytes(&[i]), 0.0, 5.0),
                phase: unit(draw(sys, b"phase").child_bytes(&[i])),
            },
            mass_class,
        });
    }
    slots
}

fn generate_belts(sys: Seed) -> Vec<BeltSlot> {
    if unit(draw(sys, b"belt")) < 0.3 {
        return Vec::new();
    }
    vec![BeltSlot {
        index: 0,
        inner_au: 2.0,
        outer_au: 3.5,
    }]
}
