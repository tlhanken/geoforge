//! Map-level system phenotype from system seed.

use geoforge_types::cosmology::{Multiplicity, SpectralClass, SystemPhenotype};
use geoforge_types::Seed;

use crate::rng::{draw, unit};

/// Derive phenotype (must match solar-system star generation).
#[must_use]
pub fn derive_phenotype(system: Seed, force_single_g: bool) -> SystemPhenotype {
    if force_single_g {
        return SystemPhenotype {
            multiplicity: Multiplicity::Single,
            primary_class: SpectralClass::G,
            secondary_class: None,
        };
    }

    let roll = unit(draw(system, b"multiplicity"));
    let multiplicity = if roll < 0.50 {
        Multiplicity::Single
    } else if roll < 0.90 {
        Multiplicity::Binary
    } else {
        Multiplicity::Trinary
    };

    let primary_class = sample_spectral(draw(system, b"primary"));
    let secondary_class = match multiplicity {
        Multiplicity::Single => None,
        Multiplicity::Binary | Multiplicity::Trinary => {
            Some(sample_spectral(draw(system, b"secondary")))
        }
    };

    SystemPhenotype {
        multiplicity,
        primary_class,
        secondary_class,
    }
}

#[must_use]
fn sample_spectral(seed: Seed) -> SpectralClass {
    let u = unit(seed);
    if u < 0.000_1 {
        SpectralClass::O
    } else if u < 0.001 {
        SpectralClass::B
    } else if u < 0.01 {
        SpectralClass::A
    } else if u < 0.05 {
        SpectralClass::F
    } else if u < 0.20 {
        SpectralClass::G
    } else if u < 0.50 {
        SpectralClass::K
    } else {
        SpectralClass::M
    }
}
