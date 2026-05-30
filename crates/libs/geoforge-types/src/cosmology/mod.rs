#![allow(missing_docs)]

//! Large-scale cosmology stubs (Stage 0) — minimal, implemented types only.

use serde::{Deserialize, Serialize};

/// Stellar spectral classification (Morgan–Keenan simplified).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SpectralClass {
    /// Hottest, bluest.
    O,
    B,
    A,
    F,
    /// Sun-like.
    G,
    K,
    /// Coolest, reddest.
    M,
}

/// Star parameters for insolation and habitability (Stage 0).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Star {
    /// Luminosity relative to the Sun.
    pub luminosity_solar: f64,
    /// Mass (solar masses).
    pub mass_solar: f64,
    /// Effective temperature (K).
    pub temperature_k: f64,
    /// Age (Gyr).
    pub age_gyr: f64,
    /// Metallicity [Fe/H] dex.
    pub metallicity_dex: f64,
    /// Spectral class.
    pub spectral_class: SpectralClass,
}

impl Star {
    /// Sun-like star.
    #[must_use]
    pub fn sun() -> Self {
        Self {
            luminosity_solar: 1.0,
            mass_solar: 1.0,
            temperature_k: 5778.0,
            age_gyr: 4.6,
            metallicity_dex: 0.0,
            spectral_class: SpectralClass::G,
        }
    }
}

/// Galaxy morphology (placeholder for future Stage 0 expansion).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum GalaxyType {
    /// Elliptical.
    Elliptical,
    /// Spiral.
    Spiral,
    /// Irregular.
    Irregular,
}
