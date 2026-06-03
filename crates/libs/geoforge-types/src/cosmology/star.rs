//! Stellar types — used in map phenotypes and solar-system zoom.

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

/// Star at solar-system zoom (barycentric AU offsets optional).
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
    /// Offset from system barycenter (AU) for binaries.
    pub position_au: BarycentricOffset,
}

/// 3D offset from the system barycenter in AU.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct BarycentricOffset {
    /// X (AU).
    pub x_au: f64,
    /// Y (AU).
    pub y_au: f64,
    /// Z (AU).
    pub z_au: f64,
}

impl BarycentricOffset {
    /// Origin at barycenter.
    #[must_use]
    pub const fn origin() -> Self {
        Self {
            x_au: 0.0,
            y_au: 0.0,
            z_au: 0.0,
        }
    }
}

impl Star {
    /// Sun-like star at barycenter.
    #[must_use]
    pub fn sun() -> Self {
        Self {
            luminosity_solar: 1.0,
            mass_solar: 1.0,
            temperature_k: 5778.0,
            age_gyr: 4.6,
            metallicity_dex: 0.0,
            spectral_class: SpectralClass::G,
            position_au: BarycentricOffset::origin(),
        }
    }
}
