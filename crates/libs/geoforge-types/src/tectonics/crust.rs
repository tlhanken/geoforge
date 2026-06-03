//! Crust composition — per-pixel fundamental property.

use serde::{Deserialize, Serialize};

/// Crustal material at a surface pixel.
///
/// Distinct from [`super::plates::PlateType`]: a plate may contain multiple crust types
/// (e.g. North America has continental interior and Atlantic ocean floor).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, Default)]
#[repr(u8)]
pub enum CrustType {
    /// Basaltic oceanic crust — dense, young, thin (~7 km).
    #[default]
    Oceanic,
    /// Granitic continental crust — buoyant, ancient (~35 km).
    Continental,
    /// Island arcs, rifted margins, transitional terranes.
    Transitional,
}

impl CrustType {
    /// Typical thickness scale (km) for isostasy models.
    #[must_use]
    pub const fn typical_thickness_km(self) -> f64 {
        match self {
            Self::Oceanic => 7.0,
            Self::Continental => 35.0,
            Self::Transitional => 20.0,
        }
    }

    /// Typical density (g/cm³) for buoyancy calculations.
    #[must_use]
    pub const fn typical_density_g_cm3(self) -> f64 {
        match self {
            Self::Oceanic => 3.0,
            Self::Continental => 2.7,
            Self::Transitional => 2.85,
        }
    }
}
