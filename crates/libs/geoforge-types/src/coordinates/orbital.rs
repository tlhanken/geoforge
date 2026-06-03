//! Keplerian orbital elements for moons and planets (Stage 0).

use serde::{Deserialize, Serialize};

/// Orbital elements for a body around a parent (star or planet).
///
/// Angles are in degrees; distances in AU unless noted in field docs.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct OrbitalElements {
    /// Orbital period in Earth days.
    pub period_days: f64,
    /// Eccentricity (0 = circle).
    pub eccentricity: f64,
    /// Semi-major axis in astronomical units.
    pub semi_major_axis_au: f64,
    /// Inclination relative to reference plane (degrees).
    pub inclination_deg: f64,
    /// Current orbital phase in `[0, 1)` (0 = periapsis).
    pub phase: f64,
}

impl OrbitalElements {
    /// Earth-like heliocentric orbit.
    #[must_use]
    pub fn earth_like() -> Self {
        Self {
            period_days: 365.25,
            eccentricity: 0.017,
            semi_major_axis_au: 1.0,
            inclination_deg: 0.0,
            phase: 0.0,
        }
    }

    /// Periapsis distance in AU.
    #[must_use]
    pub fn periapsis_au(self) -> f64 {
        self.semi_major_axis_au * (1.0 - self.eccentricity)
    }

    /// Apoapsis distance in AU.
    #[must_use]
    pub fn apoapsis_au(self) -> f64 {
        self.semi_major_axis_au * (1.0 + self.eccentricity)
    }
}
