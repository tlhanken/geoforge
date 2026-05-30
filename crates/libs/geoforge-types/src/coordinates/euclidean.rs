//! Euclidean coordinates for cosmological scales (light-years, megaparsecs as raw f64).

use serde::{Deserialize, Serialize};

use crate::coordinates::{DistanceContext, DistanceKm};

/// Position in 3D Euclidean space (units defined by caller).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct EuclideanPoint {
    /// X coordinate.
    pub x: f64,
    /// Y coordinate.
    pub y: f64,
    /// Z coordinate.
    pub z: f64,
}

impl EuclideanPoint {
    /// Euclidean distance to another point (same units as coordinates).
    #[must_use]
    pub fn distance(&self, other: &Self) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        let dz = self.z - other.z;
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}

impl DistanceKm for EuclideanPoint {
    /// Interprets coordinates as kilometers via `DistanceContext.radius_km` scale factor of 1.
    ///
    /// For cosmological use, callers should pass an appropriate effective radius or use
    /// `distance()` directly with their unit convention.
    fn distance_km(&self, other: &Self, _ctx: &DistanceContext) -> f64 {
        self.distance(other)
    }
}
