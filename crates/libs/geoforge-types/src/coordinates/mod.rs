//! Coordinate systems used across the pipeline.

mod euclidean;
mod galactic;
mod orbital;
mod planetary;

pub use euclidean::EuclideanPoint;
pub use galactic::{GalacticBounds, GalacticPoint, LoadSphere};
pub use orbital::OrbitalElements;
pub use planetary::{LatLon, MapExtent, PixelCoord};

/// Trait for coordinates that support geodesic or Euclidean distance in kilometers.
pub trait DistanceKm {
    /// Distance to another point in kilometers.
    fn distance_km(&self, other: &Self, context: &DistanceContext) -> f64;
}

/// Context required to interpret distances (e.g. planetary radius).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DistanceContext {
    /// Body radius in kilometers.
    pub radius_km: f64,
}

impl DistanceContext {
    /// Earth-like radius (WGS84 mean).
    #[must_use]
    pub const fn earth() -> Self {
        Self {
            radius_km: crate::planetary::EARTH_RADIUS_KM,
        }
    }
}
