//! Plate identifiers, seeds, and per-plate metadata.

use serde::{Deserialize, Serialize};

use crate::coordinates::{LatLon, PixelCoord};

/// Tectonic plate identifier (`0` = unassigned).
pub type PlateId = u16;

/// Dominant character of a tectonic plate (metadata, not per-pixel crust).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, Default)]
pub enum PlateType {
    /// Primarily oceanic lithosphere.
    #[default]
    Oceanic,
    /// Primarily continental lithosphere.
    Continental,
}

impl PlateType {
    /// Classify from area percentile among all plates (Earth-like heuristic).
    #[must_use]
    pub fn from_size_percentile(percentile: f64) -> Self {
        if percentile >= 0.5 {
            Self::Continental
        } else {
            Self::Oceanic
        }
    }

    /// Default crust type when initializing pixels on this plate.
    #[must_use]
    pub const fn default_crust_type(self) -> super::CrustType {
        match self {
            Self::Oceanic => super::CrustType::Oceanic,
            Self::Continental => super::CrustType::Continental,
        }
    }
}

/// Motion of a plate at its seed point.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PlateMotion {
    /// Azimuth in degrees (0 = east, 90 = north).
    pub direction_deg: f64,
    /// Speed in cm/year.
    pub speed_cm_per_year: f64,
}

/// Seed point for a tectonic plate (Voronoi generator site).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PlateSeed {
    /// Plate id.
    pub id: PlateId,
    /// Seed pixel location.
    pub pixel: PixelCoord,
    /// Geographic position of seed.
    pub latlon: LatLon,
    /// Assigned motion.
    pub motion: PlateMotion,
}

/// Aggregated statistics for one plate.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PlateStats {
    /// Pixel count.
    pub pixels: usize,
    /// Fraction of planetary surface (%).
    pub percentage: f64,
    /// Surface area (km²).
    pub area_km2: f64,
    /// Seed and motion.
    pub seed: PlateSeed,
    /// Dominant plate character.
    pub plate_type: PlateType,
}
