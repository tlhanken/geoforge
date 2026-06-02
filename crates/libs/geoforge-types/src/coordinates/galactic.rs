//! Galactic-scale positions for maps and streaming regions.

use serde::{Deserialize, Serialize};

use super::{DistanceContext, DistanceKm};

/// Position in galactic Cartesian coordinates (light-years from galaxy center).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct GalacticPoint {
    /// X component (ly).
    pub x_ly: f64,
    /// Y component (ly) — often treated as disk axis for spirals.
    pub y_ly: f64,
    /// Z component (ly).
    pub z_ly: f64,
}

impl GalacticPoint {
    /// Origin at the galactic center.
    #[must_use]
    pub const fn origin() -> Self {
        Self {
            x_ly: 0.0,
            y_ly: 0.0,
            z_ly: 0.0,
        }
    }

    /// Construct from components.
    #[must_use]
    pub const fn new(x_ly: f64, y_ly: f64, z_ly: f64) -> Self {
        Self { x_ly, y_ly, z_ly }
    }

    /// Cylindrical radius in the x–z plane (ly), assuming y is vertical axis.
    #[must_use]
    pub fn cylindrical_r_ly(self) -> f64 {
        (self.x_ly * self.x_ly + self.z_ly * self.z_ly).sqrt()
    }

    /// Azimuth in radians in the x–z plane, range `(-π, π]`.
    #[must_use]
    pub fn cylindrical_theta_rad(self) -> f64 {
        self.z_ly.atan2(self.x_ly)
    }

    /// Height along disk axis (ly).
    #[must_use]
    pub const fn height_ly(self) -> f64 {
        self.y_ly
    }

    /// Spherical radius from origin (ly).
    #[must_use]
    pub fn spherical_r_ly(self) -> f64 {
        (self.x_ly * self.x_ly + self.y_ly * self.y_ly + self.z_ly * self.z_ly).sqrt()
    }

    /// Offset by a vector (ly).
    #[must_use]
    pub fn offset(self, dx: f64, dy: f64, dz: f64) -> Self {
        Self {
            x_ly: self.x_ly + dx,
            y_ly: self.y_ly + dy,
            z_ly: self.z_ly + dz,
        }
    }
}

impl DistanceKm for GalacticPoint {
    fn distance_km(&self, other: &Self, _ctx: &DistanceContext) -> f64 {
        const LY_KM: f64 = 9.461e12;
        let dx = self.x_ly - other.x_ly;
        let dy = self.y_ly - other.y_ly;
        let dz = self.z_ly - other.z_ly;
        (dx * dx + dy * dy + dz * dz).sqrt() * LY_KM
    }
}

/// Axis-aligned bounding box in galactic space (ly).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct GalacticBounds {
    /// Minimum corner (ly).
    pub min: GalacticPoint,
    /// Maximum corner (ly).
    pub max: GalacticPoint,
}

impl GalacticBounds {
    /// Bounds from center and half-extents.
    #[must_use]
    pub fn from_center_half_extents(center: GalacticPoint, half_x: f64, half_y: f64, half_z: f64) -> Self {
        Self {
            min: GalacticPoint::new(
                center.x_ly - half_x,
                center.y_ly - half_y,
                center.z_ly - half_z,
            ),
            max: GalacticPoint::new(
                center.x_ly + half_x,
                center.y_ly + half_y,
                center.z_ly + half_z,
            ),
        }
    }

    /// True if the point lies inside (inclusive).
    #[must_use]
    pub fn contains(self, p: GalacticPoint) -> bool {
        p.x_ly >= self.min.x_ly
            && p.x_ly <= self.max.x_ly
            && p.y_ly >= self.min.y_ly
            && p.y_ly <= self.max.y_ly
            && p.z_ly >= self.min.z_ly
            && p.z_ly <= self.max.z_ly
    }
}

/// Sphere used for local region loading (ly).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct LoadSphere {
    /// Center (ly).
    pub center: GalacticPoint,
    /// Radius (ly).
    pub radius_ly: f64,
}
