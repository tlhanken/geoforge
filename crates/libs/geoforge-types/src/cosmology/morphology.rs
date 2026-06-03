//! Galaxy morphology and structural parameters.

use serde::{Deserialize, Serialize};

/// High-level galaxy classification — drives region partition scheme.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum GalaxyMorphology {
    /// Thin disk with optional spiral arms (cylindrical regions).
    SpiralDisk,
    /// 3D ellipsoidal mass distribution (spherical shell regions).
    Elliptical,
    /// Spheroidal / dwarf spheroid (spherical shells).
    Bubble,
    /// Disk with central hole (cylindrical, `r_min > 0`).
    Ring,
    /// No strong symmetry — cylindrical regions with loose validity mask.
    Irregular,
}

/// Structural parameters for partitioning and density (real-scale ly).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GalaxyProfile {
    /// Morphology class.
    pub morphology: GalaxyMorphology,
    /// Effective radius containing most stellar systems (ly).
    pub radius_ly: f64,
    /// Disk scale height or ellipsoid semi-axis b (ly).
    pub scale_height_ly: f64,
    /// For [`GalaxyMorphology::Ring`]: inner hole radius (ly).
    pub inner_radius_ly: f64,
    /// Number of radial rings (cylindrical) or shells (spherical).
    pub ring_count: u16,
    /// Azimuthal wedges per ring (θ) or per shell.
    pub wedge_count: u16,
    /// Vertical slabs (cylindrical z) or φ wedges (spherical); use 1 for thin 2D maps.
    pub depth_count: u16,
}

impl GalaxyProfile {
    /// Milky-Way-like spiral disk (order-of-magnitude scales).
    #[must_use]
    pub fn spiral_disk_default() -> Self {
        Self {
            morphology: GalaxyMorphology::SpiralDisk,
            radius_ly: 50_000.0,
            scale_height_ly: 1_000.0,
            inner_radius_ly: 0.0,
            ring_count: 64,
            wedge_count: 48,
            depth_count: 5,
        }
    }

    /// Compact spheroidal galaxy.
    #[must_use]
    pub fn bubble_default() -> Self {
        Self {
            morphology: GalaxyMorphology::Bubble,
            radius_ly: 15_000.0,
            scale_height_ly: 15_000.0,
            inner_radius_ly: 0.0,
            ring_count: 24,
            wedge_count: 16,
            depth_count: 16,
        }
    }

    /// Profile for minimal preset (single region, tiny extent).
    #[must_use]
    pub fn minimal() -> Self {
        Self {
            morphology: GalaxyMorphology::SpiralDisk,
            radius_ly: 1_000.0,
            scale_height_ly: 200.0,
            inner_radius_ly: 0.0,
            ring_count: 1,
            wedge_count: 1,
            depth_count: 1,
        }
    }
}
