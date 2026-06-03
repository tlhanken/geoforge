//! Morphology-aware galaxy regions for streaming (no full-array materialization).

use serde::{Deserialize, Serialize};

use super::morphology::{GalaxyMorphology, GalaxyProfile};
use crate::coordinates::{GalacticBounds, GalacticPoint, LoadSphere};

/// Partition cell identifier — scheme depends on galaxy morphology.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum RegionId {
    /// Thin disk / spiral / ring: annulus × azimuth × height slab.
    Cylindrical {
        /// Radial ring index (0 = innermost).
        ring: u16,
        /// Azimuthal wedge (0 .. wedge_count-1).
        wedge: u16,
        /// Vertical layer index.
        layer: u16,
    },
    /// Bubble / elliptical: shell × polar wedge × azimuthal wedge.
    Spherical {
        /// Radial shell index.
        shell: u16,
        /// Polar band (from north pole).
        theta_wedge: u16,
        /// Azimuthal wedge.
        phi_wedge: u16,
    },
}

/// Axis-aligned bounds approximating a region (ly), for culling and debug.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct RegionBounds {
    /// Conservative AABB in galactic coordinates.
    pub aabb: GalacticBounds,
}

/// List all region ids whose conservative bounds intersect a load sphere.
#[must_use]
pub fn regions_intersecting_load(profile: &GalaxyProfile, load: LoadSphere) -> Vec<RegionId> {
    let mut out = Vec::new();
    match profile.morphology {
        GalaxyMorphology::SpiralDisk | GalaxyMorphology::Ring | GalaxyMorphology::Irregular => {
            for ring in 0..profile.ring_count {
                for wedge in 0..profile.wedge_count {
                    for layer in 0..profile.depth_count {
                        let id = RegionId::Cylindrical {
                            ring,
                            wedge,
                            layer,
                        };
                        let bounds = cylindrical_region_bounds(profile, id);
                        if sphere_intersects_aabb(load, bounds.aabb) {
                            out.push(id);
                        }
                    }
                }
            }
        }
        GalaxyMorphology::Elliptical | GalaxyMorphology::Bubble => {
            for shell in 0..profile.ring_count {
                for theta in 0..profile.wedge_count {
                    for phi in 0..profile.depth_count {
                        let id = RegionId::Spherical {
                            shell,
                            theta_wedge: theta,
                            phi_wedge: phi,
                        };
                        let bounds = spherical_region_bounds(profile, id);
                        if sphere_intersects_aabb(load, bounds.aabb) {
                            out.push(id);
                        }
                    }
                }
            }
        }
    }
    out
}

/// Conservative AABB for a cylindrical region.
#[must_use]
pub fn cylindrical_region_bounds(profile: &GalaxyProfile, id: RegionId) -> RegionBounds {
    let RegionId::Cylindrical { ring, wedge, layer } = id else {
        return RegionBounds {
            aabb: GalacticBounds::from_center_half_extents(
                GalacticPoint::origin(),
                0.0,
                0.0,
                0.0,
            ),
        };
    };
    let rings = profile.ring_count.max(1) as f64;
    let wedges = profile.wedge_count.max(1) as f64;
    let layers = profile.depth_count.max(1) as f64;

    let r0 = profile.inner_radius_ly + (f64::from(ring) / rings) * profile.radius_ly;
    let r1 = profile.inner_radius_ly + (f64::from(ring + 1) / rings) * profile.radius_ly;
    let theta0 = (f64::from(wedge) / wedges) * std::f64::consts::TAU;
    let theta1 = (f64::from(wedge + 1) / wedges) * std::f64::consts::TAU;
    let z0 = (f64::from(layer) / layers - 0.5) * profile.scale_height_ly;
    let z1 = (f64::from(layer + 1) / layers - 0.5) * profile.scale_height_ly;

    let corners = [
        (r0, theta0, z0),
        (r1, theta1, z1),
    ];
    let mut min_x = f64::INFINITY;
    let mut min_y = f64::INFINITY;
    let mut min_z = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    let mut max_z = f64::NEG_INFINITY;

    for &(r, th, y) in &corners {
        let x = r * th.cos();
        let z = r * th.sin();
        min_x = min_x.min(x);
        min_y = min_y.min(y);
        min_z = min_z.min(z);
        max_x = max_x.max(x);
        max_y = max_y.max(y);
        max_z = max_z.max(z);
    }

    RegionBounds {
        aabb: GalacticBounds {
            min: GalacticPoint::new(min_x, min_y, min_z),
            max: GalacticPoint::new(max_x, max_y, max_z),
        },
    }
}

/// Conservative AABB for a spherical shell region.
#[must_use]
pub fn spherical_region_bounds(profile: &GalaxyProfile, id: RegionId) -> RegionBounds {
    let RegionId::Spherical {
        shell,
        theta_wedge,
        phi_wedge,
    } = id
    else {
        return cylindrical_region_bounds(profile, id);
    };
    let shells = profile.ring_count.max(1) as f64;
    let r0 = (f64::from(shell) / shells) * profile.radius_ly;
    let r1 = (f64::from(shell + 1) / shells) * profile.radius_ly;
    let half = r1.max(r0);
    RegionBounds {
        aabb: GalacticBounds::from_center_half_extents(
            GalacticPoint::origin(),
            half,
            half,
            half,
        ),
    }
}

/// True if `point` lies inside the galaxy's valid volume for its morphology.
#[must_use]
pub fn point_in_galaxy(profile: &GalaxyProfile, point: GalacticPoint) -> bool {
    match profile.morphology {
        GalaxyMorphology::SpiralDisk | GalaxyMorphology::Irregular => {
            let r = point.cylindrical_r_ly();
            r >= profile.inner_radius_ly
                && r <= profile.radius_ly
                && point.height_ly().abs() <= profile.scale_height_ly
        }
        GalaxyMorphology::Ring => {
            let r = point.cylindrical_r_ly();
            r >= profile.inner_radius_ly.max(profile.radius_ly * 0.3)
                && r <= profile.radius_ly
                && point.height_ly().abs() <= profile.scale_height_ly
        }
        GalaxyMorphology::Elliptical | GalaxyMorphology::Bubble => {
            let rx = point.x_ly / profile.radius_ly;
            let ry = point.y_ly / profile.scale_height_ly;
            let rz = point.z_ly / profile.radius_ly;
            (rx * rx + ry * ry + rz * rz) <= 1.0
        }
    }
}

fn sphere_intersects_aabb(load: LoadSphere, bounds: GalacticBounds) -> bool {
    let cx = load.center.x_ly.clamp(bounds.min.x_ly, bounds.max.x_ly);
    let cy = load.center.y_ly.clamp(bounds.min.y_ly, bounds.max.y_ly);
    let cz = load.center.z_ly.clamp(bounds.min.z_ly, bounds.max.z_ly);
    let dx = load.center.x_ly - cx;
    let dy = load.center.y_ly - cy;
    let dz = load.center.z_ly - cz;
    let dist_sq = dx * dx + dy * dy + dz * dz;
    dist_sq <= load.radius_ly * load.radius_ly
}
