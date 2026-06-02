//! Stellar system marker positions within a region.

use geoforge_types::cosmology::{
    point_in_galaxy, CosmologyConfig, GalaxyProfile, RegionId, StellarSystemMarker,
};
use geoforge_types::coordinates::GalacticPoint;
use geoforge_types::Seed;

use crate::density::systems_in_region;
use crate::phenotype::derive_phenotype;
use crate::rng::{draw, range};
use crate::seeds::{region_seed, system_seed};

/// Generate all markers in a region (bounded by config).
pub fn markers_in_region(
    _root: Seed,
    galaxy_index: u32,
    galaxy_seed: Seed,
    profile: &GalaxyProfile,
    region: RegionId,
    config: &CosmologyConfig,
    primary: Option<(RegionId, u32)>,
) -> Vec<StellarSystemMarker> {
    let count = systems_in_region(galaxy_seed, profile, region, config);
    let mut out = Vec::with_capacity(count as usize);

    for index_in_region in 0..count {
        let sys = system_seed(region_seed(galaxy_seed, region), index_in_region);
        let position_ly = marker_position(profile, region, index_in_region, sys);
        if !point_in_galaxy(profile, position_ly) {
            continue;
        }

        let force_g = config.preset == geoforge_types::cosmology::CosmologyPreset::Minimal
            && index_in_region == 0;
        let phenotype = derive_phenotype(sys, force_g);
        let is_primary = primary
            .map(|(r, i)| r == region && i == index_in_region)
            .unwrap_or(galaxy_index == 0 && index_in_region == 0 && region == primary_region(profile));

        out.push(StellarSystemMarker {
            galaxy_index,
            region,
            index_in_region,
            position_ly,
            phenotype,
            is_primary,
        });
    }

    out
}

fn primary_region(profile: &GalaxyProfile) -> RegionId {
    match profile.morphology {
        geoforge_types::cosmology::GalaxyMorphology::Elliptical
        | geoforge_types::cosmology::GalaxyMorphology::Bubble => RegionId::Spherical {
            shell: profile.ring_count / 2,
            theta_wedge: 0,
            phi_wedge: 0,
        },
        _ => RegionId::Cylindrical {
            ring: profile.ring_count / 3,
            wedge: 0,
            layer: 0_u16,
        },
    }
}

#[must_use]
pub(crate) fn marker_position(
    profile: &GalaxyProfile,
    region: RegionId,
    index: u32,
    sys: Seed,
) -> GalacticPoint {
    match (profile.morphology, region) {
        (
            _,
            RegionId::Cylindrical {
                ring,
                wedge,
                layer,
            },
        ) => {
            let rings = profile.ring_count.max(1) as f64;
            let wedges = profile.wedge_count.max(1) as f64;
            let layers = profile.depth_count.max(1) as f64;

            let r_inner = profile.inner_radius_ly + (f64::from(ring) / rings) * profile.radius_ly;
            let r_outer = profile.inner_radius_ly + (f64::from(ring + 1) / rings) * profile.radius_ly;
            let r = range(draw(sys, b"r"), r_inner, r_outer);

            let theta0 = (f64::from(wedge) / wedges) * std::f64::consts::TAU;
            let theta1 = (f64::from(wedge + 1) / wedges) * std::f64::consts::TAU;
            let theta = range(draw(sys, b"theta"), theta0, theta1);

            let z0 = (f64::from(layer) / layers - 0.5) * profile.scale_height_ly;
            let z1 = (f64::from(layer + 1) / layers - 0.5) * profile.scale_height_ly;
            let y = range(draw(sys, b"y"), z0.min(z1), z0.max(z1));

            GalacticPoint::new(r * theta.cos(), y, r * theta.sin())
        }
        (
            _,
            RegionId::Spherical {
                shell,
                theta_wedge,
                phi_wedge,
            },
        ) => {
            let shells = profile.ring_count.max(1) as f64;
            let r0 = (f64::from(shell) / shells) * profile.radius_ly;
            let r1 = (f64::from(shell + 1) / shells) * profile.radius_ly;
            let r = range(draw(sys, b"r"), r0, r1);
            let theta = range(draw(sys, b"theta"), 0.0, std::f64::consts::PI);
            let phi = range(
                draw(sys, b"phi"),
                (f64::from(phi_wedge) / f64::from(profile.depth_count.max(1)))
                    * std::f64::consts::TAU,
                (f64::from(phi_wedge + 1) / f64::from(profile.depth_count.max(1)))
                    * std::f64::consts::TAU,
            );
            GalacticPoint::new(
                r * theta.sin() * phi.cos(),
                r * theta.cos(),
                r * theta.sin() * phi.sin(),
            )
        }
    }
}
