//! Pipeline stage 3 — galaxy catalog from supercluster.

#![deny(clippy::todo, clippy::unimplemented)]

mod density;
mod marker;
mod phenotype;

pub use marker::{marker_position, markers_in_region};

use geoforge_types::coordinates::LoadSphere;
use geoforge_types::cosmology::galaxy_seed;
use geoforge_types::cosmology::{
    regions_intersecting_load, GalaxyCatalog, GalaxyInstance, GalaxyMorphology, GalaxyProfile,
    StellarSystemRef, Supercluster,
};
use geoforge_types::cosmology::CosmologyPreset;
use geoforge_types::seed::{draw, unit};
use geoforge_types::Seed;

/// Number of galaxies for this config.
#[must_use]
pub fn galaxy_count(config: &geoforge_types::cosmology::CosmologyConfig) -> u32 {
    match config.preset {
        CosmologyPreset::Minimal => 1,
        CosmologyPreset::Rich => 1,
        CosmologyPreset::Expansive => config.galaxy_count.min(3).max(1),
    }
}

/// Build one galaxy instance.
#[must_use]
pub fn galaxy_instance(
    root: Seed,
    index: u32,
    config: &geoforge_types::cosmology::CosmologyConfig,
) -> GalaxyInstance {
    let profile = match config.preset {
        CosmologyPreset::Minimal => GalaxyProfile::minimal(),
        CosmologyPreset::Rich => GalaxyProfile::spiral_disk_default(),
        CosmologyPreset::Expansive => {
            if index == 0 {
                GalaxyProfile::spiral_disk_default()
            } else if index == 1 {
                GalaxyProfile::bubble_default()
            } else {
                GalaxyProfile::spiral_disk_default()
            }
        }
    };

    let gseed = galaxy_seed(root, index);
    let label = format!("Galaxy-{:04}", (unit(draw(gseed, b"label")) * 10_000.0) as u32);

    GalaxyInstance {
        index,
        label,
        profile,
    }
}

/// Deterministic primary system reference for galaxy 0.
#[must_use]
pub fn primary_ref(galaxy: &GalaxyInstance) -> StellarSystemRef {
    let region = match galaxy.profile.morphology {
        GalaxyMorphology::Elliptical | GalaxyMorphology::Bubble => {
            geoforge_types::cosmology::RegionId::Spherical {
                shell: galaxy.profile.ring_count / 2,
                theta_wedge: 0,
                phi_wedge: 0,
            }
        }
        _ => geoforge_types::cosmology::RegionId::Cylindrical {
            ring: galaxy.profile.ring_count / 3,
            wedge: 0,
            layer: 0_u16,
        },
    };
    StellarSystemRef {
        galaxy_index: 0,
        region,
        index_in_region: 0,
    }
}

/// Generate galaxy catalog and markers near supercluster focus.
#[must_use]
pub fn generate(supercluster: &Supercluster) -> GalaxyCatalog {
    let count = galaxy_count(&supercluster.config);
    let galaxies: Vec<GalaxyInstance> = (0..count)
        .map(|i| galaxy_instance(supercluster.root, i, &supercluster.config))
        .collect();

    let primary_ref = primary_ref(galaxies.first().expect("at least one galaxy"));
    let primary_key = Some((primary_ref.region, primary_ref.index_in_region));

    let load = LoadSphere {
        center: supercluster.focus_ly,
        radius_ly: supercluster.config.load_radius_ly,
    };

    let mut markers_near_focus = Vec::new();
    for g in &galaxies {
        let regions = regions_intersecting_load(&g.profile, load);
        let gseed = galaxy_seed(supercluster.root, g.index);
        for region in regions {
            for m in marker::markers_in_region(
                g.index,
                gseed,
                &g.profile,
                region,
                &supercluster.config,
                primary_key,
            ) {
                markers_near_focus.push(m);
                if markers_near_focus.len() as u32 >= supercluster.config.max_markers_per_query
                {
                    return GalaxyCatalog {
                        root: supercluster.root,
                        galaxies,
                        primary_ref,
                        markers_near_focus,
                    };
                }
            }
        }
    }

    GalaxyCatalog {
        root: supercluster.root,
        galaxies,
        primary_ref,
        markers_near_focus,
    }
}
