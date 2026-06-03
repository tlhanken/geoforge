//! Galaxy instance metadata.

use geoforge_types::cosmology::{CosmologyConfig, CosmologyPreset, GalaxyInstance, GalaxyProfile};
use geoforge_types::Seed;

use crate::rng::unit;
use crate::seeds::galaxy_seed;

/// Build galaxy instance for index `g`.
#[must_use]
pub fn galaxy_instance(root: Seed, index: u32, config: &CosmologyConfig) -> GalaxyInstance {
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
    let label = format!("Galaxy-{:04}", (unit(gseed) * 10_000.0) as u32);

    GalaxyInstance {
        index,
        label,
        profile,
    }
}

/// Galaxy count for this config.
#[must_use]
pub fn galaxy_count(config: &CosmologyConfig) -> u32 {
    match config.preset {
        CosmologyPreset::Minimal => 1,
        CosmologyPreset::Rich => 1,
        CosmologyPreset::Expansive => config.galaxy_count.min(3).max(1),
    }
}
