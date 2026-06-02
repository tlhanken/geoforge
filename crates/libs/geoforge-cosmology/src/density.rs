//! Systems per region from seed and preset.

use geoforge_types::cosmology::{CosmologyConfig, CosmologyPreset, GalaxyProfile, RegionId};
use geoforge_types::Seed;

use crate::rng::{draw, unit};
use crate::seeds::region_seed;

/// Number of stellar systems to generate in a region.
#[must_use]
pub fn systems_in_region(
    galaxy_seed: Seed,
    profile: &GalaxyProfile,
    region: RegionId,
    config: &CosmologyConfig,
) -> u32 {
    match config.preset {
        CosmologyPreset::Minimal => return 1,
        CosmologyPreset::Rich | CosmologyPreset::Expansive => {}
    }

    let rseed = region_seed(galaxy_seed, region);
    let base = match profile.morphology {
        geoforge_types::cosmology::GalaxyMorphology::SpiralDisk => {
            let RegionId::Cylindrical { ring, .. } = region else {
                return 0;
            };
            let t = f64::from(ring) / f64::from(profile.ring_count.max(1));
            (8.0 * (-t * 2.0).exp() + 2.0) as u32
        }
        geoforge_types::cosmology::GalaxyMorphology::Bubble
        | geoforge_types::cosmology::GalaxyMorphology::Elliptical => 4,
        geoforge_types::cosmology::GalaxyMorphology::Ring => 3,
        geoforge_types::cosmology::GalaxyMorphology::Irregular => 2,
    };

    let jitter = (unit(draw(rseed, b"count")) * 4.0) as u32;
    let n = base + jitter;
    n.min(config.max_systems_per_region)
}
