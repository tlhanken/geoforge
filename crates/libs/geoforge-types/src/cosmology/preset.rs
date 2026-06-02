//! Generation scale presets — control region load volume, not persistence.

use serde::{Deserialize, Serialize};

/// How much of the galaxy hierarchy to materialize per query.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, Default)]
pub enum CosmologyPreset {
    /// One galaxy, one region, one primary system.
    #[default]
    Minimal,
    /// Local neighborhood: multiple regions around focus, catalog density.
    Rich,
    /// Higher caps for stress / visualization.
    Expansive,
}

/// Tunable cosmology generation parameters.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CosmologyConfig {
    /// Scale preset.
    pub preset: CosmologyPreset,
    /// Number of galaxies (1..=max_galaxies).
    pub galaxy_count: u32,
    /// Radius around focus to load regions (ly).
    pub load_radius_ly: f64,
    /// Hard cap on markers per region.
    pub max_systems_per_region: u32,
    /// Hard cap on total markers per query.
    pub max_markers_per_query: u32,
}

impl CosmologyConfig {
    /// Config for [`CosmologyPreset::Minimal`].
    #[must_use]
    pub fn minimal() -> Self {
        Self {
            preset: CosmologyPreset::Minimal,
            galaxy_count: 1,
            load_radius_ly: 500.0,
            max_systems_per_region: 1,
            max_markers_per_query: 1,
        }
    }

    /// Config for [`CosmologyPreset::Rich`].
    #[must_use]
    pub fn rich() -> Self {
        Self {
            preset: CosmologyPreset::Rich,
            galaxy_count: 1,
            load_radius_ly: 2_000.0,
            max_systems_per_region: 256,
            max_markers_per_query: 4_096,
        }
    }

    /// Config for [`CosmologyPreset::Expansive`].
    #[must_use]
    pub fn expansive() -> Self {
        Self {
            preset: CosmologyPreset::Expansive,
            galaxy_count: 3,
            load_radius_ly: 5_000.0,
            max_systems_per_region: 1_024,
            max_markers_per_query: 32_768,
        }
    }

    /// Build from preset with defaults.
    #[must_use]
    pub fn from_preset(preset: CosmologyPreset) -> Self {
        match preset {
            CosmologyPreset::Minimal => Self::minimal(),
            CosmologyPreset::Rich => Self::rich(),
            CosmologyPreset::Expansive => Self::expansive(),
        }
    }
}
