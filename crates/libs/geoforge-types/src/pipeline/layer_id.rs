//! Ordered pipeline layers for CLI and orchestration.

use serde::{Deserialize, Serialize};

/// Identifies a stage or output layer in the generation pipeline.
///
/// Ordering matches the intended generation sequence (lower = earlier).
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize,
)]
#[repr(u8)]
pub enum PipelineLayerId {
    /// Stage 0 — stellar/planetary system (future).
    PlanetarySystem = 0,
    /// Stage 1 — tectonic foundation.
    Tectonics = 10,
    /// Stage 2 — geologic provinces.
    Geology = 20,
    /// Stage 3 — elevation / bathymetry.
    Elevation = 30,
    /// Stage 4 — climate.
    Climate = 40,
    /// Stage 5 — biomes.
    Biomes = 50,
    /// Stage 6 — hydrology.
    Hydrology = 60,
    /// Stage 7+ — resources, hazards, settlements (future).
    Resources = 70,
}

impl PipelineLayerId {
    /// All layers in pipeline order.
    #[must_use]
    pub fn pipeline_order() -> &'static [Self] {
        &[
            Self::PlanetarySystem,
            Self::Tectonics,
            Self::Geology,
            Self::Elevation,
            Self::Climate,
            Self::Biomes,
            Self::Hydrology,
            Self::Resources,
        ]
    }

    /// Stable snake-case name for paths and CLI.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::PlanetarySystem => "planetary_system",
            Self::Tectonics => "tectonics",
            Self::Geology => "geology",
            Self::Elevation => "elevation",
            Self::Climate => "climate",
            Self::Biomes => "biomes",
            Self::Hydrology => "hydrology",
            Self::Resources => "resources",
        }
    }

    /// Parse from snake-case name.
    #[must_use]
    pub fn parse_str(s: &str) -> Option<Self> {
        Some(match s {
            "planetary_system" => Self::PlanetarySystem,
            "tectonics" => Self::Tectonics,
            "geology" => Self::Geology,
            "elevation" => Self::Elevation,
            "climate" => Self::Climate,
            "biomes" => Self::Biomes,
            "hydrology" => Self::Hydrology,
            "resources" => Self::Resources,
            _ => return None,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ordering_matches_pipeline() {
        let order = PipelineLayerId::pipeline_order();
        for window in order.windows(2) {
            assert!(window[0] < window[1]);
        }
    }
}
