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
    /// Stage 0 — universe constants (future).
    Universe = 0,
    /// Stage 0 — optional supercluster label.
    Supercluster = 5,
    /// Stage 0 — galaxy structure and region map.
    Galaxy = 8,
    /// Stage 0 — stellar system (stars + orbit slots).
    SolarSystem = 15,
    /// Stage 0.5 — full planetary body (surfaces; future).
    PlanetaryBody = 18,
    /// Stage 1 — tectonic foundation.
    Tectonics = 20,
    /// Stage 2 — geologic provinces.
    Geology = 30,
    /// Stage 3 — elevation / bathymetry.
    Elevation = 40,
    /// Stage 4 — climate.
    Climate = 50,
    /// Stage 5 — biomes.
    Biomes = 60,
    /// Stage 6 — hydrology.
    Hydrology = 70,
    /// Stage 7+ — resources (future).
    Resources = 80,
}

impl PipelineLayerId {
    /// All layers in pipeline order.
    #[must_use]
    pub fn pipeline_order() -> &'static [Self] {
        &[
            Self::Universe,
            Self::Supercluster,
            Self::Galaxy,
            Self::SolarSystem,
            Self::PlanetaryBody,
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
            Self::Universe => "universe",
            Self::Supercluster => "supercluster",
            Self::Galaxy => "galaxy",
            Self::SolarSystem => "solar_system",
            Self::PlanetaryBody => "planetary_body",
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
            "universe" => Self::Universe,
            "supercluster" => Self::Supercluster,
            "galaxy" => Self::Galaxy,
            "solar_system" | "planetary_system" => Self::SolarSystem,
            "planetary_body" => Self::PlanetaryBody,
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
