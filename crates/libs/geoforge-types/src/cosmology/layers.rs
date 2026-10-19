//! Aggregated outputs for cosmology pipeline stages 1–5 (pre-tectonics).

use serde::{Deserialize, Serialize};

use super::{GalaxyCatalog, SolarSystem, Supercluster, Universe};
use crate::planetary::PlanetaryBodies;
use crate::Seed;

/// Layer bundle from universe through planetary body (before tectonics).
///
/// Data only — run stages via `geoforge_export::run_pre_tectonic_layers` or call each
/// `stageN_*::generate` in order. Query methods live on this type in `geoforge-export`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CosmologyLayers {
    /// Stage 1 — universe constants and config (`universe.root` is the run seed).
    pub universe: Universe,
    /// Stage 2 — supercluster framing.
    pub supercluster: Supercluster,
    /// Stage 3 — galaxy catalog and map markers.
    pub catalog: GalaxyCatalog,
    /// Stage 4 — primary solar system zoom (when generated).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub solar_system: Option<SolarSystem>,
    /// Stage 5 — planetary parameters per slot (when generated).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub planetary: Option<PlanetaryBodies>,
}

impl CosmologyLayers {
    /// Root procedural seed.
    #[must_use]
    pub fn root(&self) -> Seed {
        self.universe.root
    }

    /// Active cosmology config.
    #[must_use]
    pub fn config(&self) -> &super::preset::CosmologyConfig {
        &self.universe.config
    }

    /// Load focus for region queries (ly).
    #[must_use]
    pub fn focus_ly(&self) -> crate::coordinates::GalacticPoint {
        self.supercluster.focus_ly
    }
}
