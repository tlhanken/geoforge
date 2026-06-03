//! Stable references for on-demand regeneration.

use serde::{Deserialize, Serialize};

use super::region::RegionId;

/// Points to one stellar system in a galaxy (enough to re-derive seeds).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct StellarSystemRef {
    /// Host galaxy index.
    pub galaxy_index: u32,
    /// Region id.
    pub region: RegionId,
    /// Index within region.
    pub index_in_region: u32,
}

/// Summary of one galaxy instance in a cosmology run.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GalaxyInstance {
    /// Index in the run (0 = primary).
    pub index: u32,
    /// Display name seed-derived label.
    pub label: String,
    /// Structural profile.
    pub profile: super::morphology::GalaxyProfile,
}
