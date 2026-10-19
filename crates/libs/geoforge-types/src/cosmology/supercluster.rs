//! Stage 2 output — local supercluster framing in galactic coordinates.

use serde::{Deserialize, Serialize};

use super::preset::CosmologyConfig;
use crate::coordinates::GalacticPoint;
use crate::Seed;

/// Supercluster-scale layer (pipeline stage 2).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Supercluster {
    /// Root seed (from stage 1).
    pub root: Seed,
    /// Cosmology config (from stage 1).
    pub config: CosmologyConfig,
    /// Display label.
    pub label: String,
    /// Nominal extent of the loaded supercluster volume (ly).
    pub extent_ly: f64,
    /// Origin offset for nested structure (ly).
    pub origin_ly: GalacticPoint,
    /// Default focus for region loading (ly).
    pub focus_ly: GalacticPoint,
}
