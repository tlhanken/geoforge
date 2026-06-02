//! Galaxy-map stellar system markers (no planetary interiors).

use serde::{Deserialize, Serialize};

use super::region::RegionId;
use super::star::SpectralClass;
use crate::coordinates::GalacticPoint;

/// Stellar multiplicity at map LOD (derived from system seed).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Multiplicity {
    /// Single star.
    Single,
    /// Binary.
    Binary,
    /// Triple or more.
    Trinary,
}

/// Compact summary for galaxy-map rendering (same seed as full system).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SystemPhenotype {
    /// Single / binary / trinary.
    pub multiplicity: Multiplicity,
    /// Primary spectral class.
    pub primary_class: SpectralClass,
    /// Companion class if binary+.
    pub secondary_class: Option<SpectralClass>,
}

/// A stellar system placeholder on the galaxy map.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StellarSystemMarker {
    /// Galaxy index (0 = primary host).
    pub galaxy_index: u32,
    /// Region containing this system.
    pub region: RegionId,
    /// Index within the region (stable for subseeds).
    pub index_in_region: u32,
    /// Absolute position for rendering (ly).
    pub position_ly: GalacticPoint,
    /// Map icon / color hints.
    pub phenotype: SystemPhenotype,
    /// True for the world's primary simulation target.
    pub is_primary: bool,
}
