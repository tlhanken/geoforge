//! Plate boundary classification.

use serde::{Deserialize, Serialize};

use super::PlateId;

/// Index into a boundary segment list.
pub type BoundarySegmentId = usize;

/// Relative motion classification at a plate boundary.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum BoundaryType {
    /// Plates moving toward each other.
    Convergent,
    /// Plates moving apart.
    Divergent,
    /// Plates sliding past one another.
    Transform,
    /// No boundary / interior pixel.
    None,
}

/// Geologic style of interaction (refinement of convergent boundaries).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum PlateBoundaryKind {
    /// Oceanic–oceanic convergence.
    OceanicOceanic,
    /// Oceanic subducts under continental.
    OceanicContinental,
    /// Continental–continental collision.
    ContinentalContinental,
    /// Oceanic spreading.
    OceanicDivergent,
    /// Continental rifting.
    ContinentalDivergent,
    /// Strike-slip.
    Transform,
}

/// Summary statistics for all boundaries in a world.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, Default)]
pub struct BoundaryStatistics {
    /// Count of convergent segments.
    pub convergent: usize,
    /// Count of divergent segments.
    pub divergent: usize,
    /// Count of transform segments.
    pub transform: usize,
    /// Total boundary length (km), if measured.
    pub total_length_km: Option<f64>,
}

/// Metadata for one boundary between two plates (not per-pixel storage).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BoundarySegmentMeta {
    /// First plate id.
    pub plate_a: PlateId,
    /// Second plate id.
    pub plate_b: PlateId,
    /// Classified boundary type.
    pub boundary_type: BoundaryType,
    /// Geologic interaction style.
    pub kind: PlateBoundaryKind,
    /// Relative convergence rate (cm/yr); negative = divergent.
    pub rate_cm_per_year: f64,
}
