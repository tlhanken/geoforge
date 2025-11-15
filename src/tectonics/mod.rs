//! Tectonic plate simulation and generation
//!
//! This module handles the generation of realistic tectonic plates using
//! electrostatic physics simulation on the sphere for natural boundaries
//! and Earth-like size variety.

pub mod generator;
pub mod plates;
pub mod electrostatic;
pub mod boundary_refinement;
pub mod island_removal;
pub mod boundary_analysis;

pub use generator::TectonicPlateGenerator;
pub use plates::{PlateSeed, PlateStats, PlateType, PlateInteraction, PlateBoundary};
pub use boundary_refinement::{BoundaryRefiner, BoundaryRefinementConfig};
pub use island_removal::{IslandRemover, IslandRemovalConfig, IslandRemovalStats};
pub use boundary_analysis::{BoundaryAnalyzer, BoundaryAnalysisConfig, BoundarySegment, BoundaryStatistics};


/// Errors that can occur during plate generation
#[derive(Debug)]
pub enum PlateError {
    SeedPlacementFailed,
    InvalidMethod(String),
    InvalidParameters(String),
}

impl std::fmt::Display for PlateError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PlateError::SeedPlacementFailed => write!(f, "Failed to place seeds with minimum distance"),
            PlateError::InvalidMethod(method) => write!(f, "Invalid generation method: {}", method),
            PlateError::InvalidParameters(msg) => write!(f, "Invalid parameters: {}", msg),
        }
    }
}

impl std::error::Error for PlateError {}