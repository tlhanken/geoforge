//! Tectonic plate simulation and generation
//!
//! This module handles the generation of realistic tectonic plates using
//! electrostatic physics simulation on the sphere for natural boundaries
//! and Earth-like size variety.

pub mod boundary_analysis;
pub mod boundary_refinement;
pub mod electrostatic;
pub mod generator;
pub mod island_removal;
pub mod motion;
pub mod plates;

pub use boundary_analysis::{
    BoundaryAnalysisConfig, BoundaryAnalyzer, BoundarySegment, BoundaryStatistics,
};
pub use boundary_refinement::{BoundaryRefinementConfig, BoundaryRefiner};
pub use generator::{GenerationMethod, TectonicPlateGenerator};
pub use island_removal::{IslandRemovalConfig, IslandRemovalStats, IslandRemover};
pub use motion::{PlateMotionAssigner, PlateMotionConfig};
pub use plates::{PlateBoundary, PlateInteraction, PlateSeed, PlateStats, PlateType};

pub use crate::error::GeoforgeError as PlateError;
