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
pub mod motion;

pub use generator::{TectonicPlateGenerator, GenerationMethod};
pub use plates::{PlateSeed, PlateStats, PlateType, PlateInteraction, PlateBoundary};
pub use boundary_refinement::{BoundaryRefiner, BoundaryRefinementConfig};
pub use island_removal::{IslandRemover, IslandRemovalConfig, IslandRemovalStats};
pub use boundary_analysis::{BoundaryAnalyzer, BoundaryAnalysisConfig, BoundarySegment, BoundaryStatistics};
pub use motion::{PlateMotionAssigner, PlateMotionConfig};

pub use crate::error::GeoforgeError as PlateError;