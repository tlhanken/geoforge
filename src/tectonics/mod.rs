//! Tectonic plate simulation and generation
//! 
//! This module handles the generation of realistic tectonic plates using various
//! algorithms including Voronoi diagrams and region growing on the sphere.

pub mod generator;
pub mod plates;
pub mod algorithms;

pub use generator::TectonicPlateGenerator;
pub use plates::{PlateSeed, PlateStats};


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