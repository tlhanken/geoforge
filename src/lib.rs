//! Geoforge - Realistic geological and climate modeling for procedural world generation
//! 
//! This library provides tools for generating realistic geological features, climate patterns,
//! and biomes for procedural world generation. It starts with tectonic plate simulation
//! and builds up through geological domains, elevation, climate, and biomes.
//! 
//! # Quick Start
//! 
//! ```rust
//! use geoforge::TectonicPlateGenerator;
//!
//! // Generate a world with 15 tectonic plates at 0.2° resolution
//! let mut generator = TectonicPlateGenerator::with_seed(1800, 900, 15, 42)?;
//! let plate_map = generator.generate("electrostatic")?;
//!
//! // Get statistics about the generated plates
//! let stats = generator.get_plate_stats();
//! for (plate_id, stat) in stats {
//!     println!("Plate {}: {:.1}% of surface", plate_id, stat.percentage);
//! }
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

// Core modules
pub mod map;
pub mod tectonics;

// Re-export main types for convenience
pub use map::{TerrainMap, PlateMap, ElevationMap, TemperatureMap, PrecipitationMap, BiomeMap};
pub use map::{SphericalPoint, PlanetaryParams, EARTH_RADIUS_KM, EARTH_SURFACE_AREA_KM2};
pub use map::{WorldMap, TectonicMetadata};
pub use tectonics::{TectonicPlateGenerator, PlateSeed, PlateStats, PlateType, PlateError};
pub use tectonics::{PlateInteraction, PlateBoundary};
pub use tectonics::{BoundaryRefiner, BoundaryRefinementConfig};
pub use tectonics::{IslandRemover, IslandRemovalConfig, IslandRemovalStats};
pub use tectonics::{BoundaryAnalyzer, BoundaryAnalysisConfig, BoundarySegment, BoundaryStatistics};
