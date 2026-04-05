//! Map data structures and geographic utilities
//!
//! This module provides the core mapping infrastructure used throughout
//! the geoforge library, including terrain data structures, coordinate
//! systems, and world generation orchestration.

pub mod spherical;
pub mod tectonics_module;
pub mod terrain;
pub mod world;

// Re-export main types for convenience
pub use spherical::{EARTH_RADIUS_KM, EARTH_SURFACE_AREA_KM2, PlanetaryParams, SphericalPoint};
pub use tectonics_module::TectonicsModule;
pub use terrain::{
    BiomeMap, ElevationMap, MapProjection, MapStats, PlateMap, PrecipitationMap, TemperatureMap,
    TerrainMap,
};
pub use world::{TectonicMetadata, WorldMap};
