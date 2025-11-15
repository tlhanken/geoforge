//! Map data structures and geographic utilities
//! 
//! This module provides the core mapping infrastructure used throughout
//! the geoforge library, including terrain data structures, coordinate
//! systems, and world generation orchestration.

pub mod terrain;
pub mod world;
pub mod spherical;

// Re-export main types for convenience
pub use terrain::{
    TerrainMap, MapProjection, MapStats,
    PlateMap, ElevationMap, TemperatureMap, PrecipitationMap, BiomeMap
};
pub use world::{WorldMap, TectonicMetadata};
pub use spherical::{SphericalPoint, PlanetaryParams, EARTH_RADIUS_KM, EARTH_SURFACE_AREA_KM2};