//! Planetary-scale physical parameters.

mod body;
mod params;

pub use body::{
    PlanetBodyRecord, PlanetaryBodies, DEFAULT_MAP_HEIGHT, DEFAULT_MAP_WIDTH,
};
pub use params::{PlanetaryParams, EARTH_RADIUS_KM, EARTH_SURFACE_AREA_KM2};
