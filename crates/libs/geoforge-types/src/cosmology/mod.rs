//! Cosmology types: galaxy regions, map markers, solar-system zoom.

mod marker;
mod morphology;
mod preset;
mod reference;
mod region;
mod solar_system;
mod star;

pub use marker::{Multiplicity, StellarSystemMarker, SystemPhenotype};
pub use morphology::{GalaxyMorphology, GalaxyProfile};
pub use preset::{CosmologyConfig, CosmologyPreset};
pub use reference::{GalaxyInstance, StellarSystemRef};
pub use region::{
    cylindrical_region_bounds, point_in_galaxy, regions_intersecting_load, RegionBounds, RegionId,
    spherical_region_bounds,
};
pub use solar_system::{BeltSlot, PlanetMassClass, PlanetSlot, SolarSystem};
pub use star::{BarycentricOffset, SpectralClass, Star};
