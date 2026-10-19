//! Cosmology types: galaxy regions, map markers, solar-system zoom.

mod catalog;
mod layers;
mod marker;
mod morphology;
mod preset;
mod reference;
mod region;
mod seeds;
mod solar_system;
mod star;
mod supercluster;
mod universe;

pub use catalog::GalaxyCatalog;
pub use layers::CosmologyLayers;
pub use marker::{Multiplicity, StellarSystemMarker, SystemPhenotype};
pub use morphology::{GalaxyMorphology, GalaxyProfile};
pub use preset::{CosmologyConfig, CosmologyPreset};
pub use reference::{GalaxyInstance, StellarSystemRef};
pub use seeds::{galaxy_seed, region_seed, system_seed};
pub use supercluster::Supercluster;
pub use universe::{
    Universe, OBSERVABLE_HUBBLE_KM_S_MPC, OBSERVABLE_OMEGA_MATTER, OBSERVABLE_UNIVERSE_LABEL,
};
pub use region::{
    cylindrical_region_bounds, point_in_galaxy, regions_intersecting_load, RegionBounds, RegionId,
    spherical_region_bounds,
};
pub use solar_system::{BeltSlot, PlanetMassClass, PlanetSlot, SolarSystem};
pub use star::{BarycentricOffset, SpectralClass, Star};
