//! Core domain types for Geoforge v2.
//!
//! This crate defines seeds, coordinates, planetary parameters, and typed descriptions
//! of each pipeline stage's outputs. It intentionally contains **no** raster storage or
//! simulation logic — those belong in `geoforge-grid` and stage crates.
//!
//! # Quick example
//!
//! ```
//! use geoforge_types::{PlanetaryParams, Seed, PipelineLayerId};
//! use geoforge_types::tectonics::CrustType;
//!
//! let root = Seed::new(42);
//! let tectonics_seed = root.child(b"tectonics");
//! let params = PlanetaryParams::earth();
//! assert_eq!(params.radius_km, 6371.0);
//! assert!(PipelineLayerId::Tectonics < PipelineLayerId::Elevation);
//! assert_ne!(CrustType::Oceanic, CrustType::Continental);
//! ```

#![deny(clippy::todo, clippy::unimplemented)]
#![warn(missing_docs)]

pub mod coordinates;
pub mod cosmology;
pub mod error;
pub mod geology;
pub mod pipeline;
pub mod planetary;
pub mod seed;
pub mod tectonics;

pub use cosmology::{
    CosmologyConfig, CosmologyPreset, GalaxyInstance, GalaxyMorphology, GalaxyProfile,
    StellarSystemMarker, StellarSystemRef, SolarSystem,
};
pub use error::TypesError;
pub use pipeline::PipelineLayerId;
pub use planetary::PlanetaryParams;
pub use seed::Seed;

#[cfg(test)]
mod lib_tests;
