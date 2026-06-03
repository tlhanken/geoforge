//! Seed-driven cosmology: region streaming, galaxy markers, solar-system zoom.
//!
//! All output is regenerable from [`geoforge_types::Seed`] + [`CosmologyConfig`].
//! No JSON persistence is required.

#![deny(clippy::todo, clippy::unimplemented)]

mod context;
mod density;
mod galaxy;
mod marker;
mod phenotype;
mod rng;
mod seeds;
mod solar_system;

pub use context::CosmologyContext;
pub use galaxy::{galaxy_count, galaxy_instance};
pub use marker::markers_in_region;
pub use solar_system::generate_solar_system;

#[cfg(test)]
mod lib_tests;
