//! Geological province generation from tectonic data
//!
//! This module implements Stage 2 of the Geoforge pipeline: converting tectonic plate
//! boundaries and motion into realistic geological provinces.
//!
//! # Stage 2.1: Orogenic Belts
//!
//! Mountain-building zones created at convergent plate boundaries:
//! - **Collision Orogens**: Continental-continental convergence (e.g., Himalayas)
//! - **Subduction Orogens**: Oceanic-continental convergence (e.g., Andes)
//! - **Accretionary Orogens**: Mixed plate convergence with terrane accretion
//!
//! # Usage
//!
//! ```rust,no_run
//! use geoforge::geology::orogenic::{OrogenicBeltGenerator, OrogenicConfig};
//! use geoforge::map::world::WorldMap;
//!
//! let mut world = WorldMap::new(1800, 900, 42).unwrap();
//! world.tectonics().generate_plates(15).unwrap();
//! world.analyze_boundaries(None).unwrap();
//!
//! // Generate orogenic belts from convergent boundaries
//! let orogens = world.generate_geology(None).unwrap();
//! println!("Generated {} orogenic belts", orogens.len());
//! ```

pub mod provinces;
pub mod orogenic;

pub use provinces::{GeologicProvince, ProvinceCharacteristics, ProvinceRegion};
pub use orogenic::{OrogenicBeltGenerator, OrogenicConfig};
