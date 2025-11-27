//! Geological province generation from tectonic data (Stage 2)
//!
//! This module implements the complete Stage 2 pipeline: converting tectonic plate
//! boundaries and motion into 19 distinct geological province types.
//!
//! # Complete Stage 2 Implementation
//!
//! **Stage 2.1: Orogenic Belts** (3 types)
//! - Collision, Subduction, and Accretionary orogens
//!
//! **Stage 2.2: Large Igneous Provinces** (3 types)
//! - Continental flood basalts, oceanic plateaus, hotspot tracks
//!
//! **Stage 2.3: Arc and Basin Systems** (4 types)
//! - Volcanic arcs, trenches, forearc/backarc basins
//!
//! **Stage 2.4: Stable Continental Regions** (3 types)
//! - Cratons/shields, platforms, extended crust
//!
//! **Stage 2.6: Oceanic Domains** (4 types)
//! - Abyssal plains, mid-ocean ridges, fracture zones, hotspot tracks
//!
//! **Deferred to Future**: Continental Rifts (Stage 2.5)
//! - Will be implemented when needed
//!
//! # Architecture
//!
//! - `provinces`: Province type definitions and characteristics
//! - `orogenic`: Specialized generator for orogenic belts
//! - `generator`: Main generator coordinating all 19 province types
//!
//! # Usage
//!
//! ```rust,no_run
//! use geoforge::map::world::WorldMap;
//!
//! let mut world = WorldMap::new(1800, 900, 42).unwrap();
//! world.tectonics().generate_plates(15).unwrap();
//! world.tectonics().analyze(None).unwrap();
//!
//! // Generate all geological provinces
//! let provinces = world.generate_geology(None).unwrap();
//! println!("Generated {} geological provinces", provinces.len());
//! ```

pub mod provinces;
pub mod orogenic;
pub mod generator;

pub use provinces::{GeologicProvince, ProvinceCharacteristics, ProvinceRegion};
pub use orogenic::{OrogenicBeltGenerator, OrogenicConfig};
pub use generator::{GeologyGenerator, GeologyConfig};
