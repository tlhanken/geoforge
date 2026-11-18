//! Geoforge v2 surface pipeline. Each layer is a [`geoforge_core::Generator`]; the input
//! is a [`geoforge_core::BodyContext`] and products of earlier layers, never cosmology.

pub mod elevation;
pub mod geology;
pub mod noise;
pub mod tectonics;

pub use elevation::Elevation;
pub use geology::Geology;
pub use tectonics::Tectonics;
