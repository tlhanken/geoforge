//! Geoforge v2 surface pipeline. Each layer is a [`geoforge_core::Generator`]; the input
//! is a [`geoforge_core::BodyContext`] and products of earlier layers, never cosmology.

pub mod biomes;
pub mod circulation;
pub mod climate;
pub mod elevation;
pub mod erosion;
pub mod geology;
pub mod hydrology;
pub mod noise;
pub mod preconditions;
pub mod tectonics;

pub use biomes::Biomes;
pub use circulation::Circulation;
pub use climate::Climate;
pub use elevation::Elevation;
pub use geology::Geology;
pub use hydrology::Hydrology;
pub use tectonics::Tectonics;
