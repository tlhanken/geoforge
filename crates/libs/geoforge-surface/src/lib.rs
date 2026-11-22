//! Geoforge v2 surface pipeline. Each layer is a [`geoforge_core::Generator`]; the input
//! is a [`geoforge_core::BodyContext`] and products of earlier layers, never cosmology.

pub mod atmosphere;
pub mod bands;
pub mod biomes;
pub mod circulation;
pub mod climate;
pub mod elevation;
pub mod erosion;
pub mod geology;
pub mod ground;
pub mod hazards;
pub mod hydrology;
pub mod koppen;
pub mod lid;
pub mod noise;
pub mod preconditions;
pub mod relief;
pub mod resources;
pub mod tectonics;

pub use bands::GiantBands;
pub use biomes::Biomes;
pub use circulation::Circulation;
pub use climate::Climate;
pub use elevation::Elevation;
pub use erosion::Erosion;
pub use geology::Geology;
pub use hydrology::Hydrology;
pub use lid::{LidElevation, LidGeology, LidTectonics};
pub use tectonics::Tectonics;
