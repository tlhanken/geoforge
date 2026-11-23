//! Geoforge surface layers. Each layer is a [`geoforge_core::Generator`]; the input
//! is a [`geoforge_core::BodyContext`] and products of earlier layers, never cosmology.

pub mod atmosphere;
pub mod bands;
pub mod biomes;
pub mod circulation;
pub mod climate;
pub mod elevation;
pub mod erosion;
pub mod geology;
pub mod hydrology;
pub mod lid;
pub mod noise;
pub mod preconditions;
pub mod relief;
pub mod shell;
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
pub use shell::{IceShellElevation, IceShellGeology, IceShellTectonics};
pub use tectonics::Tectonics;
