//! Geologic province types (Stage 2).

mod context;
mod layers;
mod province;

pub use context::TectonicContext;
pub use layers::{ElevationClass, GeologicLayerData, GeologyRaster, ProvinceInfo};
pub use province::{GeologicProvince, ProvinceCharacteristics};
