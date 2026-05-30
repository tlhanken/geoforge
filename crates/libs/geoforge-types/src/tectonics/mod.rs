//! Tectonic foundation types (Stage 1).

mod boundaries;
mod crust;
mod layers;
mod plates;

pub use boundaries::{BoundarySegmentId, BoundaryStatistics, BoundaryType, PlateBoundaryKind};
pub use crust::CrustType;
pub use layers::{TectonicLayerData, default_crust_age_ma};
pub use plates::{PlateId, PlateMotion, PlateSeed, PlateStats, PlateType};
