//! Tectonic setting for context-dependent provinces (e.g. island arc vs Andes).

use serde::{Deserialize, Serialize};

use crate::tectonics::PlateId;

/// Subduction/collision setting for a province.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum TectonicContext {
    /// Oceanic plate subducting under oceanic (island arc).
    OceanicOceanic {
        /// Subducting plate.
        subducting: PlateId,
        /// Overriding plate.
        overriding: PlateId,
    },
    /// Oceanic subducting under continental (Andes-style).
    OceanicContinental {
        /// Subducting plate.
        subducting: PlateId,
        /// Overriding plate.
        overriding: PlateId,
    },
    /// Continent–continent collision.
    ContinentalContinental {
        /// First plate.
        plate_a: PlateId,
        /// Second plate.
        plate_b: PlateId,
    },
}
