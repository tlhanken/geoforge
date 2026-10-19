//! Stage 5 output — physical parameters per planet slot (no surface rasters yet).

use serde::{Deserialize, Serialize};

use crate::cosmology::PlanetMassClass;
use crate::coordinates::OrbitalElements;
use super::PlanetaryParams;

/// Default equirectangular map width for downstream surface stages.
pub const DEFAULT_MAP_WIDTH: u32 = 360;

/// Default equirectangular map height for downstream surface stages.
pub const DEFAULT_MAP_HEIGHT: u32 = 180;

/// Per-planet summary derived from a solar-system slot (stage 5).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PlanetBodyRecord {
    /// Matching [`PlanetSlot::index`].
    pub slot_index: u8,
    /// Mass / composition class from stage 4.
    pub mass_class: PlanetMassClass,
    /// Orbital elements from stage 4.
    pub orbit: OrbitalElements,
    /// Physical and orbital parameters for simulation stages.
    pub params: PlanetaryParams,
    /// Intended raster width (surface stages not run yet).
    pub map_width: u32,
    /// Intended raster height.
    pub map_height: u32,
}

/// All planetary bodies for one solar system plus generation focus.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PlanetaryBodies {
    /// Slot index used for tectonics / surface pipeline entry.
    pub focus_slot_index: u8,
    /// One record per stage-4 planet slot.
    pub bodies: Vec<PlanetBodyRecord>,
}

impl PlanetaryBodies {
    /// Parameters for the focus planet.
    #[must_use]
    pub fn focus(&self) -> Option<&PlanetBodyRecord> {
        self.bodies
            .iter()
            .find(|b| b.slot_index == self.focus_slot_index)
    }
}
