//! Solar-system zoom — stars and orbit slots, no planetary surface interiors.

use serde::{Deserialize, Serialize};

use super::marker::SystemPhenotype;
use super::star::Star;
use crate::coordinates::{GalacticPoint, OrbitalElements};

/// Mass / composition class for a planet slot (not a full body).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum PlanetMassClass {
    /// Rocky terrestrial.
    Rocky,
    /// Super-Earth / mini-Neptune.
    SuperEarth,
    /// Gas giant.
    GasGiant,
    /// Ice giant.
    IceGiant,
    /// Dwarf planet.
    Dwarf,
}

/// Orbital slot — generated at solar-system zoom, not stored on galaxy map.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PlanetSlot {
    /// Stable index within the system.
    pub index: u8,
    /// Orbital elements around primary.
    pub orbit: OrbitalElements,
    /// Mass class.
    pub mass_class: PlanetMassClass,
}

/// Asteroid / ring belt slot.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BeltSlot {
    /// Index within system.
    pub index: u8,
    /// Inner and outer radius (AU).
    pub inner_au: f64,
    pub outer_au: f64,
}

/// Full stellar system at zoom LOD (regenerated from per-system seed).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SolarSystem {
    /// Galaxy index.
    pub galaxy_index: u32,
    /// Region id.
    pub region: super::region::RegionId,
    /// Index within region.
    pub index_in_region: u32,
    /// Absolute galactic position (ly) — matches marker.
    pub barycenter_ly: GalacticPoint,
    /// Map-level phenotype (must match marker).
    pub phenotype: SystemPhenotype,
    /// Component stars.
    pub stars: Vec<Star>,
    /// Planet orbit slots (no atmospheres / surfaces).
    pub planet_slots: Vec<PlanetSlot>,
    /// Belts and rings.
    pub belts: Vec<BeltSlot>,
}
