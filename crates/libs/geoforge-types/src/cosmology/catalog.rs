//! Stage 3 output — galaxy instances and map markers near focus.

use serde::{Deserialize, Serialize};

use super::marker::StellarSystemMarker;
use super::reference::{GalaxyInstance, StellarSystemRef};
use crate::Seed;

/// Galaxy-scale catalog (pipeline stage 3).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GalaxyCatalog {
    /// Root seed (for stage 4 zoom regeneration).
    pub root: Seed,
    /// Generated galaxies in this run.
    pub galaxies: Vec<GalaxyInstance>,
    /// Deterministic primary system reference.
    pub primary_ref: StellarSystemRef,
    /// Markers within the load sphere (capped by config).
    pub markers_near_focus: Vec<StellarSystemMarker>,
}

impl GalaxyCatalog {
    /// Galaxy instance by index.
    #[must_use]
    pub fn galaxy(&self, index: u32) -> Option<&GalaxyInstance> {
        self.galaxies.iter().find(|g| g.index == index)
    }

    /// Marker matching a reference, if present in this catalog.
    #[must_use]
    pub fn marker_for(&self, reference: &StellarSystemRef) -> Option<&StellarSystemMarker> {
        self.markers_near_focus.iter().find(|m| {
            m.galaxy_index == reference.galaxy_index
                && m.region == reference.region
                && m.index_in_region == reference.index_in_region
        })
    }

    /// Primary marker (guaranteed for minimal preset).
    #[must_use]
    pub fn primary_marker(&self) -> Option<&StellarSystemMarker> {
        self.marker_for(&self.primary_ref)
    }
}
