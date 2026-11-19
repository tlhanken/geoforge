//! The handoff from the solar system to the surface pipeline.

use super::interior::BodyInterior;
use super::system::{BodyExterior, Elements, Epoch};
use crate::key::Key;
use serde::{Deserialize, Serialize};

/// Starlight reaching the body and how it changes through the year.
///
/// For a moon the year is its parent's orbit around the star, so `eccentricity` is the
/// parent's and `obliquity_deg` is the moon's spin axis against that orbit's plane.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Irradiance {
    /// Mean flux over the year, W/m².
    pub flux_w_m2: f64,
    pub star_temperature_k: f64,
    /// Eccentricity of the path around the star.
    pub eccentricity: f64,
    /// Tilt of the spin axis against the path around the star, degrees.
    pub obliquity_deg: f64,
}

/// A body raising tides on this one.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TidalPartner {
    pub name: String,
    pub mass_kg: f64,
    pub distance_km: f64,
}

/// Everything the surface pipeline needs, and nothing else. Generated from a system, or
/// loaded from a hand-written `body.toml` (start from an [`Archetype`](super::Archetype)).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BodyContext {
    pub key: Key,
    pub name: String,
    pub exterior: BodyExterior,
    pub interior: BodyInterior,
    /// The body's own orbit (around its parent, for a moon).
    pub orbit: Elements,
    pub irradiance: Irradiance,
    #[serde(default)]
    pub tides: Vec<TidalPartner>,
    pub epoch: Epoch,
}

impl BodyContext {
    /// Surface gravity, m/s².
    #[must_use]
    pub fn gravity_m_s2(&self) -> f64 {
        crate::derived::surface_gravity_m_s2(self.exterior.mass_kg, self.exterior.radius_km)
    }
    #[must_use]
    pub fn surface_area_km2(&self) -> f64 {
        4.0 * std::f64::consts::PI * self.exterior.radius_km.powi(2)
    }
}
