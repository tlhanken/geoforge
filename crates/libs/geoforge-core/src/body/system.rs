//! The orbital tree: what each parent knows about its children.

use super::matter::Composition;
use super::SOLAR_MASS_KG;
use crate::key::{Address, Key};
use crate::seed::Seed;
use serde::{Deserialize, Serialize};

/// The world's present moment.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Epoch {
    /// Age of the system at the present, billions of years.
    pub system_age_gyr: f64,
}

/// Root of a thin world: one solar system with its own key and seed.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SystemRoot {
    pub key: Key,
    pub seed: Seed,
    pub name: String,
    pub epoch: Epoch,
    /// Filled when the system is attached to a galaxy sector. Never affects seeds.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub placement: Option<Address>,
}

/// Keplerian orbital elements. The period is derived from `a_au` and the parent's mass.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Elements {
    pub a_au: f64,
    pub e: f64,
    pub i_deg: f64,
    pub lan_deg: f64,
    pub argp_deg: f64,
    pub m0_deg: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Star {
    pub mass_solar: f64,
    pub luminosity_solar: f64,
    pub radius_solar: f64,
    pub temperature_k: f64,
    pub age_gyr: f64,
    pub metallicity_dex: f64,
    pub spectral_class: String,
}

/// Physical class of a body. "Moon" is a role (its parent is a body), not a class.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum BodyClass {
    /// Solid and massive enough to be round: rocky planets and large moons, rocky or icy.
    Terrestrial,
    GasGiant,
    IceGiant,
    /// Not round: asteroids, comets, small moons.
    SmallBody,
}

/// What a parent knows about a body.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct BodyExterior {
    pub class: BodyClass,
    pub mass_kg: f64,
    pub radius_km: f64,
    pub rotation_h: f64,
    pub obliquity_deg: f64,
}

/// Size distribution of a belt or ring population.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PowerLaw {
    pub min_km: f64,
    pub max_km: f64,
    /// Differential exponent q in dN/dD ∝ D^-q.
    pub exponent: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Population {
    pub inner_au: f64,
    pub outer_au: f64,
    pub count_estimate: f64,
    pub size_law: PowerLaw,
    pub composition: Composition,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum NodeKind {
    Barycenter,
    Star(Star),
    Body(BodyExterior),
    Belt(Population),
    Ring(Population),
}

/// A node of a system's orbital tree.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct OrbitalNode {
    pub key: Key,
    pub name: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub parent: Option<Key>,
    pub node: NodeKind,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub elements: Option<Elements>,
}

/// `system.json`: the root plus the orbital tree (exteriors only).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SolarSystem {
    pub root: SystemRoot,
    pub nodes: Vec<OrbitalNode>,
}

impl SolarSystem {
    #[must_use]
    pub fn node(&self, key: Key) -> Option<&OrbitalNode> {
        self.nodes.iter().find(|n| n.key == key)
    }
    #[must_use]
    pub fn children(&self, key: Key) -> Vec<&OrbitalNode> {
        self.nodes.iter().filter(|n| n.parent == Some(key)).collect()
    }
    #[must_use]
    pub fn by_name(&self, name: &str) -> Option<&OrbitalNode> {
        self.nodes.iter().find(|n| n.name.eq_ignore_ascii_case(name))
    }
    /// Total mass orbited by `key` (its parent's mass, or the barycentre's members' masses).
    #[must_use]
    pub fn central_mass_kg(&self, key: Key) -> f64 {
        let Some(parent) = self.node(key).and_then(|n| n.parent) else { return 0.0 };
        self.mass_of(parent)
    }
    #[must_use]
    pub fn mass_of(&self, key: Key) -> f64 {
        match self.node(key).map(|n| &n.node) {
            Some(NodeKind::Star(s)) => s.mass_solar * SOLAR_MASS_KG,
            Some(NodeKind::Body(b)) => b.mass_kg,
            Some(NodeKind::Barycenter) => self.children(key).iter().map(|c| self.mass_of(c.key)).sum(),
            _ => 0.0,
        }
    }
}
