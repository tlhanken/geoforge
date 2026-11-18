//! Bodies, orbits and solar systems.
//!
//! A parent owns each child's *exterior* (what is visible from outside: mass, radius, orbit,
//! luminosity). The child owns its *interior* (the shell stack). The surface pipeline sees
//! only a [`BodyContext`].

use crate::grid::GridSpec;
use crate::key::{Address, Key};
use crate::seed::Seed;
use serde::{Deserialize, Serialize};

pub const EARTH_MASS_KG: f64 = 5.972e24;
pub const EARTH_RADIUS_KM: f64 = 6371.0;
pub const SOLAR_MASS_KG: f64 = 1.989e30;
pub const SOLAR_CONSTANT_W_M2: f64 = 1361.0;
pub const AU_KM: f64 = 1.495_978_707e8;

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
    Terrestrial,
    GasGiant,
    IceGiant,
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Species {
    Iron,
    Nickel,
    Sulfur,
    Silicate,
    Basalt,
    Granite,
    Water,
    Ice,
    Hydrogen,
    Helium,
    Nitrogen,
    Oxygen,
    CarbonDioxide,
    Argon,
    Methane,
    Ammonia,
}

/// Species fractions by mass, sorted, summing to 1.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(try_from = "Vec<(Species, f64)>", into = "Vec<(Species, f64)>")]
pub struct Composition(Vec<(Species, f64)>);

impl Composition {
    /// # Errors
    /// Fractions must be non-negative, species unique, and sum to 1 within 1e-6.
    pub fn new(mut parts: Vec<(Species, f64)>) -> Result<Self, crate::Error> {
        parts.sort_by_key(|a| a.0);
        if parts.windows(2).any(|w| w[0].0 == w[1].0) {
            return Err(crate::Error::Invalid("composition lists a species twice".into()));
        }
        if parts.iter().any(|p| p.1 < 0.0 || !p.1.is_finite()) {
            return Err(crate::Error::Invalid("composition fractions must be >= 0".into()));
        }
        let sum: f64 = parts.iter().map(|p| p.1).sum();
        if (sum - 1.0).abs() > 1e-6 {
            return Err(crate::Error::Invalid(format!("composition sums to {sum}, not 1")));
        }
        Ok(Self(parts))
    }
    /// Normalize arbitrary positive weights.
    #[must_use]
    pub fn from_weights(parts: &[(Species, f64)]) -> Self {
        let sum: f64 = parts.iter().map(|p| p.1).sum();
        let mut v: Vec<_> = parts.iter().map(|&(s, w)| (s, w / sum)).collect();
        v.sort_by_key(|a| a.0);
        Self(v)
    }
    #[must_use]
    pub fn single(s: Species) -> Self {
        Self(vec![(s, 1.0)])
    }
    #[must_use]
    pub fn parts(&self) -> &[(Species, f64)] {
        &self.0
    }
    #[must_use]
    pub fn fraction(&self, s: Species) -> f64 {
        self.0.iter().find(|p| p.0 == s).map_or(0.0, |p| p.1)
    }
}

impl TryFrom<Vec<(Species, f64)>> for Composition {
    type Error = crate::Error;
    fn try_from(v: Vec<(Species, f64)>) -> Result<Self, Self::Error> {
        Self::new(v)
    }
}
impl From<Composition> for Vec<(Species, f64)> {
    fn from(c: Composition) -> Self {
        c.0
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ShellKind {
    Core,
    Mantle,
    Crust,
    Ocean,
    IceShell,
    Atmosphere,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Phase {
    Solid,
    Liquid,
    Gas,
    Supercritical,
}

/// One layer of a body, innermost first in [`BodyInterior::shells`].
///
/// For an ocean, `thickness_km` is the global-equivalent layer: the depth if the water were
/// spread evenly over the whole body. For an atmosphere, `pressure_bar` is surface pressure.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Shell {
    pub kind: ShellKind,
    pub phase: Phase,
    pub thickness_km: f64,
    pub composition: Composition,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pressure_bar: Option<f64>,
}

/// The body's own layer: an ordered shell stack, core outward.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BodyInterior {
    pub shells: Vec<Shell>,
}

impl BodyInterior {
    #[must_use]
    pub fn shell(&self, kind: ShellKind) -> Option<&Shell> {
        self.shells.iter().rev().find(|s| s.kind == kind)
    }
    /// The outermost solid layer people could stand on (rock crust or ice shell), if any.
    #[must_use]
    pub fn surface_shell(&self) -> Option<&Shell> {
        self.shells.iter().rev().find(|s| matches!(s.kind, ShellKind::Crust | ShellKind::IceShell))
    }
    #[must_use]
    pub fn has_solid_surface(&self) -> bool {
        // A solid shell must not lie under a deep gas envelope (gas giants).
        let Some(idx) = self.shells.iter().rposition(|s| matches!(s.kind, ShellKind::Crust | ShellKind::IceShell)) else {
            return false;
        };
        self.shells[idx + 1..].iter().all(|s| match s.kind {
            ShellKind::Atmosphere => s.pressure_bar.unwrap_or(0.0) < 1000.0,
            ShellKind::Ocean => true,
            _ => false,
        })
    }
    /// Rock crust that can host plate tectonics.
    #[must_use]
    pub fn has_lithosphere(&self) -> bool {
        self.has_solid_surface() && self.surface_shell().is_some_and(|s| s.kind == ShellKind::Crust)
    }
    #[must_use]
    pub fn has_atmosphere(&self) -> bool {
        self.shell(ShellKind::Atmosphere).is_some_and(|s| s.pressure_bar.unwrap_or(0.0) > 1e-3)
    }
    /// Surface ocean (above the crust) global-equivalent depth, km.
    #[must_use]
    pub fn surface_ocean_gel_km(&self) -> f64 {
        let crust = self.shells.iter().rposition(|s| s.kind == ShellKind::Crust);
        self.shells
            .iter()
            .enumerate()
            .filter(|(i, s)| s.kind == ShellKind::Ocean && crust.is_some_and(|c| *i > c))
            .map(|(_, s)| s.thickness_km)
            .sum()
    }
    #[must_use]
    pub fn surface_pressure_bar(&self) -> f64 {
        self.shell(ShellKind::Atmosphere).and_then(|s| s.pressure_bar).unwrap_or(0.0)
    }
}

/// Starlight reaching the body.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Irradiance {
    /// Mean flux at the body's semi-major axis, W/m².
    pub flux_w_m2: f64,
    pub star_temperature_k: f64,
}

/// A body raising tides on this one.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TidalPartner {
    pub name: String,
    pub mass_kg: f64,
    pub distance_km: f64,
}

/// Everything the surface pipeline needs, and nothing else. Generated from a system, or
/// loaded from a hand-written `body.toml`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BodyContext {
    pub key: Key,
    pub name: String,
    pub exterior: BodyExterior,
    pub interior: BodyInterior,
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
    /// An Earth-like default body, handy for tests and authoring templates.
    #[must_use]
    pub fn earth_like(key: Key, name: &str) -> Self {
        Self {
            key,
            name: name.into(),
            exterior: BodyExterior {
                class: BodyClass::Terrestrial,
                mass_kg: EARTH_MASS_KG,
                radius_km: EARTH_RADIUS_KM,
                rotation_h: 23.934,
                obliquity_deg: 23.44,
            },
            interior: BodyInterior::earth_like(),
            orbit: Elements { a_au: 1.0, e: 0.0167, i_deg: 0.0, lan_deg: 0.0, argp_deg: 102.9, m0_deg: 0.0 },
            irradiance: Irradiance { flux_w_m2: SOLAR_CONSTANT_W_M2, star_temperature_k: 5772.0 },
            tides: vec![TidalPartner { name: "Moon".into(), mass_kg: 7.342e22, distance_km: 384_400.0 }],
            epoch: Epoch { system_age_gyr: 4.57 },
        }
    }
}

impl BodyInterior {
    #[must_use]
    pub fn earth_like() -> Self {
        use Species::*;
        Self {
            shells: vec![
                Shell { kind: ShellKind::Core, phase: Phase::Solid, thickness_km: 3480.0, composition: Composition::from_weights(&[(Iron, 0.85), (Nickel, 0.05), (Sulfur, 0.10)]), pressure_bar: None },
                Shell { kind: ShellKind::Mantle, phase: Phase::Solid, thickness_km: 2856.0, composition: Composition::single(Silicate), pressure_bar: None },
                Shell { kind: ShellKind::Crust, phase: Phase::Solid, thickness_km: 35.0, composition: Composition::from_weights(&[(Basalt, 0.6), (Granite, 0.4)]), pressure_bar: None },
                Shell { kind: ShellKind::Ocean, phase: Phase::Liquid, thickness_km: 2.62, composition: Composition::single(Water), pressure_bar: None },
                Shell { kind: ShellKind::Atmosphere, phase: Phase::Gas, thickness_km: 100.0, composition: Composition::from_weights(&[(Nitrogen, 0.755), (Oxygen, 0.232), (Argon, 0.013)]), pressure_bar: Some(1.013) },
            ],
        }
    }

    #[must_use]
    pub fn gas_giant() -> Self {
        use Species::*;
        Self {
            shells: vec![
                Shell { kind: ShellKind::Core, phase: Phase::Solid, thickness_km: 10_000.0, composition: Composition::from_weights(&[(Silicate, 0.5), (Iron, 0.5)]), pressure_bar: None },
                Shell { kind: ShellKind::Mantle, phase: Phase::Liquid, thickness_km: 40_000.0, composition: Composition::from_weights(&[(Hydrogen, 0.75), (Helium, 0.25)]), pressure_bar: None },
                Shell { kind: ShellKind::Atmosphere, phase: Phase::Gas, thickness_km: 20_000.0, composition: Composition::from_weights(&[(Hydrogen, 0.86), (Helium, 0.14)]), pressure_bar: Some(1.0e6) },
            ],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn composition_checks() {
        assert!(Composition::new(vec![(Species::Water, 0.5), (Species::Ice, 0.5)]).is_ok());
        assert!(Composition::new(vec![(Species::Water, 0.5)]).is_err());
        assert!(Composition::new(vec![(Species::Water, 0.5), (Species::Water, 0.5)]).is_err());
    }

    #[test]
    fn presence_from_shells() {
        let earth = BodyInterior::earth_like();
        assert!(earth.has_solid_surface() && earth.has_lithosphere() && earth.has_atmosphere());
        assert!((earth.surface_ocean_gel_km() - 2.62).abs() < 1e-9);
        let giant = BodyInterior::gas_giant();
        assert!(!giant.has_solid_surface());
        assert!(!giant.has_lithosphere());
        // Europa: ice shell over an ocean over rock — a solid surface, but no lithosphere.
        let europa = BodyInterior {
            shells: vec![
                Shell { kind: ShellKind::Mantle, phase: Phase::Solid, thickness_km: 1400.0, composition: Composition::single(Species::Silicate), pressure_bar: None },
                Shell { kind: ShellKind::Ocean, phase: Phase::Liquid, thickness_km: 100.0, composition: Composition::single(Species::Water), pressure_bar: None },
                Shell { kind: ShellKind::IceShell, phase: Phase::Solid, thickness_km: 20.0, composition: Composition::single(Species::Ice), pressure_bar: None },
            ],
        };
        assert!(europa.has_solid_surface());
        assert!(!europa.has_lithosphere());
        assert_eq!(europa.surface_ocean_gel_km(), 0.0, "subsurface ocean is not a surface ocean");
    }
}

/// Grid used for a body's surface products, recorded in the manifest.
pub type SurfaceGrid = GridSpec;
