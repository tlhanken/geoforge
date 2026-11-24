//! The orbital tree: what each parent knows about its children.

use super::matter::Composition;
use super::SOLAR_MASS_KG;
use crate::address::Address;
use crate::key::Key;
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
    /// Where the system is in its universe (its address, to the system's level); `None` for
    /// systems made before addresses and hand-made ones. Never affects seeds.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub placement: Option<Address>,
    /// The stellar population the system was drawn from, when it was drawn from a galaxy.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub site: Option<GalacticSite>,
    /// The present as a Julian date (Terrestrial Time), for a system with a real calendar
    /// (the Solar System).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub present_jd_tt: Option<f64>,
}

/// Where a system is in its galaxy: the roster's "location".
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GalacticSite {
    /// The galaxy model's name, e.g. "Milky Way".
    pub galaxy: String,
    /// Distance from the galaxy's centre, kpc (the Sun: 8.2).
    pub radius_kpc: f64,
    /// Position in the galaxy, pc from its centre, z toward its north pole. `None` (systems
    /// made before positions): at `radius_kpc` in the plane, at an azimuth from the system's
    /// seed.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub position_pc: Option<[f64; 3]>,
    /// The seed of the galaxy's stars. `None`: drawn from the system's seed.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub galaxy_seed: Option<crate::seed::Seed>,
    /// The star the system was made at, if any.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub star: Option<crate::galactic::StarId>,
    /// How the system's frame sits in the galaxy: its rows are the system frame's axes in
    /// galactic coordinates as seen from the system (x toward the galaxy's centre, z toward
    /// its north pole). `None`: a rotation from the system's seed.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub orientation: Option<[[f64; 3]; 3]>,
    /// The universe around the galaxy. `None`: one seeded from the galaxy's seed.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub universe: Option<crate::galactic::UniverseSite>,
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

impl Elements {
    /// These elements `t_days` later under `p`: the node and periapsis turned, the mean
    /// anomaly held back by the periapsis's advance so the mean longitude keeps its pace.
    #[must_use]
    pub fn precessed(&self, p: Precession, t_days: f64) -> Self {
        let (node, peri) = (p.node_deg_day * t_days, p.periapsis_deg_day * t_days);
        Self { lan_deg: self.lan_deg + node, argp_deg: self.argp_deg + peri - node, m0_deg: self.m0_deg - peri, ..*self }
    }

    /// Position relative to the parent at mean anomaly `mean_deg`, in `a_au`'s units, in the
    /// plane the elements are measured from (x toward the reference direction).
    #[must_use]
    pub fn position_at(&self, mean_deg: f64) -> [f64; 3] {
        let mean = mean_deg.to_radians().rem_euclid(std::f64::consts::TAU);
        let ecc = self.e.clamp(0.0, 0.99);
        // Kepler's equation, E − e·sin E = M, by Newton's method.
        let mut big_e = if ecc > 0.8 { std::f64::consts::PI } else { mean };
        for _ in 0..50 {
            let step = (big_e - ecc * big_e.sin() - mean) / (1.0 - ecc * big_e.cos());
            big_e -= step;
            if step.abs() < 1e-13 {
                break;
            }
        }
        let nu = 2.0 * ((1.0 + ecc).sqrt() * (big_e / 2.0).sin()).atan2((1.0 - ecc).sqrt() * (big_e / 2.0).cos());
        let r = self.a_au * (1.0 - ecc * big_e.cos());
        let u = self.argp_deg.to_radians() + nu;
        let (node, inc) = (self.lan_deg.to_radians(), self.i_deg.to_radians());
        [r * (node.cos() * u.cos() - node.sin() * u.sin() * inc.cos()), r * (node.sin() * u.cos() + node.cos() * u.sin() * inc.cos()), r * u.sin() * inc.sin()]
    }
}

/// A star as it is now. A remnant (white dwarf, neutron star, black hole) is a star here too:
/// it is what the system's bodies orbit.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Star {
    pub mass_solar: f64,
    pub luminosity_solar: f64,
    pub radius_solar: f64,
    pub temperature_k: f64,
    pub age_gyr: f64,
    pub metallicity_dex: f64,
    pub spectral_class: String,
    #[serde(default, skip_serializing_if = "StellarStage::is_main_sequence")]
    pub stage: StellarStage,
    /// Mass on the zero-age main sequence, M☉. `None`: the star has lost no mass to speak of
    /// and this is `mass_solar`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub initial_mass_solar: Option<f64>,
}

impl Star {
    /// Mass on the zero-age main sequence, M☉.
    #[must_use]
    pub fn initial_mass(&self) -> f64 {
        self.initial_mass_solar.unwrap_or(self.mass_solar)
    }
}

/// Where a star is in its life.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum StellarStage {
    /// Burning hydrogen in its core: a dwarf (luminosity class V).
    #[default]
    MainSequence,
    /// Core hydrogen spent, swelling toward a giant (IV).
    Subgiant,
    /// A red giant: burning hydrogen in a shell, then helium (III).
    Giant,
    /// A massive star after the main sequence, blue to red (I).
    Supergiant,
    /// The carbon–oxygen core left by a star under ~8 M☉, cooling.
    WhiteDwarf,
    /// Left by a core-collapse supernova of a star of ~8–20 M☉.
    NeutronStar,
    /// Left by the collapse of a star over ~20 M☉.
    BlackHole,
    /// Too light to burn hydrogen (under ~0.08 M☉): glowing from its birth heat, cooling.
    BrownDwarf,
}

impl StellarStage {
    #[must_use]
    pub fn is_main_sequence(&self) -> bool {
        *self == Self::MainSequence
    }
    /// A white dwarf, neutron star or black hole.
    #[must_use]
    pub fn is_remnant(&self) -> bool {
        matches!(self, Self::WhiteDwarf | Self::NeutronStar | Self::BlackHole)
    }
    #[must_use]
    pub fn label(&self) -> &'static str {
        match self {
            Self::MainSequence => "main-sequence star",
            Self::Subgiant => "subgiant",
            Self::Giant => "red giant",
            Self::Supergiant => "supergiant",
            Self::WhiteDwarf => "white dwarf",
            Self::NeutronStar => "neutron star",
            Self::BlackHole => "black hole",
            Self::BrownDwarf => "brown dwarf",
        }
    }
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
    /// How the orbit turns, when it does enough to matter (the Moon's).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub precession: Option<Precession>,
}

/// An orbit turning at steady rates, degrees per day: its node (negative: regressing, as the
/// Moon's does once in 18.6 years under the Sun's pull) and its periapsis longitude. The
/// mean longitude keeps the sidereal period.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Precession {
    pub node_deg_day: f64,
    pub periapsis_deg_day: f64,
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
