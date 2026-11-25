//! Fictional materials a universe may hold: a list the user writes (`exotics.toml` in a
//! universe directory), placed as deposits and markers only. They are flagged fictional
//! wherever they show and change nothing else that is made: no body, orbit, map or number
//! moves when the list changes.
//!
//! ```toml
//! [[material]]
//! name = "Xylarium"
//! note = "a metastable superheavy alloy"
//! colour = "#7fffd4"
//!
//! [[material.found]]
//! in = "asteroid"            # asteroid | world | nebula
//! kinds = ["M"]              # which: asteroid types, world kinds or nebula kinds (empty: any)
//! stars = ["neutron_star"]   # only round these stars (empty: any)
//! chance = 0.002             # per asteroid, world or nebula of those kinds
//! ```

use crate::body::StellarStage;
use crate::Error;
use serde::{Deserialize, Serialize};

/// The list.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct Exotics {
    #[serde(default, rename = "material")]
    pub materials: Vec<Material>,
}

/// A fictional material.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Material {
    pub name: String,
    #[serde(default)]
    pub note: String,
    /// Its markers' colour, `#rrggbb` (default a violet no real map uses).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub colour: Option<String>,
    /// Where it is found.
    #[serde(default)]
    pub found: Vec<Found>,
}

/// One kind of place a material is found.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Found {
    #[serde(rename = "in")]
    pub host: Host,
    /// Which of them: asteroid types ([`ASTEROID_KINDS`]), world kinds ([`WORLD_KINDS`]) or
    /// nebula kinds ([`NEBULA_KINDS`]); empty for any.
    #[serde(default)]
    pub kinds: Vec<String>,
    /// Only in systems with a star of one of these stages; empty for any.
    #[serde(default)]
    pub stars: Vec<StellarStage>,
    /// The chance each candidate holds it, 0–1.
    pub chance: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Host {
    /// Belt members, by spectral type.
    Asteroid,
    /// Planets and moons, by kind.
    World,
    /// Nebulae, by kind.
    Nebula,
}

/// Asteroid types, as [`Found::kinds`] names them.
pub const ASTEROID_KINDS: [&str; 7] = ["S", "C", "M", "V", "P", "D", "icy"];
/// World kinds: a rock surface, an ice surface, open water, a giant, a small body.
pub const WORLD_KINDS: [&str; 5] = ["rocky", "icy", "ocean", "giant", "small"];
/// Nebula kinds.
pub const NEBULA_KINDS: [&str; 5] = ["emission", "reflection", "dark", "planetary", "supernova_remnant"];

/// The colour of a material's markers with none given.
pub const DEFAULT_COLOUR: [u8; 3] = [190, 110, 255];

impl Material {
    /// Its markers' colour.
    #[must_use]
    pub fn rgb(&self) -> [u8; 3] {
        self.colour.as_deref().and_then(parse_colour).unwrap_or(DEFAULT_COLOUR)
    }
}

impl Found {
    /// Whether it names `kind` (any, if it names none).
    #[must_use]
    pub fn admits(&self, kind: &str) -> bool {
        self.kinds.is_empty() || self.kinds.iter().any(|k| k.eq_ignore_ascii_case(kind))
    }

    /// Whether a system with stars of `stages` may hold it.
    #[must_use]
    pub fn admits_stars(&self, stages: &[StellarStage]) -> bool {
        self.stars.is_empty() || stages.iter().any(|s| self.stars.contains(s))
    }
}

impl Exotics {
    /// Checks the list: names given and distinct, chances 0–1, kinds known, colours `#rrggbb`.
    ///
    /// # Errors
    /// Naming the first fault.
    pub fn validate(&self) -> Result<(), Error> {
        let bad = |m: &Material, why: String| Err(Error::Invalid(format!("fictional material '{}': {why}", m.name)));
        let mut names = std::collections::BTreeSet::new();
        for m in &self.materials {
            if m.name.trim().is_empty() {
                return Err(Error::Invalid("a fictional material without a name".into()));
            }
            if !names.insert(m.name.trim().to_lowercase()) {
                return bad(m, "named twice".into());
            }
            if let Some(c) = &m.colour {
                if parse_colour(c).is_none() {
                    return bad(m, format!("colour '{c}' is not #rrggbb"));
                }
            }
            for f in &m.found {
                if !(0.0..=1.0).contains(&f.chance) {
                    return bad(m, format!("chance {} is not 0–1", f.chance));
                }
                let (host, known): (&str, &[&str]) = match f.host {
                    Host::Asteroid => ("asteroid", &ASTEROID_KINDS),
                    Host::World => ("world", &WORLD_KINDS),
                    Host::Nebula => ("nebula", &NEBULA_KINDS),
                };
                if let Some(k) = f.kinds.iter().find(|k| !known.iter().any(|x| x.eq_ignore_ascii_case(k))) {
                    return bad(m, format!("no {host} kind '{k}' (known: {})", known.join(", ")));
                }
            }
        }
        Ok(())
    }

    /// The material named `name` (any case).
    #[must_use]
    pub fn get(&self, name: &str) -> Option<&Material> {
        self.materials.iter().find(|m| m.name.eq_ignore_ascii_case(name.trim()))
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.materials.iter().all(|m| m.found.is_empty())
    }
}

fn parse_colour(c: &str) -> Option<[u8; 3]> {
    let h = c.strip_prefix('#')?;
    if h.len() != 6 {
        return None;
    }
    let byte = |i: usize| u8::from_str_radix(h.get(i..i + 2)?, 16).ok();
    Some([byte(0)?, byte(2)?, byte(4)?])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn checked() {
        let ok = Exotics {
            materials: vec![Material {
                name: "Xylarium".into(),
                note: String::new(),
                colour: Some("#7fffd4".into()),
                found: vec![Found { host: Host::Asteroid, kinds: vec!["m".into()], stars: vec![], chance: 0.01 }],
            }],
        };
        assert!(ok.validate().is_ok());
        assert_eq!(ok.materials[0].rgb(), [0x7f, 0xff, 0xd4]);
        assert!(ok.materials[0].found[0].admits("M") && !ok.materials[0].found[0].admits("S"));
        let mut twice = ok.clone();
        twice.materials.push(Material { name: "xylarium".into(), ..ok.materials[0].clone() });
        assert!(twice.validate().is_err());
        let mut odd = ok.clone();
        odd.materials[0].found[0].kinds = vec!["Q".into()];
        assert!(odd.validate().is_err());
        let mut chance = ok;
        chance.materials[0].found[0].chance = 2.0;
        assert!(chance.validate().is_err());
    }
}
