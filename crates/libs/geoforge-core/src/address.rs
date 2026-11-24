//! Where anything is: the path from a universe down to a galaxy, a star's system and a body.
//!
//! Every item a universe holds is generated from its seed and its address alone, so an
//! address is all it takes to find it again or to generate it for the first time. Written
//! out, the levels are separated by `/`:
//!
//! ```text
//! 7                                  universe 7
//! 7/h6:0,0,-1:17.0                   a galaxy of it
//! 7/h6:0,0,-1:17.0/f12:3,-4,0:5      the system of one of that galaxy's stars
//! 7/h6:0,0,-1:17.0/f12:3,-4,0:5/4    body 4 of that system
//! real/milky-way/sun/4               the real universe: the Milky Way, the Sun, a body
//! ```

use crate::galactic::{GalaxyId, StarId, UniverseSite};
use crate::key::Key;
use crate::seed::Seed;
use serde::{Deserialize, Serialize};
use std::fmt;
use std::str::FromStr;

/// A galaxy of a universe.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum GalaxyRef {
    /// The real Milky Way (in the real universe).
    MilkyWay,
    /// A generated galaxy.
    Generated(GalaxyId),
}

/// A star of a galaxy: the system born with it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum StarRef {
    /// The Sun (in the real Milky Way).
    Sun,
    /// A star of the galaxy's stars.
    Star(StarId),
}

/// How deep an address goes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Level {
    Universe,
    Galaxy,
    System,
    Body,
}

/// Where an item is in a universe. A star needs a galaxy and a body a star; the real
/// universe holds only the Milky Way and the Sun's system among its galaxies and named
/// systems (its other stars are generated as any galaxy's are).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(into = "String", try_from = "String")]
pub struct Address {
    pub universe: UniverseSite,
    pub galaxy: Option<GalaxyRef>,
    pub star: Option<StarRef>,
    pub body: Option<Key>,
}

impl Address {
    /// A universe seeded `seed`.
    #[must_use]
    pub fn universe(seed: Seed) -> Self {
        Self { universe: UniverseSite::Seeded(seed), galaxy: None, star: None, body: None }
    }

    /// The real universe.
    #[must_use]
    pub fn real() -> Self {
        Self { universe: UniverseSite::Real, galaxy: None, star: None, body: None }
    }

    /// The real Solar System.
    #[must_use]
    pub fn sol() -> Self {
        Self { galaxy: Some(GalaxyRef::MilkyWay), star: Some(StarRef::Sun), ..Self::real() }
    }

    /// Galaxy `g` of this address's universe.
    #[must_use]
    pub fn galaxy(self, g: GalaxyRef) -> Self {
        Self { galaxy: Some(g), star: None, body: None, ..self }
    }

    /// The system of star `s` of this address's galaxy.
    #[must_use]
    pub fn star(self, s: StarRef) -> Self {
        Self { star: Some(s), body: None, ..self }
    }

    /// Body `key` of this address's system.
    #[must_use]
    pub fn body(self, key: Key) -> Self {
        Self { body: Some(key), ..self }
    }

    #[must_use]
    pub fn level(&self) -> Level {
        match (self.galaxy, self.star, self.body) {
            (None, ..) => Level::Universe,
            (Some(_), None, _) => Level::Galaxy,
            (Some(_), Some(_), None) => Level::System,
            (Some(_), Some(_), Some(_)) => Level::Body,
        }
    }

    /// The address one level up; `None` for a universe.
    #[must_use]
    pub fn parent(&self) -> Option<Self> {
        match self.level() {
            Level::Universe => None,
            Level::Galaxy => Some(Self { galaxy: None, ..*self }),
            Level::System => Some(Self { star: None, ..*self }),
            Level::Body => Some(Self { body: None, ..*self }),
        }
    }

    /// The address cut to `level` (itself if it is no deeper).
    #[must_use]
    pub fn up_to(&self, level: Level) -> Self {
        let mut a = *self;
        while a.level() > level {
            a = a.parent().expect("deeper than a universe");
        }
        a
    }

    /// Its parts as names safe for a directory each (`:` becomes `_`).
    #[must_use]
    pub fn segments(&self) -> Vec<String> {
        self.to_string().split('/').map(|s| s.replace(':', "_")).collect()
    }

    fn check(self) -> Result<Self, String> {
        let real = self.universe == UniverseSite::Real;
        if self.star.is_some() && self.galaxy.is_none() || self.body.is_some() && self.star.is_none() {
            return Err("a star needs its galaxy and a body its star".into());
        }
        match (self.galaxy, self.star) {
            (Some(GalaxyRef::MilkyWay), _) if !real => Err("the Milky Way is in the real universe only".into()),
            (Some(GalaxyRef::Generated(_)), _) if real => Err("the real universe's galaxies are not generated: real/milky-way".into()),
            (Some(GalaxyRef::Generated(_)), Some(StarRef::Sun)) => Err("the Sun is in the real Milky Way only".into()),
            _ => Ok(self),
        }
    }
}

impl fmt::Display for Address {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.universe {
            UniverseSite::Real => write!(f, "real")?,
            UniverseSite::Seeded(s) => write!(f, "{}", s.0)?,
        }
        match self.galaxy {
            Some(GalaxyRef::MilkyWay) => write!(f, "/milky-way")?,
            Some(GalaxyRef::Generated(id)) => write!(f, "/{id}")?,
            None => {}
        }
        match self.star {
            Some(StarRef::Sun) => write!(f, "/sun")?,
            Some(StarRef::Star(id)) => write!(f, "/{id}")?,
            None => {}
        }
        if let Some(k) = self.body {
            write!(f, "/{k}")?;
        }
        Ok(())
    }
}

impl FromStr for Address {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, String> {
        let mut parts = s.trim().trim_matches('/').split('/');
        let universe = match parts.next().unwrap_or_default() {
            "real" => UniverseSite::Real,
            u => UniverseSite::Seeded(Seed(u.parse().map_err(|_| format!("not a universe: `{u}` (a seed, or `real`)"))?)),
        };
        let galaxy = parts.next().map(|g| if g == "milky-way" { Ok(GalaxyRef::MilkyWay) } else { g.parse().map(GalaxyRef::Generated) }).transpose()?;
        let star = parts.next().map(|t| if t == "sun" { Ok(StarRef::Sun) } else { t.parse().map(StarRef::Star) }).transpose()?;
        let body = parts.next().map(|b| b.parse().map(Key).map_err(|_| format!("not a body: `{b}` (its number)"))).transpose()?;
        if let Some(extra) = parts.next() {
            return Err(format!("nothing is deeper than a body: `{extra}`"));
        }
        Self { universe, galaxy, star, body }.check()
    }
}

impl From<Address> for String {
    fn from(a: Address) -> Self {
        a.to_string()
    }
}

impl TryFrom<String> for Address {
    type Error = String;
    fn try_from(s: String) -> Result<Self, String> {
        s.parse()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn addresses_read_back() {
        for s in ["7", "7/h6:0,0,-1:17.0", "7/h6:0,0,-1:17.0/f12:3,-4,0:5", "7/h6:0,0,-1:17.0/g3:1,2,3:4.5/4", "real", "real/milky-way", "real/milky-way/sun/4", "real/milky-way/a2"] {
            let a: Address = s.parse().unwrap();
            assert_eq!(a.to_string(), s);
            let json = serde_json::to_string(&a).unwrap();
            assert_eq!(serde_json::from_str::<Address>(&json).unwrap(), a);
        }
        let body: Address = "7/h6:0,0,-1:17.0/f12:3,-4,0:5/4".parse().unwrap();
        assert_eq!(body.level(), Level::Body);
        assert_eq!(body.up_to(Level::Galaxy).to_string(), "7/h6:0,0,-1:17.0");
        assert_eq!(body.parent().unwrap().level(), Level::System);
        assert_eq!(body.segments(), ["7", "h6_0,0,-1_17.0", "f12_3,-4,0_5", "4"]);
        assert_eq!(Address::sol().body(Key(4)).to_string(), "real/milky-way/sun/4");
    }

    #[test]
    fn nonsense_is_refused() {
        for s in ["x", "7/milky-way", "real/h6:0,0,0:1.0", "7/h6:0,0,0:1.0/sun", "7/h6:0,0,0:1.0/f1:0,0,0:1/2/3", "7/h6:0,0,0:1.0/f1:0,0,0:1/moon"] {
            assert!(s.parse::<Address>().is_err(), "{s}");
        }
    }
}
