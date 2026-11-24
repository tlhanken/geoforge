//! Stars in a galaxy as data: their identities, what they are, and the edits a world keeps
//! over the generated ones. Generation lives in `geoforge-cosmos::sectors`.

use crate::seed::Seed;
use serde::{Deserialize, Serialize};
use std::fmt;
use std::str::FromStr;

/// A galaxy's stable identity in its universe: its halo's mass bin, the cube of that bin's
/// grid it was generated in and its number there, and its number among the halo's galaxies
/// (0 the central, then the satellites).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(into = "String", try_from = "String")]
pub struct GalaxyId {
    pub bin: u8,
    pub cube: [i32; 3],
    pub n: u32,
    pub member: u32,
}

impl GalaxyId {
    /// The seed of the galaxy in a universe seeded `universe`.
    #[must_use]
    pub fn seed(&self, universe: Seed) -> Seed {
        let s = self.cube.iter().fold(universe.child_u64(b"halo bin", u64::from(self.bin)), |s, &v| s.child_u64(b"axis", u64::from(v as u32)));
        s.child_u64(b"n", u64::from(self.n)).child_u64(b"member", u64::from(self.member))
    }
}

impl fmt::Display for GalaxyId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let [i, j, k] = self.cube;
        write!(f, "h{}:{i},{j},{k}:{}.{}", self.bin, self.n, self.member)
    }
}

impl FromStr for GalaxyId {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, String> {
        let bad = || format!("not a galaxy identity: {s} (h<bin>:<i>,<j>,<k>:<n>.<member>)");
        let rest = s.strip_prefix('h').ok_or_else(bad)?;
        let parts: Vec<&str> = rest.split(':').collect();
        let [bin, cube, n] = parts.as_slice() else { return Err(bad()) };
        let cube: Vec<i32> = cube.split(',').map(str::parse).collect::<Result<_, _>>().map_err(|_| bad())?;
        let (n, member) = n.split_once('.').ok_or_else(bad)?;
        Ok(Self { bin: bin.parse().map_err(|_| bad())?, cube: cube.try_into().map_err(|_| bad())?, n: n.parse().map_err(|_| bad())?, member: member.parse().map_err(|_| bad())? })
    }
}

impl From<GalaxyId> for String {
    fn from(id: GalaxyId) -> Self {
        id.to_string()
    }
}

impl TryFrom<String> for GalaxyId {
    type Error = String;
    fn try_from(s: String) -> Result<Self, String> {
        s.parse()
    }
}

/// The universe a system's galaxy sits in.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum UniverseSite {
    /// The real Local Group and nearby galaxies, as catalogued (the Solar System).
    Real,
    /// A universe made from this seed; the home galaxy is the nearest one to its origin like
    /// the system's galaxy.
    Seeded(Seed),
}

/// A star's stable identity: where it was generated (or that it was added by hand). It is
/// the same for every observer and seeds the system made at the star.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(into = "String", try_from = "String")]
pub enum StarId {
    /// The `n`th field star of brightness bin `bin` in cube `cube` of the bin's level.
    Field { bin: u8, cube: [i32; 3], n: u32 },
    /// The `n`th member in bin `bin` of birth event `event` of group cube `cube`.
    Group { bin: u8, cube: [i32; 3], event: u32, n: u32 },
    /// The `n`th star added by hand.
    Added { n: u32 },
    /// The `n`th member of hand-added cluster `cluster`.
    Cluster { cluster: u32, n: u32 },
}

impl StarId {
    /// The seed of the star's system in a galaxy seeded `galaxy`.
    #[must_use]
    pub fn seed(&self, galaxy: Seed) -> Seed {
        let s = galaxy.child(b"star");
        let cube = |s: Seed, c: [i32; 3]| c.iter().fold(s, |s, &v| s.child_u64(b"axis", u64::from(v as u32)));
        match *self {
            Self::Field { bin, cube: c, n } => cube(s.child_u64(b"field", u64::from(bin)), c).child_u64(b"n", u64::from(n)),
            Self::Group { bin, cube: c, event, n } => cube(s.child_u64(b"group", u64::from(bin)), c).child_u64(b"event", u64::from(event)).child_u64(b"n", u64::from(n)),
            Self::Added { n } => s.child_u64(b"added", u64::from(n)),
            Self::Cluster { cluster, n } => s.child_u64(b"cluster", u64::from(cluster)).child_u64(b"n", u64::from(n)),
        }
    }
}

impl fmt::Display for StarId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Field { bin, cube: [i, j, k], n } => write!(f, "f{bin}:{i},{j},{k}:{n}"),
            Self::Group { bin, cube: [i, j, k], event, n } => write!(f, "g{bin}:{i},{j},{k}:{event}.{n}"),
            Self::Added { n } => write!(f, "a{n}"),
            Self::Cluster { cluster, n } => write!(f, "c{cluster}:{n}"),
        }
    }
}

impl FromStr for StarId {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, String> {
        let bad = || format!("not a star identity: {s} (f<bin>:<i>,<j>,<k>:<n>, g<bin>:<i>,<j>,<k>:<event>.<n>, a<n> or c<cluster>:<n>)");
        let num = |t: &str| t.parse::<u32>().map_err(|_| bad());
        let cube = |t: &str| -> Result<[i32; 3], String> {
            let v: Vec<i32> = t.split(',').map(str::parse).collect::<Result<_, _>>().map_err(|_| bad())?;
            v.try_into().map_err(|_| bad())
        };
        let (tag, rest) = s.split_at_checked(1).ok_or_else(bad)?;
        let parts: Vec<&str> = rest.split(':').collect();
        match (tag, parts.as_slice()) {
            ("f", [bin, c, n]) => Ok(Self::Field { bin: bin.parse().map_err(|_| bad())?, cube: cube(c)?, n: num(n)? }),
            ("g", [bin, c, en]) => {
                let (e, n) = en.split_once('.').ok_or_else(bad)?;
                Ok(Self::Group { bin: bin.parse().map_err(|_| bad())?, cube: cube(c)?, event: num(e)?, n: num(n)? })
            }
            ("a", [n]) => Ok(Self::Added { n: num(n)? }),
            ("c", [c, n]) => Ok(Self::Cluster { cluster: num(c)?, n: num(n)? }),
            _ => Err(bad()),
        }
    }
}

impl From<StarId> for String {
    fn from(id: StarId) -> Self {
        id.to_string()
    }
}

impl TryFrom<String> for StarId {
    type Error = String;
    fn try_from(s: String) -> Result<Self, String> {
        s.parse()
    }
}

/// What a star is, as it was born: enough to evolve it to any age.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct StarKind {
    /// Initial mass, M☉.
    pub initial_mass: f64,
    pub age_gyr: f64,
    pub metallicity_dex: f64,
    /// An unresolved companion's initial mass, M☉.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub companion_mass: Option<f64>,
}

/// A star as an edit holds it: where it is in its galaxy (pc, from the centre, z toward the
/// north pole) and what it is.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StarRecord {
    pub id: StarId,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub name: Option<String>,
    pub position_pc: [f64; 3],
    pub kind: StarKind,
}

/// A cluster added by hand: its members are generated from its own seed.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ClusterEdit {
    /// Its number, which its members' identities carry.
    pub n: u32,
    pub name: String,
    pub centre_pc: [f64; 3],
    /// Plummer radius, pc.
    pub radius_pc: f64,
    /// Stars, of every brightness.
    pub members: u32,
    pub age_gyr: f64,
    pub metallicity_dex: f64,
}

/// A world's edits over its galaxy's generated stars, applied in every query.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct StarEdits {
    /// Generated stars taken away.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub remove: Vec<StarId>,
    /// Generated stars changed (moved, renamed or made something else), whole.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub change: Vec<StarRecord>,
    /// Stars added by hand (`StarId::Added`).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub add: Vec<StarRecord>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub clusters: Vec<ClusterEdit>,
}

impl StarEdits {
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.remove.is_empty() && self.change.is_empty() && self.add.is_empty() && self.clusters.is_empty()
    }

    /// The identity the next added star takes.
    #[must_use]
    pub fn next_added(&self) -> StarId {
        let n = self.add.iter().filter_map(|r| if let StarId::Added { n } = r.id { Some(n + 1) } else { None }).max().unwrap_or(0);
        StarId::Added { n }
    }

    /// The number the next added cluster takes.
    #[must_use]
    pub fn next_cluster(&self) -> u32 {
        self.clusters.iter().map(|c| c.n + 1).max().unwrap_or(0)
    }

    /// The edited record of `id`, if it was added or changed.
    #[must_use]
    pub fn record(&self, id: &StarId) -> Option<&StarRecord> {
        self.add.iter().chain(&self.change).find(|r| r.id == *id)
    }

    /// Store `record`: an added star is replaced in place, a generated one becomes (or
    /// replaces) a change, and a removed star comes back.
    pub fn put(&mut self, record: StarRecord) {
        self.remove.retain(|id| *id != record.id);
        let list = if matches!(record.id, StarId::Added { .. }) { &mut self.add } else { &mut self.change };
        match list.iter_mut().find(|r| r.id == record.id) {
            Some(r) => *r = record,
            None => list.push(record),
        }
    }

    /// Take `id` away: an added star is deleted, a generated one removed.
    pub fn take_away(&mut self, id: StarId) {
        self.add.retain(|r| r.id != id);
        self.change.retain(|r| r.id != id);
        if !matches!(id, StarId::Added { .. }) && !self.remove.contains(&id) {
            self.remove.push(id);
        }
    }
}

/// A named figure of lines between stars, as one observer sees them.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Constellation {
    pub name: String,
    pub lines: Vec<[StarId; 2]>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn galaxy_ids_read_back() {
        let id = GalaxyId { bin: 6, cube: [-3, 0, 12], n: 2, member: 17 };
        assert_eq!(id.to_string().parse::<GalaxyId>(), Ok(id));
        assert_eq!(serde_json::from_str::<GalaxyId>(&serde_json::to_string(&id).unwrap()).unwrap(), id);
        assert!("h6:1,2:3.4".parse::<GalaxyId>().is_err());
        assert_ne!(id.seed(Seed(1)), GalaxyId { member: 18, ..id }.seed(Seed(1)));
    }

    #[test]
    fn identities_round_trip_as_text() {
        let ids = [StarId::Field { bin: 17, cube: [1025, -3, 0], n: 4 }, StarId::Group { bin: 5, cube: [-2, 0, 14], event: 9, n: 3 }, StarId::Added { n: 12 }, StarId::Cluster { cluster: 2, n: 15 }];
        for id in ids {
            assert_eq!(id.to_string().parse::<StarId>(), Ok(id));
            let json = serde_json::to_string(&id).unwrap();
            assert_eq!(serde_json::from_str::<StarId>(&json).unwrap(), id);
        }
        assert_eq!(ids[0].to_string(), "f17:1025,-3,0:4");
        assert!("x1".parse::<StarId>().is_err() && "f1:2,3:4".parse::<StarId>().is_err());
        // Each identity seeds its own system.
        let seeds: std::collections::HashSet<u64> = ids.iter().map(|id| id.seed(Seed(1)).0).collect();
        assert_eq!(seeds.len(), 4);
    }

    #[test]
    fn edits_keep_one_record_per_star() {
        let kind = StarKind { initial_mass: 1.0, age_gyr: 4.6, metallicity_dex: 0.0, companion_mass: None };
        let field = StarId::Field { bin: 9, cube: [0, 0, 0], n: 1 };
        let mut e = StarEdits::default();
        e.take_away(field);
        e.put(StarRecord { id: field, name: Some("Moved".into()), position_pc: [1.0, 0.0, 0.0], kind });
        e.put(StarRecord { id: field, name: Some("Moved again".into()), position_pc: [2.0, 0.0, 0.0], kind });
        assert!(e.remove.is_empty() && e.change.len() == 1 && e.change[0].position_pc[0] == 2.0);
        let added = e.next_added();
        e.put(StarRecord { id: added, name: None, position_pc: [0.0; 3], kind });
        assert_eq!(e.next_added(), StarId::Added { n: 1 });
        e.take_away(added);
        e.take_away(field);
        assert!(e.add.is_empty() && e.change.is_empty() && e.remove == vec![field]);
    }
}
