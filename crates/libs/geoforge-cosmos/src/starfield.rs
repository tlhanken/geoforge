//! A galaxy's stars with a world's edits over them, and where a system sits among them.
//!
//! Every query reads the generated stars ([`crate::sectors`]) and applies the edits
//! ([`StarEdits`]): removed stars go, changed ones are replaced by their stored records
//! wherever those now lie, added stars and the members of added clusters join. An added
//! cluster's members are drawn from its own seed, alive at its age, inside a Plummer sphere.

use crate::draw::Draw;
use crate::galaxy::{presets, Galaxy};
use crate::sectors::tables::Table;
use crate::sectors::{GalacticStar, Sectors};
use geoforge_core::body::SolarSystem;
use geoforge_core::galactic::{ClusterEdit, StarEdits, StarId, StarKind};
use geoforge_core::seed::Seed;

/// Most members an added cluster may have.
const MAX_MEMBERS: u32 = 200_000;

/// A galaxy's stars as a world sees them.
pub struct Starfield {
    sectors: Sectors,
    edits: StarEdits,
    /// Members of the added clusters.
    members: Vec<GalacticStar>,
}

impl Starfield {
    /// The stars of `galaxy` seeded `seed`, with `edits` over them.
    #[must_use]
    pub fn new(galaxy: &Galaxy, seed: Seed, edits: StarEdits) -> Self {
        let sectors = Sectors::new(galaxy, seed);
        let members = edits.clusters.iter().flat_map(|c| cluster_members(galaxy, seed, c)).collect();
        Self { sectors, edits, members }
    }

    /// The stars around `system`: its site's galaxy (or `galaxy` if given) and galaxy seed,
    /// with `edits`.
    #[must_use]
    pub fn of_system(system: &SolarSystem, galaxy: Option<Galaxy>, edits: StarEdits) -> Self {
        let site = system.root.site.as_ref();
        let galaxy = galaxy.unwrap_or_else(|| site.map_or_else(Galaxy::milky_way, |s| presets::by_name(&s.galaxy)));
        let seed = site.and_then(|s| s.galaxy_seed).unwrap_or_else(|| system.root.seed.child(b"galaxy"));
        Self::new(&galaxy, seed, edits)
    }

    #[must_use]
    pub fn sectors(&self) -> &Sectors {
        &self.sectors
    }

    #[must_use]
    pub fn galaxy(&self) -> &Galaxy {
        self.sectors.galaxy()
    }

    #[must_use]
    pub fn edits(&self) -> &StarEdits {
        &self.edits
    }

    #[must_use]
    pub fn galaxy_seed(&self) -> Seed {
        self.sectors.seed()
    }

    /// Stars within `radius_pc` of `p` brighter than absolute V `faintest`.
    #[must_use]
    pub fn near(&self, p: [f64; 3], radius_pc: f64, faintest: f64) -> Vec<GalacticStar> {
        self.apply(self.sectors.near(p, radius_pc, faintest), |s| s.distance_pc(p) <= radius_pc && s.absolute_v <= faintest)
    }

    /// Stars that could be brighter than `limit` seen from `p` before dust, out to `max_pc`.
    #[must_use]
    pub fn visible(&self, p: [f64; 3], limit: f64, max_pc: f64) -> Vec<GalacticStar> {
        self.apply(self.sectors.visible(p, limit, max_pc), |s| {
            let d = s.distance_pc(p).max(1e-3);
            d <= max_pc && s.absolute_v + 5.0 * (d / 10.0).log10() <= limit
        })
    }

    /// The star `id` as the world has it, if it exists.
    #[must_use]
    pub fn star(&self, id: &StarId) -> Option<GalacticStar> {
        if self.edits.remove.contains(id) {
            return None;
        }
        if let Some(r) = self.edits.record(id) {
            return Some(GalacticStar::new(r.id, r.position_pc, r.kind, None));
        }
        match id {
            StarId::Cluster { .. } => self.members.iter().find(|m| m.id == *id).cloned(),
            StarId::Added { .. } => None,
            _ => self.sectors.star(id),
        }
    }

    /// The name of star `id`: the edited one, else one drawn from its identity.
    #[must_use]
    pub fn name_of(&self, id: &StarId) -> String {
        if let Some(name) = self.edits.record(id).and_then(|r| r.name.clone()) {
            return name;
        }
        if let StarId::Cluster { cluster, n } = id {
            if let Some(c) = self.edits.clusters.iter().find(|c| c.n == *cluster) {
                return format!("{} {}", c.name, n + 1);
            }
        }
        geoforge_names::name(id.seed(self.galaxy_seed()).0)
    }

    /// Generated stars with the edits applied; `keep` selects the edited and added ones.
    fn apply(&self, generated: Vec<GalacticStar>, keep: impl Fn(&GalacticStar) -> bool) -> Vec<GalacticStar> {
        let e = &self.edits;
        let replaced = |id: &StarId| e.remove.contains(id) || e.change.iter().any(|r| r.id == *id);
        let mut out: Vec<GalacticStar> = generated.into_iter().filter(|s| !replaced(&s.id)).collect();
        let edited = e.change.iter().chain(&e.add).map(|r| GalacticStar::new(r.id, r.position_pc, r.kind, None));
        let members = self.members.iter().filter(|m| !replaced(&m.id)).cloned();
        out.extend(edited.chain(members).filter(|s| keep(s)));
        out
    }
}

/// The members of added cluster `c` in `galaxy` seeded `seed`.
fn cluster_members(galaxy: &Galaxy, seed: Seed, c: &ClusterEdit) -> Vec<GalacticStar> {
    let seed = seed.child(b"added cluster").child_u64(b"n", u64::from(c.n));
    let table = Table::at_age(&galaxy.imf, c.age_gyr, &|_, s| c.metallicity_dex + 0.05 * s, seed.child(b"kinds"), 6000);
    let mut d = Draw::new(seed.child(b"members"));
    (0..c.members.min(MAX_MEMBERS))
        .filter_map(|n| {
            let (_, k) = table.pick_any(d.u(), d.u())?;
            let r = c.radius_pc / (d.u().max(1e-12).powf(-2.0 / 3.0) - 1.0).max(1e-12).sqrt();
            let u = d.direction();
            let position = [0, 1, 2].map(|i| c.centre_pc[i] + u[i] * r.min(20.0 * c.radius_pc));
            let kind = StarKind { initial_mass: k.initial_mass, age_gyr: c.age_gyr, metallicity_dex: c.metallicity_dex, companion_mass: k.companion };
            Some(GalacticStar::new(StarId::Cluster { cluster: c.n, n }, position, kind, None))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use geoforge_core::galactic::StarRecord;

    #[test]
    fn edits_apply_to_every_query() {
        let g = Galaxy::milky_way();
        let sun = [8200.0, 0.0, 25.0];
        let plain = Starfield::new(&g, Seed(3), StarEdits::default());
        let stars = plain.near(sun, 12.0, 14.0);
        let (gone, moved) = (stars[0].clone(), stars[1].clone());
        let mut e = StarEdits::default();
        e.take_away(gone.id);
        let far = [sun[0] + 500.0, sun[1], sun[2]];
        e.put(StarRecord { id: moved.id, name: Some("Wanderer".into()), position_pc: far, kind: moved.kind });
        let added = e.next_added();
        e.put(StarRecord {
            id: added,
            name: Some("Beacon".into()),
            position_pc: [sun[0] + 3.0, sun[1], sun[2]],
            kind: StarKind { initial_mass: 2.0, age_gyr: 0.3, metallicity_dex: 0.0, companion_mass: None },
        });
        e.clusters.push(ClusterEdit { n: 0, name: "Hearth".into(), centre_pc: [sun[0], sun[1] + 6.0, sun[2]], radius_pc: 1.0, members: 300, age_gyr: 0.1, metallicity_dex: 0.0 });
        let f = Starfield::new(&g, Seed(3), e);
        let now = f.near(sun, 12.0, 14.0);
        assert!(!now.iter().any(|s| s.id == gone.id || s.id == moved.id));
        assert!(now.iter().any(|s| s.id == added) && f.name_of(&added) == "Beacon");
        let members = f.near(sun, 12.0, 30.0).iter().filter(|s| matches!(s.id, StarId::Cluster { .. })).count();
        assert!((250..=300).contains(&members), "{members} members");
        assert!(now.iter().any(|s| matches!(s.id, StarId::Cluster { .. })));
        assert!(f.near(far, 1.0, 20.0).iter().any(|s| s.id == moved.id), "the moved star is at its new place");
        assert_eq!(f.star(&moved.id).map(|s| s.position_pc), Some(far));
        assert!(f.star(&gone.id).is_none() && plain.star(&gone.id).is_some());
        assert_eq!(f.name_of(&StarId::Cluster { cluster: 0, n: 4 }), "Hearth 5");
    }
}
