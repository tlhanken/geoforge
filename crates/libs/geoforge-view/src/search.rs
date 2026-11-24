//! Finding things by name or address.
//!
//! A universe has no end to search, so names are looked for where they are already known:
//! the bodies made and their systems' bodies, the stars within 15 pc of each made system,
//! and the galaxies within 10 Mpc of the home galaxy. Anything else is reached by typing its
//! address.

use crate::{Result, Viewer};
use geoforge_core::address::{GalaxyRef, Level, StarRef};
use geoforge_core::body::NodeKind;
use geoforge_core::galactic::UniverseSite;
use serde::Serialize;

/// Hits returned at most.
const MOST: usize = 40;

/// Something found.
#[derive(Debug, Clone, Serialize, PartialEq)]
pub struct Hit {
    pub name: String,
    /// galaxy, star, body or address.
    pub kind: String,
    pub address: String,
    /// What it is, in a few words.
    pub note: String,
}

impl Viewer {
    /// What `query` names: the address it is, if it is one, then the known names containing
    /// it (case aside), whole names first, then those it begins.
    ///
    /// # Errors
    /// On unreadable worlds.
    pub fn search(&self, query: &str) -> Result<Vec<Hit>> {
        let q = query.trim().to_lowercase();
        if q.is_empty() {
            return Ok(Vec::new());
        }
        let mut hits = Vec::new();
        if let Ok(a) = self.parse(query.trim()) {
            hits.push(Hit { name: a.to_string(), kind: "address".into(), address: a.to_string(), note: format!("{:?}", a.level()).to_lowercase() });
        }
        let mut known = {
            let mut index = self.index.lock().expect("index");
            if index.is_none() {
                *index = Some(self.known()?);
            }
            index.clone().unwrap_or_default()
        };
        known.retain(|h| h.name.to_lowercase().contains(&q));
        known.sort_by_cached_key(|h| {
            let n = h.name.to_lowercase();
            (n != q, !n.starts_with(&q), n)
        });
        known.dedup_by(|a, b| a.address == b.address);
        hits.extend(known.into_iter().take(MOST));
        Ok(hits)
    }

    /// Every name the viewer knows now.
    fn known(&self) -> Result<Vec<Hit>> {
        let mut out = Vec::new();
        let (universe, _) = self.start();
        if let Ok(crate::scene::Scene::Universe(u)) = self.scene(&universe, false) {
            out.extend(u.galaxies.into_iter().filter_map(|g| Some(Hit { name: g.name, kind: "galaxy".into(), address: g.address?, note: g.kind })));
        }
        if universe.universe == UniverseSite::Real {
            out.push(Hit { name: "Sun".into(), kind: "star".into(), address: universe.galaxy(GalaxyRef::MilkyWay).star(StarRef::Sun).to_string(), note: "the Solar System".into() });
        }
        let mut systems = Vec::new();
        for w in self.worlds()? {
            let Some(a) = w.manifest.address else { continue };
            if let Some(sys) = w.load_system()? {
                out.push(Hit { name: sys.root.name.clone(), kind: "star".into(), address: a.up_to(Level::System).to_string(), note: "a system made here".into() });
                for n in &sys.nodes {
                    if let NodeKind::Body(_) = n.node {
                        let at = a.up_to(Level::System).body(n.key);
                        let made = at == a;
                        out.push(Hit { name: n.name.clone(), kind: "body".into(), address: at.to_string(), note: if made { "made".into() } else { format!("in {}", sys.root.name) } });
                    }
                }
            }
            systems.push(a);
        }
        for a in systems {
            if let Ok(crate::scene::Scene::Stars(s)) = self.scene(&a.up_to(Level::System), true) {
                out.extend(s.stars.into_iter().map(|st| Hit { name: st.name, kind: "star".into(), address: st.address, note: st.spectral }));
            }
        }
        Ok(out)
    }
}
