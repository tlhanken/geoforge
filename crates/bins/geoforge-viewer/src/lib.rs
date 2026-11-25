//! The viewer: a universe directory seen at any level, generated as it is looked at.
//!
//! [`Viewer`] answers by address: the scene there ([`scene`]: the universe's galaxies, a
//! galaxy's shape and the systems made in it, the stars around one, a system's orbits, a
//! body's maps), the maps of a body as images (true colour for every body; every layer and
//! view of one generated here), a search over what is made or near it ([`search`]), and
//! making a body or generating its layers. [`server`] puts it behind a local web server for
//! the page in `web/` (three.js draws it); everything the page shows is worked out here.

pub mod scene;
pub mod search;
pub mod server;

use geoforge_core::address::{Address, GalaxyRef, Level, StarRef};
use geoforge_cosmos::resolve::{self, GalaxyAt};
use geoforge_cosmos::starfield::Starfield;
use geoforge_store::{UniverseDir, World};
use std::collections::HashMap;
use std::path::Path;
use std::sync::{Arc, Mutex, RwLock};

/// Why the viewer cannot answer.
#[derive(Debug)]
pub struct ViewError(pub String);

impl std::fmt::Display for ViewError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for ViewError {}

macro_rules! from_error {
    ($($t:ty),+) => {$(
        impl From<$t> for ViewError {
            fn from(e: $t) -> Self {
                Self(e.to_string())
            }
        }
    )+};
}
from_error!(geoforge_store::StoreError, geoforge_pipeline::PipelineError, resolve::Unresolved, image::ImageError, String);

pub type Result<T> = std::result::Result<T, ViewError>;

/// A galaxy and its stars (with the universe's edits), worked out once.
pub struct Galaxy {
    pub at: GalaxyAt,
    pub field: Starfield,
    /// Its shape as drawn stars, old and young, once asked for ([`scene`]).
    pub cloud: std::sync::OnceLock<(Vec<f32>, Vec<f32>)>,
}

/// A universe directory being viewed.
pub struct Viewer {
    universe: RwLock<UniverseDir>,
    galaxies: Mutex<HashMap<Address, Arc<Galaxy>>>,
    /// Generation of bodies' layers under way or done: what it says now.
    pub jobs: Mutex<HashMap<Address, String>>,
    /// The names search looks among, gathered on first need ([`search`]).
    index: Mutex<Option<Vec<search::Hit>>>,
}

impl Viewer {
    /// # Errors
    /// If `dir` holds no universe.
    pub fn open(dir: &Path) -> Result<Self> {
        Ok(Self { universe: RwLock::new(UniverseDir::open(dir)?), galaxies: Mutex::new(HashMap::new()), jobs: Mutex::new(HashMap::new()), index: Mutex::new(None) })
    }

    /// The universe (an address at its level) and the focus body, if any.
    #[must_use]
    pub fn start(&self) -> (Address, Option<Address>) {
        let u = self.universe.read().expect("universe lock");
        (u.doc.universe, u.doc.focus)
    }

    /// `text` as an address of this universe (its universe part may be left out).
    ///
    /// # Errors
    /// If it is no address, or one of another universe.
    pub fn parse(&self, text: &str) -> Result<Address> {
        let universe = self.start().0;
        let a = text.parse::<Address>().or_else(|e| format!("{universe}/{}", text.trim_start_matches('/')).parse::<Address>().map_err(|_| e))?;
        if a.universe != universe.universe {
            return Err(ViewError(format!("{a} is not in this universe ({universe})")));
        }
        Ok(a)
    }

    /// The galaxy at `address` with its stars, as the universe keeps them.
    ///
    /// # Errors
    /// If the address names no galaxy.
    pub fn galaxy(&self, address: &Address) -> Result<Arc<Galaxy>> {
        let key = address.up_to(Level::Galaxy);
        if let Some(g) = self.galaxies.lock().expect("galaxy cache").get(&key) {
            return Ok(Arc::clone(g));
        }
        let u = self.universe.read().expect("universe lock");
        let mut at = resolve::galaxy(&key)?;
        if let Some(model) = u.galaxy_model(&key)? {
            at.model = model;
        }
        let field = at.starfield(u.star_edits(&key)?);
        drop(u);
        let g = Arc::new(Galaxy { at, field, cloud: std::sync::OnceLock::new() });
        self.galaxies.lock().expect("galaxy cache").insert(key, Arc::clone(&g));
        Ok(g)
    }

    /// The world made for body `address`, if any.
    ///
    /// # Errors
    /// On unreadable manifests.
    pub fn world(&self, address: &Address) -> Result<Option<World>> {
        Ok(self.universe.read().expect("universe lock").world(address)?)
    }

    /// Every body made in the universe.
    ///
    /// # Errors
    /// On unreadable directories.
    pub fn worlds(&self) -> Result<Vec<World>> {
        Ok(self.universe.read().expect("universe lock").worlds()?)
    }

    /// Make the world of body `address` (completed to a body by the default rule), at a
    /// grid of `width` cells round, and make it the universe's focus.
    ///
    /// # Errors
    /// If the address names nothing or the world cannot be written.
    pub fn make(&self, address: &Address, width: u32) -> Result<World> {
        let mut u = self.universe.write().expect("universe lock");
        let w = geoforge_pipeline::create_at(&mut u, address, None, geoforge_core::grid::GridSpec::Equirect { width, height: width / 2 })?;
        *self.index.lock().expect("index") = None;
        Ok(w)
    }

    /// Generate every layer of the world of body `address` (made first if need be), noting
    /// progress in [`Self::jobs`].
    ///
    /// # Errors
    /// If the body cannot be made or a layer fails.
    pub fn generate(&self, address: &Address, width: u32) -> Result<()> {
        let note = |s: String| {
            self.jobs.lock().expect("jobs").insert(*address, s);
        };
        let mut w = match self.world(address)? {
            Some(w) => w,
            None => self.make(address, width)?,
        };
        note("starting".into());
        let r = geoforge_pipeline::generate(&mut w, &geoforge_pipeline::GenerateOptions::default(), |r| note(format!("{} done", r.layer.name())));
        note(match &r {
            Ok(_) => "done".into(),
            Err(e) => format!("failed: {e}"),
        });
        r.map(|_| ()).map_err(Into::into)
    }

    /// The system at `address`: as stored, if a body of it is made; else generated from its
    /// star, with the universe's edits.
    ///
    /// # Errors
    /// If the address stops above a system or names nothing.
    pub fn system(&self, address: &Address) -> Result<geoforge_cosmos::system::Generated> {
        if let Some(w) = self.made_in(address)? {
            if let Some(system) = w.load_system()? {
                return Ok(geoforge_cosmos::system::Generated { focus: w.manifest.focus, system });
            }
        }
        let g = self.galaxy(address)?;
        Ok(resolve::system(address, &g.field)?)
    }

    /// Any world made in `address`'s system.
    fn made_in(&self, address: &Address) -> Result<Option<World>> {
        let sys = address.up_to(Level::System);
        Ok(self.worlds()?.into_iter().find(|w| w.manifest.address.map(|a| a.up_to(Level::System)) == Some(sys)))
    }

    /// The universe's default galaxy, or the focus body's.
    #[must_use]
    pub fn home(&self) -> Address {
        let (universe, focus) = self.start();
        focus.map_or_else(
            || match universe.universe {
                geoforge_core::galactic::UniverseSite::Real => universe.galaxy(GalaxyRef::MilkyWay),
                geoforge_core::galactic::UniverseSite::Seeded(s) => universe.galaxy(GalaxyRef::Generated(geoforge_cosmos::universe::Universe::new(s).home(&geoforge_cosmos::Galaxy::milky_way()).id)),
            },
            |f| f.up_to(Level::Galaxy),
        )
    }

    /// The star a galaxy view centres on: the focus body's, else the default star's.
    #[must_use]
    pub fn home_star(&self, galaxy: &Address) -> Option<Address> {
        let (_, focus) = self.start();
        focus
            .filter(|f| f.up_to(Level::Galaxy) == galaxy.up_to(Level::Galaxy))
            .map(|f| f.up_to(Level::System))
            .or_else(|| resolve::complete(&galaxy.up_to(Level::Galaxy)).ok().map(|a| a.up_to(Level::System)))
    }
}

/// The star of a system address, if it is one of the galaxy's generated or added stars.
#[must_use]
pub fn star_id(address: &Address) -> Option<geoforge_core::galactic::StarId> {
    match address.star {
        Some(StarRef::Star(id)) => Some(id),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scene::Scene;
    use geoforge_core::Seed;

    /// Every level of universe 7 around its default body: the scenes, a map, the way down,
    /// search, and an address given without its universe.
    #[test]
    fn a_universe_seen_level_by_level() {
        let dir = std::env::temp_dir().join(format!("geoforge-viewer-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        let mut u = UniverseDir::create(&dir, Address::universe(Seed(7))).unwrap();
        let w = geoforge_pipeline::create_at(&mut u, &Address::universe(Seed(7)), None, geoforge_core::grid::GridSpec::Equirect { width: 64, height: 32 }).unwrap();
        let body = w.manifest.address.unwrap();
        let v = Viewer::open(&dir).unwrap();
        assert_eq!(v.home(), body.up_to(Level::Galaxy));

        let Scene::Universe(un) = v.scene(&body.up_to(Level::Universe), false).unwrap() else { panic!("universe") };
        assert!(un.galaxies.iter().any(|g| g.home) && un.galaxies.len() > 10, "{} galaxies", un.galaxies.len());
        let Scene::Galaxy(g) = v.scene(&body.up_to(Level::Galaxy), false).unwrap() else { panic!("galaxy") };
        assert!(g.cloud.len() == 6 * 20_000 && !g.young.is_empty());
        assert!(g.systems.iter().any(|s| s.made && s.address == body.up_to(Level::System).to_string()));
        let Scene::Stars(st) = v.scene(&body, true).unwrap() else { panic!("stars") };
        assert_eq!(st.stars.iter().filter(|s| s.centre).count(), 1);
        let Scene::System(sys) = v.scene(&body.up_to(Level::System), false).unwrap() else { panic!("system") };
        assert!(sys.bodies.iter().any(|b| b.made && b.focus && b.orbit_au.len() == 96));
        let Scene::Body(b) = v.scene(&body, false).unwrap() else { panic!("body") };
        assert!(b.made && b.layers == ["appearance"], "{:?}", b.layers);

        let png = v.texture(&body, "appearance", 128).unwrap();
        let img = image::load_from_memory(&png).unwrap();
        assert_eq!((img.width(), img.height()), (128, 64));
        assert!(v.texture(&body, "koppen", 128).is_err(), "not generated yet");
        let other = sys.bodies.iter().find(|x| !x.made).unwrap();
        assert!(v.texture(&other.address.parse().unwrap(), "appearance", 64).is_ok(), "any body has a face");

        let crumbs = v.crumbs(&body).unwrap();
        assert_eq!(crumbs.len(), 5);
        assert_eq!(crumbs[4].name, w.manifest.focus_name);
        let hits = v.search(&w.manifest.focus_name).unwrap();
        assert_eq!(hits[0].address, body.to_string());
        let short = body.to_string().split_once('/').unwrap().1.to_string();
        assert_eq!(v.parse(&short).unwrap(), body);
        assert!(v.parse("8").is_err(), "another universe");
        let _ = std::fs::remove_dir_all(&dir);
    }
}
