//! The universe and world directories. Your sources and edits are kept; generated products
//! are a cache that can always be rebuilt.
//!
//! A universe directory holds what of a universe has been looked at or edited; everything
//! else regenerates from its seed. Each body made there is a world, at its address:
//!
//! ```text
//! universe/
//!   universe.json                    the universe (its seed, or real) and the focus body
//!   galaxies/<galaxy>/
//!     galaxy.json                    the galaxy's model, if edited (else generated)
//!     stars.json                     edits over its generated stars
//!   bodies/<galaxy>/<star>/<body>/   a world (below), one per body made
//! ```
//!
//! A world (on its own, or a body of a universe):
//!
//! ```text
//! world/
//!   world.toml                       manifest: seed, grid, focus, address, sources, provenance
//!   system.json                      the solar system (absent for a hand-written body)
//!   bodies/<body>/
//!     body.json                      BodyContext (the handoff to the surface pipeline)
//!     products/<product>.json|.bin   canonical stored products
//!     edits/<product>.*              authored data and patches (never auto-overwritten)
//!   galaxy/                          a world on its own: its galaxy.json and stars.json
//!   cache/baselines/                 last edit-view export of each product (for diffs)
//! ```

pub mod cache;
pub mod codec;
pub mod edit;
pub mod polar;
pub mod raster;
pub mod render;
pub mod report;

use geoforge_core::address::{Address, Level};
use geoforge_core::body::{BodyContext, SolarSystem};
use geoforge_core::grid::GridSpec;
use geoforge_core::key::Key;
use geoforge_core::product::{Hash, ProductData, ProductId, Source};
use geoforge_core::Seed;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

#[derive(Debug, thiserror::Error)]
pub enum StoreError {
    #[error("io: {0}")]
    Io(#[from] std::io::Error),
    #[error("json: {0}")]
    Json(#[from] serde_json::Error),
    #[error("toml: {0}")]
    TomlDe(#[from] toml::de::Error),
    #[error("toml: {0}")]
    TomlSer(#[from] toml::ser::Error),
    #[error("image: {0}")]
    Image(#[from] image::ImageError),
    #[error("{0}")]
    Core(#[from] geoforge_core::Error),
    #[error("{0}")]
    Invalid(String),
}

pub type Result<T> = std::result::Result<T, StoreError>;

/// In a galaxy's directory: its model and the edits over its stars.
const GALAXY_MODEL: &str = "galaxy.json";
const STAR_EDITS: &str = "stars.json";
const CONSTELLATIONS: &str = "constellations.json";
const UNIVERSE: &str = "universe.json";
/// In a universe directory: its fictional materials, written by hand.
const EXOTICS: &str = "exotics.toml";
/// What a new universe's `exotics.toml` says: how to write one, with an example left off.
const EXOTICS_TEMPLATE: &str = r##"# Fictional materials of this universe. They appear only as deposits and markers, flagged
# fictional, and change nothing else that is made. Each is drawn from its own stream (keyed
# by its name), so adding or editing one moves no other. Remove the '#' marks to try this one.
#
# [[material]]
# name = "Xylarium"
# note = "a metastable superheavy alloy"
# colour = "#7fffd4"
#
# [[material.found]]
# in = "asteroid"             # asteroid | world | nebula
# kinds = ["M"]               # asteroid: S C M V P D icy; world: rocky icy ocean giant small;
#                             # nebula: emission reflection dark planetary supernova_remnant
# stars = ["neutron_star"]    # only round these stars (main_sequence, subgiant, giant,
#                             # supergiant, white_dwarf, neutron_star, black_hole, brown_dwarf)
# chance = 0.002              # per asteroid, world or nebula of those kinds
"##;

/// `world.toml`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Manifest {
    pub schema_version: u32,
    pub code_version: String,
    pub root_seed: Seed,
    /// Key and name of the body whose surface this world builds.
    pub focus: Key,
    pub focus_name: String,
    pub grid: GridSpec,
    /// Layer parameters by layer name (defaults are written out so they can be edited).
    #[serde(default)]
    pub params: BTreeMap<String, serde_json::Value>,
    /// Sources by product name; a product not listed is generated.
    #[serde(default)]
    pub sources: BTreeMap<String, Source>,
    /// For each stored product, the hash of everything it was built from.
    #[serde(default)]
    pub provenance: BTreeMap<String, Hash>,
    /// Layers that do not apply to the focus body, with the reason.
    #[serde(default)]
    pub not_applicable: BTreeMap<String, String>,
    /// Where the focus body is in its universe; `None` for a world on its own.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub address: Option<Address>,
}

impl Manifest {
    #[must_use]
    pub fn source(&self, p: ProductId) -> Source {
        self.sources.get(p.name()).cloned().unwrap_or_else(Source::generate)
    }
}

/// An open world directory.
#[derive(Debug, Clone)]
pub struct World {
    pub root: PathBuf,
    pub manifest: Manifest,
    /// The universe directory it is a body of, if any.
    pub universe: Option<PathBuf>,
}

/// `universe.json`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct UniverseDoc {
    /// The universe (an address at its level).
    pub universe: Address,
    /// The body commands work on when given the universe directory.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub focus: Option<Address>,
}

/// An open universe directory.
#[derive(Debug, Clone)]
pub struct UniverseDir {
    pub root: PathBuf,
    pub doc: UniverseDoc,
}

impl UniverseDir {
    /// Start a universe directory for `universe` (any address in it: it is cut to the
    /// universe).
    ///
    /// # Errors
    /// If the directory already holds a universe or a world, or on IO errors.
    pub fn create(root: &Path, universe: Address) -> Result<Self> {
        if root.join(UNIVERSE).exists() || root.join("world.toml").exists() {
            return Err(StoreError::Invalid(format!("{} already holds a universe or a world", root.display())));
        }
        std::fs::create_dir_all(root)?;
        let u = Self { root: root.to_path_buf(), doc: UniverseDoc { universe: universe.up_to(Level::Universe), focus: None } };
        u.save()?;
        if !root.join(EXOTICS).exists() {
            std::fs::write(root.join(EXOTICS), EXOTICS_TEMPLATE)?;
        }
        Ok(u)
    }

    /// The universe's fictional materials (`exotics.toml`; none if there is no such file).
    ///
    /// # Errors
    /// If the file is unreadable, not TOML of the list's shape, or fails its checks.
    pub fn exotics(&self) -> Result<geoforge_core::exotics::Exotics> {
        let path = self.root.join(EXOTICS);
        if !path.exists() {
            return Ok(geoforge_core::exotics::Exotics::default());
        }
        let list: geoforge_core::exotics::Exotics = toml::from_str(&std::fs::read_to_string(&path)?).map_err(|e| StoreError::Invalid(format!("{}: {e}", path.display())))?;
        list.validate()?;
        Ok(list)
    }

    /// # Errors
    /// If there is no universe at `root`.
    pub fn open(root: &Path) -> Result<Self> {
        let text = std::fs::read_to_string(root.join(UNIVERSE)).map_err(|e| StoreError::Invalid(format!("no universe at {}: {e}", root.display())))?;
        Ok(Self { root: root.to_path_buf(), doc: serde_json::from_str(&text)? })
    }

    /// Whether `root` holds a universe.
    #[must_use]
    pub fn is_at(root: &Path) -> bool {
        root.join(UNIVERSE).exists()
    }

    fn save(&self) -> Result<()> {
        std::fs::write(self.root.join(UNIVERSE), serde_json::to_string_pretty(&self.doc)?)?;
        Ok(())
    }

    /// Make body `address` the one commands work on.
    ///
    /// # Errors
    /// On IO errors.
    pub fn set_focus(&mut self, address: Address) -> Result<()> {
        self.doc.focus = Some(address);
        self.save()
    }

    /// Where the world of body `address` is (whether or not it has been made).
    #[must_use]
    pub fn body_root(&self, address: &Address) -> PathBuf {
        address.segments().iter().skip(1).fold(self.root.join("bodies"), |p, s| p.join(s))
    }

    /// Where the documents of `address`'s galaxy are.
    #[must_use]
    pub fn galaxy_dir(&self, address: &Address) -> PathBuf {
        galaxy_dir_in(&self.root, address)
    }

    /// The edited model of `address`'s galaxy, if it has one.
    ///
    /// # Errors
    /// On IO or parse errors.
    pub fn galaxy_model<T: serde::de::DeserializeOwned>(&self, address: &Address) -> Result<Option<T>> {
        read_json(&self.galaxy_dir(address).join(GALAXY_MODEL))
    }

    /// The edits over the stars of `address`'s galaxy.
    ///
    /// # Errors
    /// On IO or parse errors.
    pub fn star_edits(&self, address: &Address) -> Result<geoforge_core::galactic::StarEdits> {
        Ok(read_json(&self.galaxy_dir(address).join(STAR_EDITS))?.unwrap_or_default())
    }

    /// The world of body `address`, if it has been made.
    ///
    /// # Errors
    /// If its manifest is unreadable.
    pub fn world(&self, address: &Address) -> Result<Option<World>> {
        let root = self.body_root(address);
        if root.join("world.toml").exists() {
            World::open(&root).map(Some)
        } else {
            Ok(None)
        }
    }

    /// Every body made so far, by address.
    ///
    /// # Errors
    /// On unreadable directories or manifests.
    pub fn worlds(&self) -> Result<Vec<World>> {
        fn walk(dir: &Path, out: &mut Vec<PathBuf>) -> std::io::Result<()> {
            if dir.join("world.toml").exists() {
                out.push(dir.to_path_buf());
                return Ok(());
            }
            if dir.is_dir() {
                for e in std::fs::read_dir(dir)? {
                    walk(&e?.path(), out)?;
                }
            }
            Ok(())
        }
        let mut roots = Vec::new();
        walk(&self.root.join("bodies"), &mut roots)?;
        roots.sort();
        roots.iter().map(|r| World::open(r)).collect()
    }
}

fn galaxy_dir_in(universe: &Path, address: &Address) -> PathBuf {
    universe.join("galaxies").join(address.segments().get(1).map_or("_", String::as_str))
}

fn read_json<T: serde::de::DeserializeOwned>(p: &Path) -> Result<Option<T>> {
    if !p.exists() {
        return Ok(None);
    }
    Ok(Some(serde_json::from_str(&std::fs::read_to_string(p)?)?))
}

fn write_json<T: Serialize>(p: &Path, doc: &T) -> Result<()> {
    if let Some(dir) = p.parent() {
        std::fs::create_dir_all(dir)?;
    }
    std::fs::write(p, serde_json::to_string_pretty(doc)?)?;
    Ok(())
}

/// The universe directory holding the world at `root`: the nearest ancestor with a
/// `universe.json`.
fn universe_above(root: &Path) -> Option<PathBuf> {
    root.ancestors().skip(1).find(|a| UniverseDir::is_at(a)).map(Path::to_path_buf)
}

/// Lower-case, filesystem-safe form of a name.
#[must_use]
pub fn slug(s: &str) -> String {
    let mut out: String = s.chars().map(|c| if c.is_ascii_alphanumeric() { c.to_ascii_lowercase() } else { '-' }).collect();
    while out.contains("--") {
        out = out.replace("--", "-");
    }
    out.trim_matches('-').to_string()
}

/// A world's body, grid and stored products, loaded once for a command.
pub struct Snapshot {
    pub body: BodyContext,
    pub grid: geoforge_core::grid::AnyGrid,
    pub products: geoforge_core::product::Products,
}

/// blake3 of arbitrary parts, as a [`Hash`].
#[must_use]
pub fn hash_parts(parts: &[&[u8]]) -> Hash {
    let mut h = blake3::Hasher::new();
    for p in parts {
        h.update(&(p.len() as u64).to_le_bytes());
        h.update(p);
    }
    Hash(*h.finalize().as_bytes())
}

impl World {
    /// Create a new world directory.
    ///
    /// # Errors
    /// If the directory already holds a world, or on IO errors.
    pub fn create(root: &Path, manifest: Manifest, system: Option<&SolarSystem>, body: &BodyContext) -> Result<Self> {
        if root.join("world.toml").exists() {
            return Err(StoreError::Invalid(format!("{} already holds a world", root.display())));
        }
        std::fs::create_dir_all(root)?;
        let universe = manifest.address.and_then(|_| universe_above(root));
        let w = Self { root: root.to_path_buf(), manifest, universe };
        if let Some(s) = system {
            std::fs::write(w.root.join("system.json"), serde_json::to_string_pretty(s)?)?;
        }
        std::fs::create_dir_all(w.body_dir().join("products"))?;
        std::fs::create_dir_all(w.body_dir().join("edits"))?;
        w.save_body(body)?;
        w.save_manifest()?;
        Ok(w)
    }

    /// The world at `root`; at a universe directory, its focus body's.
    ///
    /// # Errors
    /// If there is no world at `root` or its manifest is invalid.
    pub fn open(root: &Path) -> Result<Self> {
        if !root.join("world.toml").exists() && UniverseDir::is_at(root) {
            let u = UniverseDir::open(root)?;
            let focus = u.doc.focus.ok_or_else(|| StoreError::Invalid(format!("the universe at {} has no body yet: geoforge new {} --at <address>", root.display(), root.display())))?;
            return Self::open(&u.body_root(&focus));
        }
        let text = std::fs::read_to_string(root.join("world.toml")).map_err(|e| StoreError::Invalid(format!("no world at {}: {e}", root.display())))?;
        let manifest: Manifest = toml::from_str(&text)?;
        if manifest.schema_version != geoforge_core::SCHEMA_VERSION {
            return Err(StoreError::Invalid(format!("world schema {} but this build reads {}", manifest.schema_version, geoforge_core::SCHEMA_VERSION)));
        }
        let universe = manifest.address.and_then(|_| universe_above(root));
        Ok(Self { root: root.to_path_buf(), manifest, universe })
    }

    /// # Errors
    /// On IO errors.
    pub fn save_manifest(&self) -> Result<()> {
        std::fs::write(self.root.join("world.toml"), toml::to_string_pretty(&self.manifest)?)?;
        Ok(())
    }

    #[must_use]
    pub fn body_dir(&self) -> PathBuf {
        self.root.join("bodies").join(slug(&self.manifest.focus_name))
    }
    #[must_use]
    pub fn products_dir(&self) -> PathBuf {
        self.body_dir().join("products")
    }
    #[must_use]
    pub fn edits_dir(&self) -> PathBuf {
        self.body_dir().join("edits")
    }
    #[must_use]
    pub fn baseline_dir(&self) -> PathBuf {
        self.root.join("cache").join("baselines")
    }

    /// # Errors
    /// On IO or parse errors.
    pub fn load_body(&self) -> Result<BodyContext> {
        Ok(serde_json::from_str(&std::fs::read_to_string(self.body_dir().join("body.json"))?)?)
    }
    /// # Errors
    /// On IO errors.
    pub fn save_body(&self, b: &BodyContext) -> Result<()> {
        std::fs::create_dir_all(self.body_dir())?;
        std::fs::write(self.body_dir().join("body.json"), serde_json::to_string_pretty(b)?)?;
        Ok(())
    }
    /// The world's galaxy documents: the galaxy model (`galaxy.json`, if edited or the world
    /// has its own) and the edits over its stars (`stars.json`). A body of a universe shares
    /// its galaxy's with every other body there.
    #[must_use]
    pub fn galaxy_dir(&self) -> PathBuf {
        match (&self.universe, self.manifest.address) {
            (Some(u), Some(a)) => galaxy_dir_in(u, &a),
            _ => self.root.join("galaxy"),
        }
    }

    /// The galaxy model the world keeps, if any (else its galaxy's is generated).
    ///
    /// # Errors
    /// On IO or parse errors.
    pub fn load_galaxy_model<T: serde::de::DeserializeOwned>(&self) -> Result<Option<T>> {
        read_json(&self.galaxy_dir().join(GALAXY_MODEL))
    }

    /// # Errors
    /// On IO errors.
    pub fn save_galaxy_model<T: Serialize>(&self, model: &T) -> Result<()> {
        write_json(&self.galaxy_dir().join(GALAXY_MODEL), model)
    }

    /// A JSON document at `rel` under the world, if it exists.
    ///
    /// # Errors
    /// On IO or parse errors (a missing file is `Ok(None)`).
    pub fn load_doc<T: serde::de::DeserializeOwned>(&self, rel: &str) -> Result<Option<T>> {
        read_json(&self.root.join(rel))
    }

    /// Write `doc` as JSON at `rel` under the world.
    ///
    /// # Errors
    /// On IO errors.
    pub fn save_doc<T: Serialize>(&self, rel: &str, doc: &T) -> Result<()> {
        write_json(&self.root.join(rel), doc)
    }

    /// The edits over the galaxy's generated stars (none if the world has made none).
    ///
    /// # Errors
    /// On IO or parse errors.
    pub fn load_star_edits(&self) -> Result<geoforge_core::galactic::StarEdits> {
        Ok(read_json(&self.galaxy_dir().join(STAR_EDITS))?.unwrap_or_default())
    }

    /// # Errors
    /// On IO errors.
    pub fn save_star_edits(&self, edits: &geoforge_core::galactic::StarEdits) -> Result<()> {
        write_json(&self.galaxy_dir().join(STAR_EDITS), edits)
    }

    /// The constellations the focus's sky is drawn with.
    ///
    /// # Errors
    /// On IO or parse errors.
    pub fn load_constellations(&self) -> Result<Vec<geoforge_core::galactic::Constellation>> {
        Ok(self.load_doc(CONSTELLATIONS)?.unwrap_or_default())
    }

    /// # Errors
    /// On IO errors.
    pub fn save_constellations(&self, c: &[geoforge_core::galactic::Constellation]) -> Result<()> {
        self.save_doc(CONSTELLATIONS, &c)
    }

    /// The focus world's own calendar (`calendar.json` beside its `body.json`), if it has one.
    ///
    /// # Errors
    /// On IO or parse errors (a missing file is `Ok(None)`).
    pub fn load_calendar(&self) -> Result<Option<geoforge_core::time::CalendarDef>> {
        let p = self.body_dir().join("calendar.json");
        if !p.exists() {
            return Ok(None);
        }
        Ok(Some(serde_json::from_str(&std::fs::read_to_string(p)?)?))
    }
    /// # Errors
    /// On parse errors (a missing file is `Ok(None)`).
    pub fn load_system(&self) -> Result<Option<SolarSystem>> {
        let p = self.root.join("system.json");
        if !p.exists() {
            return Ok(None);
        }
        Ok(Some(serde_json::from_str(&std::fs::read_to_string(p)?)?))
    }

    /// Is a product stored?
    #[must_use]
    pub fn has_product(&self, p: ProductId) -> bool {
        self.products_dir().join(format!("{}.json", p.name())).exists()
    }

    /// # Errors
    /// On IO or format errors.
    pub fn save_product(&self, data: &ProductData, provenance: Hash) -> Result<Hash> {
        codec::write_canonical(&self.products_dir(), data, provenance)
    }

    /// Every stored product (no generation).
    ///
    /// # Errors
    /// On unreadable products.
    pub fn load_products(&self) -> Result<geoforge_core::product::Products> {
        let mut out = geoforge_core::product::Products::new();
        for p in ProductId::ALL {
            if let Some(d) = self.load_product(p)? {
                out.insert(d);
            }
        }
        Ok(out)
    }

    /// The body, its grid and every stored product, loaded together.
    ///
    /// # Errors
    /// On unreadable files or an unbuildable grid.
    pub fn snapshot(&self) -> Result<Snapshot> {
        Ok(Snapshot { body: self.load_body()?, grid: self.manifest.grid.build().map_err(|e| StoreError::Invalid(e.to_string()))?, products: self.load_products()? })
    }

    /// # Errors
    /// On IO or format errors.
    pub fn load_product(&self, p: ProductId) -> Result<Option<ProductData>> {
        if !self.has_product(p) {
            return Ok(None);
        }
        codec::read_canonical(&self.products_dir(), p, self.manifest.grid).map(Some)
    }

    /// The provenance recorded inside a stored product.
    #[must_use]
    pub fn stored_provenance(&self, p: ProductId) -> Option<Hash> {
        codec::read_provenance(&self.products_dir(), p)
    }

    /// # Errors
    /// On IO errors.
    pub fn remove_product(&self, p: ProductId) -> Result<()> {
        for ext in ["json", "bin"] {
            let f = self.products_dir().join(format!("{}.{ext}", p.name()));
            if f.exists() {
                std::fs::remove_file(f)?;
            }
        }
        Ok(())
    }
}
