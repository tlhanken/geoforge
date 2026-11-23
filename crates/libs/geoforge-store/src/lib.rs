//! The world directory. Your sources and edits are kept; generated products are a cache
//! that can always be rebuilt.
//!
//! ```text
//! world/
//!   world.toml                       manifest: seed, grid, focus, sources, provenance
//!   system.json                      the solar system (absent for a hand-written body)
//!   bodies/<body>/
//!     body.json                      BodyContext (the handoff to the surface pipeline)
//!     products/<product>.json|.bin   canonical stored products
//!     edits/<product>.*              authored data and patches (never auto-overwritten)
//!   cache/baselines/                 last edit-view export of each product (for diffs)
//! ```

pub mod cache;
pub mod codec;
pub mod edit;
pub mod polar;
pub mod raster;
pub mod render;
pub mod report;

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

/// Where a world keeps its galaxy model, its star edits and its constellations.
pub const GALAXY_MODEL: &str = "galaxy/galaxy.json";
const STAR_EDITS: &str = "galaxy/stars.json";
const CONSTELLATIONS: &str = "constellations.json";

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
        let w = Self { root: root.to_path_buf(), manifest };
        if let Some(s) = system {
            std::fs::write(w.root.join("system.json"), serde_json::to_string_pretty(s)?)?;
        }
        std::fs::create_dir_all(w.body_dir().join("products"))?;
        std::fs::create_dir_all(w.body_dir().join("edits"))?;
        w.save_body(body)?;
        w.save_manifest()?;
        Ok(w)
    }

    /// # Errors
    /// If there is no world at `root` or its manifest is invalid.
    pub fn open(root: &Path) -> Result<Self> {
        let text = std::fs::read_to_string(root.join("world.toml")).map_err(|e| StoreError::Invalid(format!("no world at {}: {e}", root.display())))?;
        let manifest: Manifest = toml::from_str(&text)?;
        if manifest.schema_version != geoforge_core::SCHEMA_VERSION {
            return Err(StoreError::Invalid(format!("world schema {} but this build reads {}", manifest.schema_version, geoforge_core::SCHEMA_VERSION)));
        }
        Ok(Self { root: root.to_path_buf(), manifest })
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
    /// The world's galaxy documents: the galaxy model (`galaxy.json`, if the world has its
    /// own) and the edits over its stars (`stars.json`).
    #[must_use]
    pub fn galaxy_dir(&self) -> PathBuf {
        self.root.join("galaxy")
    }

    /// A JSON document at `rel` under the world, if it exists.
    ///
    /// # Errors
    /// On IO or parse errors (a missing file is `Ok(None)`).
    pub fn load_doc<T: serde::de::DeserializeOwned>(&self, rel: &str) -> Result<Option<T>> {
        let p = self.root.join(rel);
        if !p.exists() {
            return Ok(None);
        }
        Ok(Some(serde_json::from_str(&std::fs::read_to_string(p)?)?))
    }

    /// Write `doc` as JSON at `rel` under the world.
    ///
    /// # Errors
    /// On IO errors.
    pub fn save_doc<T: Serialize>(&self, rel: &str, doc: &T) -> Result<()> {
        let p = self.root.join(rel);
        if let Some(dir) = p.parent() {
            std::fs::create_dir_all(dir)?;
        }
        std::fs::write(p, serde_json::to_string_pretty(doc)?)?;
        Ok(())
    }

    /// The edits over the galaxy's generated stars (none if the world has made none).
    ///
    /// # Errors
    /// On IO or parse errors.
    pub fn load_star_edits(&self) -> Result<geoforge_core::galactic::StarEdits> {
        Ok(self.load_doc(STAR_EDITS)?.unwrap_or_default())
    }

    /// # Errors
    /// On IO errors.
    pub fn save_star_edits(&self, edits: &geoforge_core::galactic::StarEdits) -> Result<()> {
        self.save_doc(STAR_EDITS, edits)
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
