//! Runs the surface pipeline over a world directory.
//!
//! - **Plan**: for each layer, which implementation fills it on this body (see [`registry`])
//!   or why none applies, and for each product its source (generate / patched / authored)
//!   and state (fresh / stale / missing).
//! - **Provenance**: a product's hash covers the implementation and its version, the root
//!   seed, the body, the grid, the implementation's parameters, its products' sources and
//!   edit files, and the provenance of every input it reads. A stored product whose recorded provenance matches
//!   is reused; anything else is rebuilt. Change one thing and exactly what depends on it
//!   goes stale.
//! - **Edits**: authored products replace the generator's output; patched products are
//!   generated, patched, and then fed back so the rest of their layer is consistent with
//!   the edit. Downstream layers always see the edited version.

pub mod place;
pub mod registry;
pub mod sky;
pub mod stats;
pub mod weather;
pub mod views;

pub use registry::Implementation;
pub use views::View;

use geoforge_core::grid::{AnyGrid, Cell, Grid, GridSpec};
use geoforge_core::key::Key;
use geoforge_core::product::{Hash, LayerId, ProductId, Products, RunContext, Source};
use geoforge_core::{BodyContext, Seed};
use geoforge_store::edit::{self, BodyRef, Patch};
use geoforge_store::{codec, hash_parts, Manifest, World};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use std::time::Instant;

#[derive(Debug, thiserror::Error)]
pub enum PipelineError {
    #[error("{0}")]
    Store(#[from] geoforge_store::StoreError),
    #[error("{0}")]
    Core(#[from] geoforge_core::Error),
    #[error("json: {0}")]
    Json(#[from] serde_json::Error),
    #[error("io: {0}")]
    Io(#[from] std::io::Error),
    #[error("{0}")]
    Invalid(String),
}

pub type Result<T> = std::result::Result<T, PipelineError>;

/// Version of this build, recorded in new worlds.
pub const CODE_VERSION: &str = env!("CARGO_PKG_VERSION");

// ---------------------------------------------------------------- worlds

/// Options for a new world.
#[derive(Debug, Clone)]
pub struct NewWorld {
    pub seed: Seed,
    pub grid: GridSpec,
    /// How to make the solar system: drawn from a galaxy (the default, the Milky Way), or the
    /// habitable recipe, a Sun-like star with an Earth-like world.
    pub recipe: geoforge_cosmos::Recipe,
    /// Use this body instead of generating a solar system.
    pub body: Option<BodyContext>,
    /// Focus on this body of the generated system instead of its habitable world: a body's
    /// name, or `"<belt name> #<n>"` for the n-th largest member of a belt.
    pub focus: Option<String>,
}

/// Create a world: a generated solar system and its focus body, or a hand-written body.
///
/// # Errors
/// If the directory already holds a world, the grid is invalid, the focus name is unknown, or
/// a system drawn from a galaxy has nothing solid in it.
pub fn create_world(root: &Path, opts: &NewWorld) -> Result<World> {
    opts.grid.build()?;
    let (system, body) = match &opts.body {
        Some(b) => (None, b.clone()),
        None => {
            let g = opts.recipe.generate(opts.seed).map_err(|b| PipelineError::Invalid(format!("{b}; try another seed (`census` lists some), or the habitable recipe")))?;
            let body = match &opts.focus {
                Some(name) => geoforge_cosmos::find_body(&g.system, name).ok_or_else(|| PipelineError::Invalid(format!("no body named '{name}' in the system")))?,
                None => geoforge_cosmos::any_body(&g.system, g.focus).expect("the focus is a body"),
            };
            (Some(g.system), body)
        }
    };
    // Defaults for the implementations this body uses; others fall back to defaults if the
    // body later changes.
    let params = LayerId::ALL
        .into_iter()
        .filter_map(|l| registry::select(l, &body).ok())
        .map(|i| (i.name().to_string(), i.default_params()))
        .collect();
    let manifest = Manifest {
        schema_version: geoforge_core::SCHEMA_VERSION,
        code_version: CODE_VERSION.into(),
        root_seed: opts.seed,
        focus: body.key,
        focus_name: body.name.clone(),
        grid: opts.grid,
        params,
        sources: BTreeMap::new(),
        provenance: BTreeMap::new(),
        not_applicable: BTreeMap::new(),
    };
    Ok(World::create(root, manifest, system.as_ref(), &body)?)
}

// ---------------------------------------------------------------- plan

/// State of one product in a plan.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum State {
    Fresh,
    Stale,
    Missing,
    NotApplicable,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProductPlan {
    pub product: ProductId,
    pub source: String,
    pub state: State,
    pub provenance: Hash,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LayerPlan {
    pub layer: LayerId,
    /// The generator that fills the layer. `None` when none applies, or when every product
    /// of the layer is authored and no generator is needed.
    pub implementation: Option<Implementation>,
    /// `None` when the layer applies; otherwise the reason it does not.
    pub not_applicable: Option<String>,
    pub products: Vec<ProductPlan>,
}

impl LayerPlan {
    #[must_use]
    pub fn needs_run(&self) -> bool {
        self.not_applicable.is_none() && self.products.iter().any(|p| p.state != State::Fresh)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Plan {
    pub layers: Vec<LayerPlan>,
}

fn json_bytes<T: Serialize>(v: &T) -> Result<Vec<u8>> {
    Ok(serde_json::to_vec(v)?)
}

/// Hash of a product's edit file (authored data or patch), if it has one.
fn edit_hash(world: &World, p: ProductId, src: &Source) -> Result<Option<Hash>> {
    Ok(match src {
        Source::Generate { .. } => None,
        Source::Authored { .. } => Some(
            codec::read_content_hash(&world.edits_dir(), p)
                .ok_or_else(|| PipelineError::Invalid(format!("{p} is marked authored but its edit file is missing or unreadable")))?,
        ),
        Source::Patched { patch, .. } => {
            let bytes = std::fs::read(world.body_dir().join(patch))
                .map_err(|e| PipelineError::Invalid(format!("{p} is marked patched but {patch} is unreadable: {e}")))?;
            Some(hash_parts(&[b"patch-file", &bytes]))
        }
    })
}

/// Parameters of an implementation, as the world sets them (defaults filled in).
///
/// # Errors
/// If the stored parameters do not parse.
pub fn params_of(world: &World, i: Implementation) -> Result<serde_json::Value> {
    i.resolve_params(world.manifest.params.get(i.name()))
}

/// Work out what is fresh, stale or missing.
///
/// # Errors
/// On unreadable world files.
pub fn plan(world: &World) -> Result<Plan> {
    let body = world.load_body()?;
    let body_hash = hash_parts(&[b"body", &json_bytes(&body)?]);
    let grid_bytes = json_bytes(&world.manifest.grid)?;
    let mut available: BTreeMap<ProductId, Hash> = BTreeMap::new();
    let mut out = Vec::new();
    for l in LayerId::ALL {
        let sources: Vec<(ProductId, Source)> = l.products().into_iter().map(|p| (p, world.manifest.source(p))).collect();
        let all_user = sources.iter().all(|(_, s)| matches!(s, Source::Authored { .. }));
        let selected = registry::select(l, &body);
        let reason = match (&selected, all_user) {
            (_, true) => None,
            (Err(pre), false) => Some(pre.0.clone()),
            (Ok(i), false) => i.needs().iter().find(|n| !available.contains_key(n)).map(|n| format!("needs {n}, which is not available")),
        };
        let implementation = selected.ok();
        let mut lp = LayerPlan { layer: l, implementation, not_applicable: reason.clone(), products: Vec::new() };
        if reason.is_some() {
            for (p, s) in sources {
                lp.products.push(ProductPlan { product: p, source: s.label().into(), state: State::NotApplicable, provenance: Hash([0; 32]) });
            }
            out.push(lp);
            continue;
        }
        // Layer provenance. A fully authored layer with no applicable generator hashes only
        // its edits.
        let (name, version, params) = match implementation {
            Some(i) => (i.name(), i.version(), params_of(world, i)?),
            None => ("authored", 0, serde_json::Value::Null),
        };
        let mut parts: Vec<Vec<u8>> = vec![
            b"layer".to_vec(),
            l.name().as_bytes().to_vec(),
            name.as_bytes().to_vec(),
            version.to_le_bytes().to_vec(),
            world.manifest.root_seed.0.to_le_bytes().to_vec(),
            body_hash.0.to_vec(),
            grid_bytes.clone(),
            json_bytes(&params)?,
        ];
        for (p, s) in &sources {
            parts.push(p.name().as_bytes().to_vec());
            parts.push(json_bytes(s)?);
            if let Some(h) = edit_hash(world, *p, s)? {
                parts.push(h.0.to_vec());
            }
        }
        let mut inputs: Vec<ProductId> = implementation.map_or(Vec::new(), |i| i.needs().iter().chain(i.wants()).copied().collect());
        inputs.sort();
        inputs.dedup();
        for i in inputs {
            parts.push(i.name().as_bytes().to_vec());
            match available.get(&i) {
                Some(h) => parts.push(h.0.to_vec()),
                None => parts.push(b"absent".to_vec()),
            }
        }
        let refs: Vec<&[u8]> = parts.iter().map(Vec::as_slice).collect();
        let layer_prov = hash_parts(&refs);
        for (p, s) in sources {
            let prov = hash_parts(&[&layer_prov.0, p.name().as_bytes()]);
            let state = match world.stored_provenance(p) {
                None => State::Missing,
                Some(h) if h == prov => State::Fresh,
                Some(_) => State::Stale,
            };
            available.insert(p, prov);
            lp.products.push(ProductPlan { product: p, source: s.label().into(), state, provenance: prov });
        }
        out.push(lp);
    }
    Ok(Plan { layers: out })
}

// ---------------------------------------------------------------- generate

#[derive(Debug, Clone, Default)]
pub struct GenerateOptions {
    /// Rebuild everything, ignoring stored products.
    pub force: bool,
    /// Stop after this layer.
    pub until: Option<LayerId>,
    /// Rebuild this layer and everything after it, even if fresh.
    pub from: Option<LayerId>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "action", rename_all = "snake_case")]
pub enum LayerAction {
    Ran { seconds: f64 },
    Reused,
    NotApplicable { reason: String },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LayerReport {
    pub layer: LayerId,
    pub action: LayerAction,
    pub notes: Vec<String>,
}

/// Result of applying an edit during the last run, kept in `cache/status.json`.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct EditStatus {
    pub rebased: bool,
    pub orphaned: usize,
}

fn status_path(world: &World) -> PathBuf {
    world.root.join("cache").join("status.json")
}

/// Rebased/orphaned edits from the last run.
#[must_use]
pub fn edit_status(world: &World) -> BTreeMap<String, EditStatus> {
    std::fs::read_to_string(status_path(world)).ok().and_then(|s| serde_json::from_str(&s).ok()).unwrap_or_default()
}

fn load_patch(world: &World, rel: &str) -> Result<Patch> {
    Ok(serde_json::from_str(&std::fs::read_to_string(world.body_dir().join(rel))?)?)
}

/// Fill rows for partition keys that an edit introduced without attributes.
pub fn filler<'a>(root_seed: Seed, body: &'a BodyContext, grid: &'a AnyGrid, ctx: &'a Products) -> impl Fn(ProductId, Key, &[Cell]) -> Option<serde_json::Value> + 'a {
    move |p, key, cells| match p {
        ProductId::Plates => {
            let pts: Vec<_> = cells.iter().map(|&c| grid.center(c)).collect();
            let seed_point = geoforge_core::geom::UnitVec3::centroid(pts.iter().copied())?;
            let s = root_seed.for_product(ProductId::Plates, key);
            let dir = geoforge_core::seed::range(s.child(b"dir"), 0.0, 360.0);
            let speed = geoforge_core::seed::range(s.child(b"speed"), 1.0, 10.0);
            let (pole, rate) = geoforge_surface::tectonics::motion::euler_pole_for(seed_point, dir, speed, body.exterior.radius_km);
            serde_json::to_value(geoforge_core::surface::Plate { name: format!("Plate {key}"), euler_pole: pole, rate_deg_per_myr: rate, seed_point }).ok()
        }
        ProductId::Lakes => {
            let e = ctx.elevation()?;
            let top = cells.iter().map(|&c| e.values[c as usize]).fold(f32::NEG_INFINITY, f32::max);
            serde_json::to_value(geoforge_core::surface::Lake { surface_elevation_m: f64::from(top), endorheic: false }).ok()
        }
        _ => None,
    }
}

/// Run the pipeline. Fresh products are loaded; everything else is (re)built and stored.
///
/// # Errors
/// On generator errors or unreadable world files.
pub fn generate(world: &mut World, opts: &GenerateOptions, mut log: impl FnMut(&LayerReport)) -> Result<(Products, Vec<LayerReport>)> {
    let body = world.load_body()?;
    let grid = world.manifest.grid.build()?;
    let plan = plan(world)?;
    let mut products = Products::new();
    let mut reports = Vec::new();
    let mut status: BTreeMap<String, EditStatus> = edit_status(world);
    world.manifest.not_applicable.clear();
    for lp in &plan.layers {
        let l = lp.layer;
        let mut notes = Vec::new();
        if let Some(reason) = &lp.not_applicable {
            world.manifest.not_applicable.insert(l.name().into(), reason.clone());
            for p in l.products() {
                world.remove_product(p)?;
                world.manifest.provenance.remove(p.name());
            }
            let r = LayerReport { layer: l, action: LayerAction::NotApplicable { reason: reason.clone() }, notes };
            log(&r);
            reports.push(r);
            continue;
        }
        let rebuild = opts.force || opts.from.is_some_and(|f| l >= f) || lp.needs_run();
        if !rebuild {
            for pp in &lp.products {
                if let Some(d) = world.load_product(pp.product)? {
                    products.insert(d);
                }
            }
            let r = LayerReport { layer: l, action: LayerAction::Reused, notes };
            log(&r);
            reports.push(r);
            if opts.until == Some(l) {
                break;
            }
            continue;
        }
        let t0 = Instant::now();
        let run = |inputs: &Products| -> Result<Products> {
            let i = lp.implementation.ok_or_else(|| PipelineError::Invalid(format!("no generator fills {}", l.name())))?;
            i.run(&RunContext { root_seed: world.manifest.root_seed, body: &body, grid: &grid }, inputs, &params_of(world, i)?)
        };
        let mut inputs = products.clone();
        let mut patched = Vec::new();
        for pp in &lp.products {
            match world.manifest.source(pp.product) {
                Source::Authored { .. } => {
                    let d = codec::read_canonical(&world.edits_dir(), pp.product, world.manifest.grid)?;
                    inputs.insert(d);
                }
                Source::Patched { patch, .. } => patched.push((pp.product, patch)),
                Source::Generate { .. } => {}
            }
        }
        let all_user = lp.products.iter().all(|pp| inputs.contains(pp.product));
        let mut out = if all_user { Products::new() } else { run(&inputs)? };
        if !patched.is_empty() {
            for (p, rel) in &patched {
                let mut d = out.remove(*p).ok_or_else(|| PipelineError::Invalid(format!("{p}: the layer did not produce a base to patch")))?;
                // Store the generated base so later imports can be diffed against it.
                codec::write_canonical(&world.root.join("cache").join("bases"), &d, Hash([0; 32]))?;
                let patch = load_patch(world, rel)?;
                let fill = filler(world.manifest.root_seed, &body, &grid, &products);
                let a = edit::apply(&mut d, &patch, &fill)?;
                codec::validate(&d)?;
                if a.rebased {
                    notes.push(format!("{p}: base changed since the edit was made; edits re-applied on top"));
                }
                if a.orphaned > 0 {
                    notes.push(format!("{p}: {} edited cells could not be applied", a.orphaned));
                }
                status.insert(p.name().into(), EditStatus { rebased: a.rebased, orphaned: a.orphaned });
                inputs.insert(d);
            }
            // Rebuild the rest of the layer so it agrees with the edited products.
            let siblings_generated = lp.products.iter().any(|pp| !inputs.contains(pp.product));
            out = if siblings_generated { run(&inputs)? } else { Products::new() };
        }
        for pp in &lp.products {
            if !out.contains(pp.product) {
                if let Some(d) = inputs.remove(pp.product) {
                    out.insert(d);
                }
            }
            if !matches!(world.manifest.source(pp.product), Source::Patched { .. }) {
                status.remove(pp.product.name());
            }
        }
        for pp in &lp.products {
            let d = out.get(pp.product).ok_or_else(|| PipelineError::Invalid(format!("layer {} did not produce {}", l.name(), pp.product)))?;
            world.save_product(d, pp.provenance)?;
            world.manifest.provenance.insert(pp.product.name().into(), pp.provenance);
        }
        products.extend(out);
        let r = LayerReport { layer: l, action: LayerAction::Ran { seconds: t0.elapsed().as_secs_f64() }, notes };
        log(&r);
        reports.push(r);
        if opts.until == Some(l) {
            break;
        }
    }
    world.save_manifest()?;
    std::fs::create_dir_all(world.root.join("cache"))?;
    std::fs::write(status_path(world), serde_json::to_string_pretty(&status)?)?;
    Ok((products, reports))
}

/// Load every stored product (no generation).
///
/// # Errors
/// On unreadable products.
pub fn load_all(world: &World) -> Result<Products> {
    let mut out = Products::new();
    for p in ProductId::ALL {
        if let Some(d) = world.load_product(p)? {
            out.insert(d);
        }
    }
    Ok(out)
}

// ---------------------------------------------------------------- edit views

/// Default folder for edit views.
#[must_use]
pub fn edit_folder(world: &World) -> PathBuf {
    world.root.join("edit")
}

/// Export edit views of stored products. Returns the files written per product.
///
/// # Errors
/// If a product is missing (run `generate` first) or on IO errors.
pub fn export(world: &World, ps: &[ProductId], to: &Path) -> Result<Vec<(ProductId, Vec<String>)>> {
    let body = world.load_body()?;
    let bref = BodyRef { name: body.name.clone(), radius_km: body.exterior.radius_km };
    let mut out = Vec::new();
    for &p in ps {
        let d = world.load_product(p)?.ok_or_else(|| PipelineError::Invalid(format!("{p} is not generated yet; run generate first")))?;
        let files = edit::export_view(to, &d, &bref)?;
        codec::write_canonical(&world.baseline_dir(), &d, Hash([0; 32]))?;
        out.push((p, files));
    }
    Ok(out)
}

#[derive(Debug, Clone)]
pub struct ImportReport {
    pub product: ProductId,
    pub mode: &'static str,
    pub changed_cells: usize,
    pub snapped: usize,
    pub rows_changed: usize,
}

/// Import an edited view as a patch (default) or as authored data.
///
/// # Errors
/// If the view is invalid, or a patch is requested without a baseline.
pub fn import(world: &mut World, p: ProductId, from: &Path, authored: bool) -> Result<ImportReport> {
    let grid = world.manifest.grid;
    let view = edit::read_view(from, p, grid)?;
    let body = world.load_body()?;
    let g = grid.build()?;
    let ctx = load_all(world)?;
    let fill = filler(world.manifest.root_seed, &body, &g, &ctx);
    let edits = world.edits_dir();
    std::fs::create_dir_all(&edits)?;
    let patch_file = edits.join(format!("{}.patch.json", p.name()));
    let use_patch = !authored && matches!(p.shape(), geoforge_core::product::Shape::Field | geoforge_core::product::Shape::Partition);
    if !use_patch {
        let d = edit::view_to_product(&view, grid, &fill)?;
        let h = codec::write_canonical(&edits, &d, Hash([0; 32]))?;
        let _ = std::fs::remove_file(&patch_file);
        world.manifest.sources.insert(p.name().into(), Source::Authored { file: format!("edits/{}.json", p.name()) });
        world.save_manifest()?;
        let _ = h;
        return Ok(ImportReport { product: p, mode: "authored", changed_cells: grid.cell_count(), snapped: view.snapped, rows_changed: view.table.len() });
    }
    let baseline = codec::read_canonical(&world.baseline_dir(), p, grid)
        .map_err(|_| PipelineError::Invalid(format!("no export baseline for {p}; run `export {p}` before editing, or import with --authored")))?;
    let fresh = edit::diff(&baseline, &view)?;
    let (patch, params) = match world.manifest.source(p) {
        Source::Patched { patch: rel, params, .. } => (edit::merge(&load_patch(world, &rel)?, &fresh)?, params),
        Source::Generate { params } => (fresh.clone(), params),
        Source::Authored { .. } => {
            // Editing an authored product: fold the edit into the authored data.
            let mut d = codec::read_canonical(&edits, p, grid)?;
            edit::apply(&mut d, &fresh, &fill)?;
            codec::write_canonical(&edits, &d, Hash([0; 32]))?;
            return Ok(ImportReport { product: p, mode: "authored (updated)", changed_cells: fresh.cell_count(), snapped: view.snapped, rows_changed: fresh.rows.len() });
        }
    };
    std::fs::write(&patch_file, serde_json::to_string(&patch)?)?;
    world.manifest.sources.insert(p.name().into(), Source::Patched { params, base: patch.base, patch: format!("edits/{}.patch.json", p.name()) });
    world.save_manifest()?;
    Ok(ImportReport { product: p, mode: "patch", changed_cells: fresh.cell_count(), snapped: view.snapped, rows_changed: fresh.rows.len() })
}

/// Drop a product's edits; it goes back to being generated.
///
/// # Errors
/// On IO errors.
pub fn reset(world: &mut World, p: ProductId) -> Result<()> {
    for f in [format!("{}.patch.json", p.name()), format!("{}.json", p.name()), format!("{}.bin", p.name())] {
        let path = world.edits_dir().join(f);
        if path.exists() {
            std::fs::remove_file(path)?;
        }
    }
    world.manifest.sources.remove(p.name());
    world.save_manifest()?;
    Ok(())
}

/// Products whose layer is in `layers` (all if empty).
#[must_use]
pub fn products_of(layers: &[LayerId]) -> Vec<ProductId> {
    let set: BTreeSet<LayerId> = layers.iter().copied().collect();
    ProductId::ALL.into_iter().filter(|p| set.is_empty() || set.contains(&p.layer())).collect()
}

/// Write preview renders of the stored raster `products` and derived `views` into `to`.
///
/// # Errors
/// On IO errors.
pub fn render_all(world: &World, stored: &Products, to: &Path, products: &[ProductId], views: &[View]) -> Result<Vec<PathBuf>> {
    std::fs::create_dir_all(to)?;
    let body = world.load_body()?;
    let grid = world.manifest.grid.build()?;
    let (radius, sea) = (body.exterior.radius_km, body.interior.has_surface_liquid_water());
    let mut out = Vec::new();
    for &p in products {
        if let Some(img) = geoforge_store::render::render(p, stored, &grid, radius, sea, body.exterior.class)? {
            let f = to.join(format!("{}.png", p.name()));
            img.save(&f).map_err(geoforge_store::StoreError::from)?;
            out.push(f);
        }
    }
    let source = views::Source::new(stored, &grid, &body);
    let sea_cells = source.sea();
    for &v in views {
        let relief = Some((stored, radius));
        let img = match views::compute(v, &source) {
            Some(views::Map::Classes(c)) => geoforge_store::render::render_classes(&c.colors(), &grid, relief)?,
            Some(views::Map::Chances(c)) => geoforge_store::render::render_chances(&c.cells, &sea_cells, &grid, relief)?,
            Some(views::Map::Quantity(q)) => geoforge_store::render::render_quantity(&q.cells, &sea_cells, q.contour, 6.0 * q.contour, &grid, relief)?,
            None => continue,
        };
        let f = to.join(format!("{}.png", v.name()));
        img.save(&f).map_err(geoforge_store::StoreError::from)?;
        out.push(f);
    }
    Ok(out)
}

/// Export derived views with world files: a class map as an exact-colour PNG and a legend
/// (`<view>.legend.json`: code, name, colour and land share of each class present); a chance
/// map as an 8-bit grayscale PNG and a note of what it means (`<view>.json`). Returns the
/// files written per view; views whose products are missing are skipped.
///
/// # Errors
/// On IO errors.
pub fn export_views(world: &World, views: &[View], to: &Path) -> Result<Vec<(View, Vec<String>)>> {
    std::fs::create_dir_all(to)?;
    let body = world.load_body()?;
    let grid = world.manifest.grid.build()?;
    let stored = load_all(world)?;
    let bref = BodyRef { name: body.name.clone(), radius_km: body.exterior.radius_km };
    let source = views::Source::new(&stored, &grid, &body);
    let mut out = Vec::new();
    for &v in views {
        let name = v.name();
        let png = to.join(format!("{name}.png"));
        let sidecar = match views::compute(v, &source) {
            Some(views::Map::Classes(c)) => {
                geoforge_store::render::render_classes(&c.colors(), &grid, None)?.save(&png).map_err(geoforge_store::StoreError::from)?;
                std::fs::write(to.join(format!("{name}.legend.json")), serde_json::to_string_pretty(&c.legend)?)?;
                "legend.json"
            }
            Some(views::Map::Chances(c)) => {
                geoforge_store::render::gray_chances(&c.cells, &grid)?.save(&png).map_err(geoforge_store::StoreError::from)?;
                let note = serde_json::json!({ "view": name, "value": "0–255 for 0–1", "one_means": c.meaning });
                std::fs::write(to.join(format!("{name}.json")), serde_json::to_string_pretty(&note)?)?;
                "json"
            }
            Some(views::Map::Quantity(q)) => {
                let (img, min, max) = geoforge_store::render::gray_quantity(&q.cells, &grid)?;
                img.save(&png).map_err(geoforge_store::StoreError::from)?;
                let note = serde_json::json!({ "view": name, "units": q.units, "min": min, "max": max, "value": "min + v/65535 × (max − min)" });
                std::fs::write(to.join(format!("{name}.json")), serde_json::to_string_pretty(&note)?)?;
                "json"
            }
            None => continue,
        };
        edit::write_georef(to, name, world.manifest.grid, &bref)?;
        out.push((v, ["png", sidecar, "pgw", "prj"].map(|e| format!("{name}.{e}")).to_vec()));
    }
    Ok(out)
}
