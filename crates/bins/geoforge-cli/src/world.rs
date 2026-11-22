//! Making a world and running its layers: `new`, `generate`, `plan`, `status`, `stats`,
//! `seed`, and the names commands take.

use crate::system::show;
use crate::R;
use geoforge_core::grid::GridSpec;
use geoforge_core::product::{LayerId, ProductId};
use geoforge_core::{BodyContext, Seed};
use geoforge_cosmos::{GalaxyPreset, Recipe};
use geoforge_pipeline::{self as pipe, GenerateOptions, LayerAction, NewWorld, State, View};
use geoforge_store::World;
use std::path::{Path, PathBuf};

#[allow(clippy::too_many_arguments)]
pub fn new_world(dir: &Path, seed: u64, galaxy: GalaxyPreset, habitable: bool, width: u32, body: Option<PathBuf>, archetype: Option<geoforge_core::Archetype>, focus: Option<String>) -> R<()> {
    if seed > i64::MAX as u64 {
        return Err("seed must be at most 9223372036854775807".into());
    }
    let mut body = body.as_deref().map(read_body).transpose()?;
    if let Some(a) = archetype {
        body = Some(a.context(geoforge_core::Key(1), a.name()));
    }
    let recipe = if habitable { Recipe::Habitable } else { Recipe::Galaxy(galaxy.galaxy()) };
    let w = pipe::create_world(dir, &NewWorld { seed: Seed(seed), recipe, grid: GridSpec::Equirect { width, height: width / 2 }, body, focus })?;
    show(&w)?;
    println!("\ncreated {} — next: geoforge generate {}", dir.display(), dir.display());
    Ok(())
}

pub fn generate(dir: &Path, force: bool, until: Option<String>, from: Option<String>, render: bool) -> R<()> {
    let mut w = World::open(dir)?;
    let until = until.map(|u| LayerId::from_name(&u).ok_or(format!("unknown layer `{u}`"))).transpose()?;
    let from = from.map(|u| LayerId::from_name(&u).ok_or(format!("unknown layer `{u}`"))).transpose()?;
    let (products, _) = pipe::generate(&mut w, &GenerateOptions { force, until, from }, |r| {
        let what = match &r.action {
            LayerAction::Ran { seconds } => format!("generated in {seconds:.1}s"),
            LayerAction::Reused => "up to date".into(),
            LayerAction::NotApplicable { reason } => format!("does not apply: {reason}"),
        };
        println!("{:<10} {what}", r.layer.name());
        for n in &r.notes {
            println!("           note: {n}");
        }
    })?;
    if render {
        let files = pipe::render_all(&w, &products, &w.root.join("renders"), &ProductId::ALL, &View::all())?;
        println!("renders: {} files in {}", files.len(), w.root.join("renders").display());
    }
    Ok(())
}

pub fn plan(dir: &Path) -> R<()> {
    let w = World::open(dir)?;
    let plan = pipe::plan(&w)?;
    for lp in &plan.layers {
        let by = lp.implementation.map_or(String::new(), |i| format!(" ({})", i.name()));
        let action = match (&lp.not_applicable, lp.needs_run()) {
            (Some(r), _) => format!("skip — {r}"),
            (None, true) => format!("run{by}"),
            (None, false) => format!("reuse{by}"),
        };
        println!("{:<11} {action}", lp.layer.name());
        for p in &lp.products {
            println!("  {:<18} {:<9} {:?}", p.product.name(), p.source, p.state);
        }
    }
    Ok(())
}

pub fn status(dir: &Path) -> R<()> {
    let w = World::open(dir)?;
    let plan = pipe::plan(&w)?;
    let edits = pipe::edit_status(&w);
    println!("{:<18} {:<10} {:<10} edit", "product", "source", "state");
    for lp in &plan.layers {
        for p in &lp.products {
            let e = edits.get(p.product.name()).map_or(String::new(), |e| {
                let mut s = Vec::new();
                if e.rebased {
                    s.push("rebased".to_string());
                }
                if e.orphaned > 0 {
                    s.push(format!("{} orphaned", e.orphaned));
                }
                s.join(", ")
            });
            let state = match p.state {
                State::Fresh => "fresh",
                State::Stale => "stale",
                State::Missing => "missing",
                State::NotApplicable => "n/a",
            };
            println!("{:<18} {:<10} {:<10} {e}", p.product.name(), p.source, state);
        }
    }
    Ok(())
}

pub fn seed(dir: &Path, name: &str, owner: Option<u64>) -> R<()> {
    let w = World::open(dir)?;
    let p = product(name)?;
    let key = owner.map_or(w.manifest.focus, geoforge_core::Key);
    println!("{}", w.manifest.root_seed.for_product(p, key).0);
    Ok(())
}

pub fn stats(dir: &Path) -> R<()> {
    let w = World::open(dir)?;
    let products = w.load_products()?;
    let body = w.load_body()?;
    let grid = w.manifest.grid.build()?;
    let s = geoforge_pipeline::stats::summarize(&products, &grid, body.exterior.radius_km, body.interior.has_surface_liquid_water());
    println!("{}", serde_json::to_string_pretty(&s)?);
    Ok(())
}

pub fn parse_galaxy(name: &str) -> Result<GalaxyPreset, String> {
    GalaxyPreset::from_name(name).ok_or_else(|| format!("unknown galaxy `{name}` (one of: {})", GalaxyPreset::ALL.map(GalaxyPreset::name).join(", ")))
}

/// Products and derived views by name (or a layer's products); everything for "all" or none.
pub fn parse_targets(what: &[String]) -> R<(Vec<ProductId>, Vec<View>)> {
    if what.is_empty() || what.iter().any(|w| w == "all") {
        return Ok((ProductId::ALL.to_vec(), View::all()));
    }
    let (mut ps, mut vs) = (Vec::new(), Vec::new());
    for w in what {
        if let Some(p) = ProductId::from_name(w) {
            ps.push(p);
        } else if let Some(l) = LayerId::from_name(w) {
            ps.extend(l.products());
        } else if let Some(v) = View::from_name(w) {
            vs.push(v);
        } else if let Some(group) = View::group(w) {
            vs.extend(group);
        } else {
            let views = View::all().into_iter().map(View::name).collect::<Vec<_>>().join(", ");
            return Err(format!(
                "unknown product, layer or view `{w}` (views: {views}, or the groups climate_maps, oceans, resources and hazards; products: {})",
                ProductId::ALL.map(ProductId::name).join(", ")
            )
            .into());
        }
    }
    ps.sort();
    ps.dedup();
    vs.sort();
    vs.dedup();
    Ok((ps, vs))
}

pub fn product(name: &str) -> R<ProductId> {
    ProductId::from_name(name).ok_or_else(|| format!("unknown product `{name}`").into())
}

fn read_body(path: &Path) -> R<BodyContext> {
    let text = std::fs::read_to_string(path)?;
    Ok(if path.extension().is_some_and(|e| e == "json") { serde_json::from_str(&text)? } else { toml::from_str(&text)? })
}
