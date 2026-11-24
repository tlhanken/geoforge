//! Making a world and running its layers: `new`, `generate`, `plan`, `status`, `stats`,
//! `seed`, and the names commands take.

use crate::system::show;
use crate::R;
use geoforge_core::grid::GridSpec;
use geoforge_core::product::{LayerId, ProductId};
use geoforge_core::Address;
use geoforge_core::{BodyContext, Seed};
use geoforge_cosmos::{GalaxyPreset, Recipe};
use geoforge_pipeline::{self as pipe, GenerateOptions, LayerAction, NewWorld, State, View};
use geoforge_store::{UniverseDir, World};
use std::path::{Path, PathBuf};

pub fn new_world(dir: &Path, seed: u64, at: Option<&str>, width: u32, body: Option<PathBuf>, archetype: Option<geoforge_core::Archetype>, focus: Option<&str>) -> R<()> {
    if seed > i64::MAX as u64 {
        return Err("seed must be at most 9223372036854775807".into());
    }
    let grid = GridSpec::Equirect { width, height: width / 2 };
    if body.is_some() || archetype.is_some() {
        // A world on its own: a real body with the real Solar System at its present, or a
        // hand-made body.
        let mut body = body.as_deref().map(read_body).transpose()?;
        let mut system = None;
        if let Some(a) = archetype {
            let (sys, b) = geoforge_cosmos::sol::archetype(Seed(seed), a);
            (system, body) = (Some(sys), Some(b));
        }
        let w = pipe::create_world(dir, &NewWorld { seed: Seed(seed), recipe: Recipe::Habitable, grid, body, system, focus: None })?;
        show(&w)?;
        println!("\ncreated {} — next: geoforge generate {}", dir.display(), dir.display());
        return Ok(());
    }
    let existing = UniverseDir::is_at(dir).then(|| UniverseDir::open(dir)).transpose()?;
    let address = match at {
        Some(text) => parse_address(text, existing.as_ref())?,
        None => existing.as_ref().map_or(Address::universe(Seed(seed)), |u| u.doc.universe),
    };
    let made = existing.is_none();
    let mut u = match existing {
        Some(u) => u,
        None => UniverseDir::create(dir, address)?,
    };
    let w = match pipe::create_at(&mut u, &address, focus, grid) {
        Ok(w) => w,
        Err(e) => {
            if made {
                // Leave no empty universe behind.
                std::fs::remove_file(dir.join("universe.json"))?;
                let _ = std::fs::remove_dir(dir);
            }
            return Err(e.into());
        }
    };
    show(&w)?;
    let at = w.manifest.address.map(|a| a.to_string()).unwrap_or_default();
    println!("\n{} at {at} in {} — next: geoforge generate {}", w.manifest.focus_name, dir.display(), dir.display());
    Ok(())
}

/// The bodies made in universe `dir`, or focus it on the one `at`.
pub fn focus(dir: &Path, at: Option<&str>) -> R<()> {
    let mut u = UniverseDir::open(dir)?;
    if let Some(text) = at {
        let a = parse_address(text, Some(&u))?;
        let w = u.world(&a)?.ok_or_else(|| format!("{a} is not made yet: geoforge new {} --at {a}", dir.display()))?;
        u.set_focus(a)?;
        println!("focus: {} ({a})", w.manifest.focus_name);
        return Ok(());
    }
    println!("universe {}", u.doc.universe);
    for w in u.worlds()? {
        let a = w.manifest.address;
        let mark = if a == u.doc.focus { "*" } else { " " };
        println!("{mark} {:<20} {}", w.manifest.focus_name, a.map(|a| a.to_string()).unwrap_or_default());
    }
    Ok(())
}

/// An address as typed; in a universe's directory the universe part may be left out.
pub fn parse_address(text: &str, universe: Option<&UniverseDir>) -> R<Address> {
    match (text.parse::<Address>(), universe) {
        (Ok(a), _) => Ok(a),
        (Err(e), Some(u)) => format!("{}/{}", u.doc.universe, text.trim_start_matches('/')).parse::<Address>().map_err(|_| e.into()),
        (Err(e), None) => Err(e.into()),
    }
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

pub fn report(dir: &Path, to: Option<PathBuf>, places: &[String], day: f64, days: usize, width: u32) -> R<()> {
    let w = World::open(dir)?;
    let places = places.iter().map(|p| parse_place(p)).collect::<R<Vec<_>>>()?;
    let opts = pipe::report::Options { width, places, day, days, system_text: Some(crate::system::show_text(&w)?) };
    let r = pipe::report::build(&w, &opts)?;
    let to = to.unwrap_or_else(|| w.root.join("report"));
    let files = geoforge_store::report::write(&to, &r.html, &r.images)?;
    println!("{}: {} files in {} (open index.html)", r.title, files.len(), to.display());
    Ok(())
}

pub fn compare(real: &Path, models: &[PathBuf], labels: &[String], probes: &[String], sources: &[String], to: Option<PathBuf>, width: u32) -> R<()> {
    let reference = World::open(real)?;
    let worlds = models.iter().map(|m| World::open(m)).collect::<Result<Vec<_>, _>>()?;
    let models: Vec<pipe::compare::Model> = worlds
        .iter()
        .zip(models)
        .enumerate()
        .map(|(i, (w, dir))| pipe::compare::Model {
            label: labels.get(i).cloned().unwrap_or_else(|| dir.file_name().map_or_else(|| dir.display().to_string(), |n| n.to_string_lossy().into_owned())),
            world: w,
        })
        .collect();
    let probes = probes.iter().map(|p| parse_probe(p)).collect::<R<Vec<_>>>()?;
    let sources = sources
        .iter()
        .map(|s| s.split_once('=').map(|(p, t)| (p.trim().to_string(), t.trim().to_string())).ok_or_else(|| format!("source `{s}`: expected product=text")))
        .collect::<Result<Vec<_>, _>>()?;
    let c = pipe::compare::build(&reference, &models, &pipe::compare::Options { width, probes, sources })?;
    let to = to.unwrap_or_else(|| reference.root.join("compare"));
    let files = geoforge_store::report::write(&to, &c.html, &c.images)?;
    println!("{}: {} sections, {} files in {} (open index.html)", c.title, c.sections.len(), files.len(), to.display());
    Ok(())
}

/// "Lima:-12.05,-77.05" or "Rockies:37..50,-115..-104".
fn parse_probe(s: &str) -> R<pipe::compare::Probe> {
    let (name, at) = s.rsplit_once(':').ok_or_else(|| format!("probe `{s}`: expected name:lat,lon or name:lat0..lat1,lon0..lon1"))?;
    let (lat, lon) = at.split_once(',').ok_or_else(|| format!("probe `{s}`: expected lat,lon"))?;
    let range = |t: &str| -> R<(f64, f64)> {
        Ok(match t.split_once("..") {
            Some((a, b)) => {
                let (a, b): (f64, f64) = (a.trim().parse()?, b.trim().parse()?);
                (a.min(b), a.max(b))
            }
            None => {
                let v: f64 = t.trim().parse()?;
                (v, v)
            }
        })
    };
    Ok(pipe::compare::Probe { name: name.to_string(), lat: range(lat)?, lon: range(lon)? })
}

/// "London:51.5,-0.1" (or "51.5,-0.1", named by its coordinates).
fn parse_place(s: &str) -> R<pipe::report::Place> {
    let (name, at) = s.rsplit_once(':').map_or((None, s), |(n, a)| (Some(n), a));
    let (lat, lon) = at.split_once(',').ok_or_else(|| format!("place `{s}`: expected name:lat,lon"))?;
    let (lat, lon): (f64, f64) = (lat.trim().parse()?, lon.trim().parse()?);
    Ok(pipe::report::Place { name: name.map_or_else(|| format!("{lat:.1}°, {lon:.1}°"), str::to_string), lat, lon })
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
