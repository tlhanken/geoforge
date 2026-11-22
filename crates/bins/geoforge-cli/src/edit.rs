//! Editing: `export` views to paint, `import` them back, `reset` edits, `render` previews.

use crate::world::{parse_targets, product};
use crate::R;
use geoforge_core::product::ProductId;
use geoforge_pipeline as pipe;
use geoforge_store::World;
use std::path::{Path, PathBuf};

pub fn export(dir: &Path, what: &[String], to: Option<PathBuf>, scale: u32, caps: bool) -> R<()> {
    let w = World::open(dir)?;
    let to = to.unwrap_or_else(|| pipe::edit_folder(&w));
    let (ps, vs) = parse_targets(what)?;
    let ps: Vec<ProductId> = ps.into_iter().filter(|p| w.has_product(*p)).collect();
    for (p, files) in pipe::export(&w, &ps, &to, &geoforge_store::edit::ViewOptions { scale, caps })? {
        println!("{:<18} {}", p.name(), files.join(" "));
    }
    for (v, files) in pipe::export_views(&w, &vs, &to)? {
        println!("{:<18} {} (derived; not imported)", v.name(), files.join(" "));
    }
    if !ps.is_empty() {
        println!("edit views in {} — edit, then: geoforge import {} <product>", to.display(), dir.display());
    }
    Ok(())
}

pub fn import(dir: &Path, name: &str, from: Option<PathBuf>, authored: bool) -> R<()> {
    let mut w = World::open(dir)?;
    let from = from.unwrap_or_else(|| pipe::edit_folder(&w));
    let r = pipe::import(&mut w, product(name)?, &from, authored)?;
    println!("{}: imported as {} ({} cells changed, {} rows changed)", r.product.name(), r.mode, r.changed_cells, r.rows_changed);
    if let Some((iw, ih)) = r.resampled_from {
        println!("  the image was {iw}x{ih}; resampled to the grid");
    }
    if r.cap_cells > 0 {
        println!("  {} cells came from the polar cap views", r.cap_cells);
    }
    if r.snapped > 0 {
        println!("  {} pixels were not exact palette colours and were snapped to the nearest one", r.snapped);
    }
    println!("next: geoforge generate {} (downstream layers will update)", dir.display());
    Ok(())
}

pub fn reset(dir: &Path, name: &str) -> R<()> {
    let mut w = World::open(dir)?;
    pipe::reset(&mut w, product(name)?)?;
    println!("{name}: edits removed; it will be generated");
    Ok(())
}

pub fn render(dir: &Path, what: &[String], to: Option<PathBuf>) -> R<()> {
    let w = World::open(dir)?;
    let products = w.load_products()?;
    let (ps, vs) = parse_targets(what)?;
    let to = to.unwrap_or_else(|| w.root.join("renders"));
    for f in pipe::render_all(&w, &products, &to, &ps, &vs)? {
        println!("{}", f.display());
    }
    Ok(())
}
