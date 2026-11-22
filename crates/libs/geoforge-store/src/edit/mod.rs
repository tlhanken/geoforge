//! Edit views: products as files that mapping and paint software can open, and the way back.
//!
//! | shape                 | view                                                         |
//! |-----------------------|--------------------------------------------------------------|
//! | continuous field      | 16-bit greyscale PNG, value = offset + scale × pixel         |
//! | categorical field     | RGB PNG in the category palette (legend in the sidecar)      |
//! | partition             | RGB PNG, one colour per key, plus `<p>.table.json`           |
//! | catalog (rock, biome) | as a partition; each row's own `color` is its palette colour |
//! | timeline, network     | `<p>.json` (the document itself)                             |
//!
//! Every raster comes with `<p>.edit.json` (how pixels map to values), a `.pgw` world file
//! and a `.prj` sphere definition, so QGIS and similar tools place it on the globe. Rasters
//! can be exported larger than the grid (`scale`) and with polar cap views (`<p>.north.png`,
//! `<p>.south.png`, see [`crate::polar`]); edited rasters come back at any 2:1 size.
//!
//! Import turns an edited view back into a product. Against the baseline saved at export
//! time it becomes a **patch** (only the cells you changed: deltas for continuous fields,
//! values for categories and owners), which survives regenerating the layer underneath;
//! or, with `authored`, it replaces the product outright.

mod import;
mod patch;

pub use import::{default_encoding, read_view, view_to_product, Fill, View, ViewValues};
pub use patch::{apply, diff, merge, no_fill, Applied, Patch, PatchValues};

use crate::codec::{self, decode_cat};
use crate::polar::Cap;
use crate::raster::{Depth, Pixels};
use crate::{hash_parts, Result, StoreError};
use geoforge_core::grid::{Cell, GridSpec};
use geoforge_core::key::Key;
use geoforge_core::product::{Hash, ProductData, ProductId, Shape};
use geoforge_core::shapes::Field;
use geoforge_core::surface::{hex_color, Category, CrustType};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

/// How a view's pixels encode values.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum Encoding {
    /// `value = offset + scale * pixel` (16-bit).
    Gray16 { offset: f64, scale: f64, unit: String },
    /// Vector field: red = east, green = north, each `offset + scale * channel` (16-bit RGB).
    Vector16 { offset: f64, scale: f64, unit: String },
    /// Exact colours for each value (category u8 or partition key).
    Palette { entries: Vec<PaletteEntry> },
    /// The product document as JSON.
    Json,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PaletteEntry {
    pub value: u32,
    pub color: String,
    pub name: String,
}

/// `<p>.edit.json`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Sidecar {
    pub product: ProductId,
    pub grid: GridSpec,
    pub encoding: Encoding,
    pub note: String,
}

/// One row of `<p>.table.json`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TableRow {
    pub key: Key,
    pub color: String,
    /// The row's attributes. May be omitted for a new key when the layer can fill it in.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub row: Option<serde_json::Value>,
}

/// Default 16-bit range of each continuous field.
#[must_use]
pub fn gray_range(p: ProductId) -> Option<(f64, f64)> {
    Some(match p {
        ProductId::Elevation | ProductId::RoughElevation => (-12_000.0, 12_000.0),
        ProductId::CrustAge => (0.0, 4_500.0),
        ProductId::Temperature => (-90.0, 70.0),
        ProductId::TemperatureRange => (0.0, 80.0),
        ProductId::Precipitation => (0.0, 12_000.0),
        ProductId::RockAge => (0.0, 4_500.0),
        ProductId::SstAnomaly => (-20.0, 20.0),
        ProductId::PrecipitationWarmShare => (0.0, 1.0),
        ProductId::WindJan | ProductId::WindJul => (-40.0, 40.0),
        ProductId::OceanCurrents => (-3.0, 3.0),
        _ => return None,
    })
}

fn hex(c: [u8; 3]) -> String {
    hex_color::format(c)
}

fn parse_hex(s: &str) -> Result<[u8; 3]> {
    hex_color::parse(s).map_err(StoreError::Invalid)
}

/// A stable, well-spread colour for partition key `k`.
#[must_use]
pub fn key_color(k: u32) -> [u8; 3] {
    let h = (f64::from(k) * 0.618_033_988_749_895).fract();
    let s = 0.55 + 0.35 * ((f64::from(k) * 0.381_966).fract());
    let v = 0.70 + 0.28 * ((f64::from(k) * 0.723_606).fract());
    hsv(h, s, v)
}

pub(crate) fn hsv(h: f64, s: f64, v: f64) -> [u8; 3] {
    let i = (h * 6.0).floor();
    let f = h * 6.0 - i;
    let (p, q, t) = (v * (1.0 - s), v * (1.0 - f * s), v * (1.0 - (1.0 - f) * s));
    let (r, g, b) = match (i as i64).rem_euclid(6) {
        0 => (v, t, p),
        1 => (q, v, p),
        2 => (p, v, t),
        3 => (p, q, v),
        4 => (t, p, v),
        _ => (v, p, q),
    };
    [(r * 255.0).round() as u8, (g * 255.0).round() as u8, (b * 255.0).round() as u8]
}

fn category_palette<C: Category>() -> Vec<PaletteEntry> {
    C::ALL.iter().map(|c| PaletteEntry { value: u32::from(c.to_u8()), color: hex(c.color()), name: c.name().into() }).collect()
}

/// Colours for a partition's keys: a row's own colour where it has one, otherwise a stable
/// generated colour, unique (nudged on collision).
pub fn partition_palette(keys: &BTreeMap<u32, Option<[u8; 3]>>) -> BTreeMap<u32, [u8; 3]> {
    let mut used: BTreeSet<[u8; 3]> = keys.values().flatten().copied().collect();
    used.insert([0u8, 0, 0]);
    let mut out = BTreeMap::new();
    for (&k, own) in keys {
        if let Some(c) = own {
            out.insert(k, *c);
            continue;
        }
        let mut c = key_color(k);
        let mut n = 0u32;
        while used.contains(&c) {
            n += 1;
            c = key_color(k.wrapping_add(n.wrapping_mul(7919)));
        }
        used.insert(c);
        out.insert(k, c);
    }
    out
}

fn equirect(grid: GridSpec) -> Result<(u32, u32)> {
    match grid {
        GridSpec::Equirect { width, height } => Ok((width, height)),
        GridSpec::Goldberg { .. } => Err(StoreError::Invalid("edit views need an equirectangular grid (Goldberg views will resample; not implemented yet)".into())),
    }
}

/// Write the world file (`<stem>.pgw`) and projection (`<stem>.prj`) that place an
/// equirectangular PNG on the body's sphere.
///
/// # Errors
/// On IO errors, or for grids without a raster form.
pub fn write_georef(dir: &Path, stem: &str, grid: GridSpec, body: &BodyRef) -> Result<()> {
    let (w, h) = equirect(grid)?;
    georef_sized(dir, stem, w, h, body)
}

/// [`write_georef`] for an image `w` × `h` covering the globe.
fn georef_sized(dir: &Path, stem: &str, w: u32, h: u32, body: &BodyRef) -> Result<()> {
    let (dx, dy) = (360.0 / f64::from(w), 180.0 / f64::from(h));
    let pgw = format!("{dx:.12}\n0.0\n0.0\n{:.12}\n{:.12}\n{:.12}\n", -dy, -180.0 + dx / 2.0, 90.0 - dy / 2.0);
    std::fs::write(dir.join(format!("{stem}.pgw")), pgw)?;
    let n = body.name.replace('"', "");
    let prj = format!("GEOGCS[\"{n} sphere\",DATUM[\"D_{n}\",SPHEROID[\"{n}\",{:.1},0]],PRIMEM[\"Reference_Meridian\",0],UNIT[\"Degree\",0.0174532925199433]]", body.radius_km * 1000.0);
    std::fs::write(dir.join(format!("{stem}.prj")), prj)?;
    Ok(())
}

/// Body facts the georeference files need.
#[derive(Debug, Clone)]
pub struct BodyRef {
    pub name: String,
    pub radius_km: f64,
}

/// How rasters are exported.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ViewOptions {
    /// Pixels per cell along each axis (each cell repeated, so unpainted areas come back
    /// exactly).
    pub scale: u32,
    /// Also write polar cap views.
    pub caps: bool,
}

impl Default for ViewOptions {
    fn default() -> Self {
        Self { scale: 1, caps: false }
    }
}

/// A product as a raster view at grid size: its pixels, how they encode values, a note for
/// the sidecar, and (for partitions) the table file's rows.
#[derive(Debug, Clone)]
pub struct Raster {
    pub pixels: Pixels,
    pub encoding: Encoding,
    pub grid: GridSpec,
    pub note: String,
    pub table: Option<Vec<TableRow>>,
}

/// The raster view of `data`, or `None` for products kept as JSON documents.
///
/// # Errors
/// For grids without a raster form.
pub fn raster(data: &ProductData) -> Result<Option<Raster>> {
    let p = data.id();
    let plain = |pixels, encoding, grid, note: String| Ok(Some(Raster { pixels, encoding, grid, note, table: None }));
    match data {
        _ if data.vector().is_some() => {
            let f = data.vector().expect("vector");
            equirect(f.grid)?;
            let (lo, hi) = gray_range(p).expect("vector field");
            let scale = (hi - lo) / 65535.0;
            let px = Pixels::Rgb16(f.values.iter().map(|v| [to_u16(f64::from(v[0]), lo, scale), to_u16(f64::from(v[1]), lo, scale), 0]).collect());
            let note = format!("16-bit RGB: red = eastward, green = northward [{}], each offset + scale * channel; blue unused.", p.unit());
            plain(px, Encoding::Vector16 { offset: lo, scale, unit: p.unit().into() }, f.grid, note)
        }
        _ if data.scalar().is_some() => {
            let f = data.scalar().expect("scalar");
            equirect(f.grid)?;
            let (lo, hi) = gray_range(p).expect("continuous field");
            let scale = (hi - lo) / 65535.0;
            let px = Pixels::Gray16(f.values.iter().map(|&v| to_u16(f64::from(v), lo, scale)).collect());
            let note = format!("16-bit greyscale: value [{}] = offset + scale * pixel. Values outside [{lo}, {hi}] are clipped.", p.unit());
            plain(px, Encoding::Gray16 { offset: lo, scale, unit: p.unit().into() }, f.grid, note)
        }
        ProductData::CrustType(f) => {
            equirect(f.grid)?;
            let px = Pixels::Rgb8(f.values.iter().map(|c| c.color()).collect());
            plain(px, Encoding::Palette { entries: category_palette::<CrustType>() }, f.grid, "Paint only with the palette colours; other colours snap to the nearest one.".into())
        }
        ProductData::Lithology(pt) => partition_raster(p, &pt.owner, &pt.table.rows, |r| Some(r.color)).map(Some),
        ProductData::Biome(pt) => partition_raster(p, &pt.owner, &pt.table.rows, |r| Some(r.color)).map(Some),
        ProductData::Plates(pt) => partition_raster(p, &pt.owner, &pt.table.rows, |_| None).map(Some),
        ProductData::Provinces(pt) => partition_raster(p, &pt.owner, &pt.table.rows, |_| None).map(Some),
        ProductData::Lakes(pt) => partition_raster(p, &pt.owner, &pt.table.rows, |_| None).map(Some),
        ProductData::Watersheds(pt) => partition_raster(p, &pt.owner, &pt.table.rows, |_| None).map(Some),
        _ => Ok(None),
    }
}

/// Write the edit view of a product into `dir`; returns the files written.
///
/// # Errors
/// On IO or encoding errors, or for grids without a raster view.
pub fn export_view(dir: &Path, data: &ProductData, body: &BodyRef, opts: &ViewOptions) -> Result<Vec<String>> {
    std::fs::create_dir_all(dir)?;
    let p = data.id();
    let name = p.name();
    let Some(r) = raster(data)? else {
        let (_, doc) = codec::encode(data)?;
        std::fs::write(dir.join(format!("{name}.json")), serde_json::to_string_pretty(&doc)?)?;
        return Ok(vec![format!("{name}.json")]);
    };
    let (w, h) = equirect(r.grid)?;
    let k = opts.scale.max(1);
    let mut files = vec![format!("{name}.png")];
    r.pixels.enlarge(w, h, k).save(w * k, h * k, &dir.join(format!("{name}.png")))?;
    if let Some(table) = &r.table {
        std::fs::write(dir.join(format!("{name}.table.json")), serde_json::to_string_pretty(table)?)?;
        files.push(format!("{name}.table.json"));
    }
    georef_sized(dir, name, w * k, h * k, body)?;
    files.extend([format!("{name}.pgw"), format!("{name}.prj")]);
    if opts.caps {
        for north in [true, false] {
            let cap = Cap::for_grid(north, h, k);
            cap.render(&r.pixels, w, h).save(cap.side, cap.side, &dir.join(format!("{name}.{}.png", cap.name())))?;
            cap.write_georef(dir, name, body)?;
            files.extend(["png", "pgw", "prj"].map(|e| format!("{name}.{}.{e}", cap.name())));
        }
    }
    let side = Sidecar { product: p, grid: r.grid, encoding: r.encoding, note: r.note };
    std::fs::write(dir.join(format!("{name}.edit.json")), serde_json::to_string_pretty(&side)?)?;
    files.push(format!("{name}.edit.json"));
    Ok(files)
}

fn to_u16(v: f64, lo: f64, scale: f64) -> u16 {
    ((v - lo) / scale).round().clamp(0.0, 65535.0) as u16
}

fn partition_raster<E: Serialize>(p: ProductId, owner: &Field<u32>, rows: &BTreeMap<Key, E>, own_color: impl Fn(&E) -> Option<[u8; 3]>) -> Result<Raster> {
    equirect(owner.grid)?;
    let keys: BTreeMap<u32, Option<[u8; 3]>> = rows.iter().map(|(k, r)| (k.0 as u32, own_color(r))).collect();
    let pal = partition_palette(&keys);
    let pixels = Pixels::Rgb8(owner.values.iter().map(|o| if *o == 0 { [0, 0, 0] } else { pal.get(o).copied().unwrap_or([0, 0, 0]) }).collect());
    let table: Vec<TableRow> = rows.iter().map(|(k, r)| Ok(TableRow { key: *k, color: hex(pal[&(k.0 as u32)]), row: Some(serde_json::to_value(r)?) })).collect::<Result<_>>()?;
    let name = |t: &TableRow| t.row.as_ref().and_then(|r| r.get("name")).and_then(|n| n.as_str()).map_or_else(|| format!("key {}", t.key), str::to_string);
    let entries = table.iter().map(|t| PaletteEntry { value: t.key.0 as u32, color: t.color.clone(), name: name(t) }).collect();
    let note = if p.is_catalog() {
        "One colour per kind (the row's own colour). To add a kind, add a row with a new key, its colour and its attributes to the table file, then paint with that colour."
    } else {
        "One colour per key; black = no owner. To add a region, add a row with a new key and colour to the table file and paint with that colour."
    };
    Ok(Raster { pixels, encoding: Encoding::Palette { entries }, grid: owner.grid, note: note.into(), table: Some(table) })
}
