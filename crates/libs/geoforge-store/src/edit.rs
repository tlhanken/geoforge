//! Edit views: products as files that mapping and paint software can open, and the way back.
//!
//! | shape                 | view                                                         |
//! |-----------------------|--------------------------------------------------------------|
//! | continuous field      | 16-bit greyscale PNG, value = offset + scale × pixel         |
//! | categorical field     | RGB PNG in the category palette (legend in the sidecar)      |
//! | partition             | RGB PNG, one colour per key, plus `<p>.table.json`           |
//! | timeline, network     | `<p>.json` (the document itself)                             |
//!
//! Every raster comes with `<p>.edit.json` (how pixels map to values), a `.pgw` world file
//! and a `.prj` sphere definition, so QGIS and similar tools place it on the globe.
//!
//! Import turns an edited view back into a product. Against the baseline saved at export
//! time it becomes a **patch** (only the cells you changed: deltas for continuous fields,
//! values for categories and owners), which survives regenerating the layer underneath;
//! or, with `authored`, it replaces the product outright.

use crate::codec::{self, decode_cat};
use crate::{hash_parts, Result, StoreError};
use geoforge_core::grid::{Cell, GridSpec};
use geoforge_core::key::Key;
use geoforge_core::product::{Hash, ProductData, ProductId, Shape};
use geoforge_core::shapes::Field;
use geoforge_core::surface::{Biome, Category, CrustType, Lithology};
use image::{ImageBuffer, Luma, Rgb};
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
        ProductId::Elevation => (-12_000.0, 12_000.0),
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
    format!("#{:02x}{:02x}{:02x}", c[0], c[1], c[2])
}

fn parse_hex(s: &str) -> Result<[u8; 3]> {
    let t = s.trim_start_matches('#');
    if t.len() != 6 {
        return Err(StoreError::Invalid(format!("bad colour {s}")));
    }
    let b = |i: usize| u8::from_str_radix(&t[i..i + 2], 16).map_err(|_| StoreError::Invalid(format!("bad colour {s}")));
    Ok([b(0)?, b(2)?, b(4)?])
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

/// Colours for a partition's keys, unique (nudged on collision).
fn partition_palette(keys: &BTreeSet<u32>) -> BTreeMap<u32, [u8; 3]> {
    let mut used = BTreeSet::new();
    used.insert([0u8, 0, 0]);
    let mut out = BTreeMap::new();
    for &k in keys {
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

fn write_georef(dir: &Path, p: ProductId, grid: GridSpec, body_name: &str, radius_km: f64) -> Result<()> {
    let (w, h) = equirect(grid)?;
    let (dx, dy) = (360.0 / f64::from(w), 180.0 / f64::from(h));
    let pgw = format!("{dx:.12}\n0.0\n0.0\n{:.12}\n{:.12}\n{:.12}\n", -dy, -180.0 + dx / 2.0, 90.0 - dy / 2.0);
    std::fs::write(dir.join(format!("{}.pgw", p.name())), pgw)?;
    let n = body_name.replace('"', "");
    let prj = format!(
        "GEOGCS[\"{n} sphere\",DATUM[\"D_{n}\",SPHEROID[\"{n}\",{:.1},0]],PRIMEM[\"Reference_Meridian\",0],UNIT[\"Degree\",0.0174532925199433]]",
        radius_km * 1000.0
    );
    std::fs::write(dir.join(format!("{}.prj", p.name())), prj)?;
    Ok(())
}

/// Body facts the georeference files need.
#[derive(Debug, Clone)]
pub struct BodyRef {
    pub name: String,
    pub radius_km: f64,
}

/// Write the edit view of a product into `dir`; returns the files written.
///
/// # Errors
/// On IO or encoding errors, or for grids without a raster view.
pub fn export_view(dir: &Path, data: &ProductData, body: &BodyRef) -> Result<Vec<String>> {
    std::fs::create_dir_all(dir)?;
    let p = data.id();
    let mut files = Vec::new();
    let name = p.name();
    let png = dir.join(format!("{name}.png"));
    let (encoding, grid, note) = match data {
        _ if data.vector().is_some() => {
            let f = data.vector().expect("vector");
            let (w, h) = equirect(f.grid)?;
            let (lo, hi) = gray_range(p).expect("vector field");
            let scale = (hi - lo) / 65535.0;
            let img: ImageBuffer<Rgb<u16>, Vec<u16>> =
                ImageBuffer::from_raw(w, h, f.values.iter().flat_map(|v| [to_u16(f64::from(v[0]), lo, scale), to_u16(f64::from(v[1]), lo, scale), 0]).collect()).expect("size");
            img.save(&png)?;
            (
                Encoding::Vector16 { offset: lo, scale, unit: p.unit().into() },
                f.grid,
                format!("16-bit RGB: red = eastward, green = northward [{}], each offset + scale * channel; blue unused.", p.unit()),
            )
        }
        _ if data.scalar().is_some() => {
            let f = data.scalar().expect("scalar");
            let (w, h) = equirect(f.grid)?;
            let (lo, hi) = gray_range(p).expect("continuous field");
            let scale = (hi - lo) / 65535.0;
            let img: ImageBuffer<Luma<u16>, Vec<u16>> =
                ImageBuffer::from_raw(w, h, f.values.iter().map(|&v| to_u16(f64::from(v), lo, scale)).collect()).expect("size");
            img.save(&png)?;
            (Encoding::Gray16 { offset: lo, scale, unit: p.unit().into() }, f.grid, format!("16-bit greyscale: value [{}] = offset + scale * pixel. Values outside [{lo}, {hi}] are clipped.", p.unit()))
        }
        ProductData::CrustType(f) => (write_category(&png, f)?, f.grid, "Paint only with the palette colours; other colours snap to the nearest one.".into()),
        ProductData::Biome(f) => (write_category(&png, f)?, f.grid, "Paint only with the palette colours; other colours snap to the nearest one.".into()),
        ProductData::Lithology(f) => (write_category(&png, f)?, f.grid, "Paint only with the palette colours; other colours snap to the nearest one.".into()),
        ProductData::Plates(pt) => write_partition(dir, p, &pt.owner, &pt.table.rows)?,
        ProductData::Provinces(pt) => write_partition(dir, p, &pt.owner, &pt.table.rows)?,
        ProductData::Lakes(pt) => write_partition(dir, p, &pt.owner, &pt.table.rows)?,
        ProductData::Watersheds(pt) => write_partition(dir, p, &pt.owner, &pt.table.rows)?,
        _ => {
            let (_, doc) = codec::encode(data)?;
            std::fs::write(dir.join(format!("{name}.json")), serde_json::to_string_pretty(&doc)?)?;
            files.push(format!("{name}.json"));
            return Ok(files);
        }
    };
    files.push(format!("{name}.png"));
    if p.shape() == Shape::Partition {
        files.push(format!("{name}.table.json"));
    }
    write_georef(dir, p, grid, &body.name, body.radius_km)?;
    files.push(format!("{name}.pgw"));
    files.push(format!("{name}.prj"));
    let side = Sidecar { product: p, grid, encoding, note };
    std::fs::write(dir.join(format!("{name}.edit.json")), serde_json::to_string_pretty(&side)?)?;
    files.push(format!("{name}.edit.json"));
    Ok(files)
}

fn to_u16(v: f64, lo: f64, scale: f64) -> u16 {
    ((v - lo) / scale).round().clamp(0.0, 65535.0) as u16
}

fn write_category<C: Category>(png: &Path, f: &Field<C>) -> Result<Encoding> {
    let (w, h) = equirect(f.grid)?;
    let img: ImageBuffer<Rgb<u8>, Vec<u8>> = ImageBuffer::from_raw(w, h, f.values.iter().flat_map(|c| c.color()).collect()).expect("size");
    img.save(png)?;
    Ok(Encoding::Palette { entries: category_palette::<C>() })
}

fn write_partition<E: Serialize>(dir: &Path, p: ProductId, owner: &Field<u32>, rows: &BTreeMap<Key, E>) -> Result<(Encoding, GridSpec, String)> {
    let (w, h) = equirect(owner.grid)?;
    let keys: BTreeSet<u32> = rows.keys().map(|k| k.0 as u32).collect();
    let pal = partition_palette(&keys);
    let img: ImageBuffer<Rgb<u8>, Vec<u8>> = ImageBuffer::from_raw(
        w,
        h,
        owner.values.iter().flat_map(|o| if *o == 0 { [0, 0, 0] } else { pal.get(o).copied().unwrap_or([0, 0, 0]) }).collect(),
    )
    .expect("size");
    img.save(dir.join(format!("{}.png", p.name())))?;
    let table: Vec<TableRow> = rows
        .iter()
        .map(|(k, r)| Ok(TableRow { key: *k, color: hex(pal[&(k.0 as u32)]), row: Some(serde_json::to_value(r)?) }))
        .collect::<Result<_>>()?;
    std::fs::write(dir.join(format!("{}.table.json", p.name())), serde_json::to_string_pretty(&table)?)?;
    let entries = table.iter().map(|t| PaletteEntry { value: t.key.0 as u32, color: t.color.clone(), name: format!("key {}", t.key) }).collect();
    Ok((
        Encoding::Palette { entries },
        owner.grid,
        "One colour per key; black = no owner. To add a region, add a row with a new key and colour to the table file and paint with that colour.".into(),
    ))
}

/// The raw content of an edited view, before it becomes a product.
#[derive(Debug, Clone)]
pub enum ViewValues {
    F32(Vec<f32>),
    V2(Vec<[f32; 2]>),
    U8(Vec<u8>),
    U32(Vec<u32>),
    Doc(serde_json::Value),
}

#[derive(Debug, Clone)]
pub struct View {
    pub product: ProductId,
    pub values: ViewValues,
    /// Partition rows by key (None = row to be filled in by the layer).
    pub table: BTreeMap<Key, Option<serde_json::Value>>,
    /// Pixels whose colour was not an exact palette colour.
    pub snapped: usize,
}

/// The encoding a view of `p` uses when no sidecar is present (a map drawn from scratch).
///
/// # Errors
/// For partitions, which need their table file.
pub fn default_encoding(p: ProductId) -> Result<Encoding> {
    match p {
        ProductId::CrustType => Ok(Encoding::Palette { entries: category_palette::<CrustType>() }),
        ProductId::Biome => Ok(Encoding::Palette { entries: category_palette::<Biome>() }),
        ProductId::Lithology => Ok(Encoding::Palette { entries: category_palette::<Lithology>() }),
        _ if p.value_type() == Some(geoforge_core::product::ValueType::V2F32) => {
            let (lo, hi) = gray_range(p).expect("vector field");
            Ok(Encoding::Vector16 { offset: lo, scale: (hi - lo) / 65535.0, unit: p.unit().into() })
        }
        _ if p.shape() == Shape::Field => {
            let (lo, hi) = gray_range(p).expect("continuous field");
            Ok(Encoding::Gray16 { offset: lo, scale: (hi - lo) / 65535.0, unit: p.unit().into() })
        }
        _ if p.shape() == Shape::Partition => Ok(Encoding::Palette { entries: Vec::new() }),
        _ => Ok(Encoding::Json),
    }
}

/// Read an edited view from `dir`.
///
/// # Errors
/// If files are missing, the image size does not match the grid, or colours are unknown.
pub fn read_view(dir: &Path, p: ProductId, grid: GridSpec) -> Result<View> {
    let name = p.name();
    if matches!(p.shape(), Shape::Timeline | Shape::Network) {
        let doc: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(dir.join(format!("{name}.json")))?)?;
        return Ok(View { product: p, values: ViewValues::Doc(doc), table: BTreeMap::new(), snapped: 0 });
    }
    let (w, h) = equirect(grid)?;
    let side_path = dir.join(format!("{name}.edit.json"));
    let encoding = if side_path.exists() {
        let s: Sidecar = serde_json::from_str(&std::fs::read_to_string(&side_path)?)?;
        if s.product != p {
            return Err(StoreError::Invalid(format!("{} describes {}, not {name}", side_path.display(), s.product)));
        }
        s.encoding
    } else {
        default_encoding(p)?
    };
    let img = image::open(dir.join(format!("{name}.png")))?;
    if img.width() != w || img.height() != h {
        return Err(StoreError::Invalid(format!(
            "{name}.png is {}x{} but the world grid is {w}x{h} (resize it in your editor, keeping the full-globe extent)",
            img.width(),
            img.height()
        )));
    }
    let mut table = BTreeMap::new();
    let (values, snapped) = match encoding {
        Encoding::Gray16 { offset, scale, .. } => {
            let g = img.into_luma16();
            (ViewValues::F32(g.into_raw().into_iter().map(|v| (offset + scale * f64::from(v)) as f32).collect()), 0)
        }
        Encoding::Vector16 { offset, scale, .. } => {
            let g = img.into_rgb16();
            (ViewValues::V2(g.pixels().map(|px| [(offset + scale * f64::from(px.0[0])) as f32, (offset + scale * f64::from(px.0[1])) as f32]).collect()), 0)
        }
        Encoding::Palette { mut entries } => {
            if p.shape() == Shape::Partition {
                let rows: Vec<TableRow> = serde_json::from_str(&std::fs::read_to_string(dir.join(format!("{name}.table.json")))?)?;
                entries = rows.iter().map(|r| PaletteEntry { value: r.key.0 as u32, color: r.color.clone(), name: String::new() }).collect();
                for r in rows {
                    if r.key.0 == 0 || r.key.0 > u64::from(u32::MAX) {
                        return Err(StoreError::Invalid(format!("{name}: key {} is out of range (1..=4294967295)", r.key)));
                    }
                    table.insert(r.key, r.row);
                }
            }
            let rgb = img.into_rgb8();
            let mut pal: Vec<([u8; 3], u32)> = entries.iter().map(|e| Ok((parse_hex(&e.color)?, e.value))).collect::<Result<_>>()?;
            if p.shape() == Shape::Partition {
                pal.push(([0, 0, 0], 0));
            }
            let exact: BTreeMap<[u8; 3], u32> = pal.iter().copied().collect();
            let mut snapped = 0;
            let mut cache: BTreeMap<[u8; 3], u32> = BTreeMap::new();
            let vals: Vec<u32> = rgb
                .pixels()
                .map(|px| {
                    let c = px.0;
                    if let Some(&v) = exact.get(&c) {
                        return v;
                    }
                    snapped += 1;
                    *cache.entry(c).or_insert_with(|| {
                        pal.iter()
                            .min_by_key(|(q, _)| (0..3).map(|i| (i32::from(c[i]) - i32::from(q[i])).pow(2)).sum::<i32>())
                            .map_or(0, |(_, v)| *v)
                    })
                })
                .collect();
            if p.shape() == Shape::Partition {
                (ViewValues::U32(vals), snapped)
            } else {
                (ViewValues::U8(vals.into_iter().map(|v| v as u8).collect()), snapped)
            }
        }
        Encoding::Json => return Err(StoreError::Invalid(format!("{name}: a raster product cannot use a JSON encoding"))),
    };
    Ok(View { product: p, values, table, snapped })
}

/// Turn a complete view into a product. `fill` supplies rows for partition keys that have
/// none (it gets the key and its cells); a key it cannot fill is an error.
///
/// # Errors
/// If the view does not form a valid product.
pub fn view_to_product(view: &View, grid: GridSpec, fill: &dyn Fn(ProductId, Key, &[Cell]) -> Option<serde_json::Value>) -> Result<ProductData> {
    let p = view.product;
    let n = grid.cell_count();
    match &view.values {
        ViewValues::Doc(d) => codec::decode(p, grid, None, Some(d.clone())),
        ViewValues::F32(v) => {
            let b: Vec<u8> = v.iter().flat_map(|x| x.to_le_bytes()).collect();
            codec::decode(p, grid, Some(&b), None)
        }
        ViewValues::V2(v) => {
            let b: Vec<u8> = v.iter().flat_map(|x| x[0].to_le_bytes().into_iter().chain(x[1].to_le_bytes())).collect();
            codec::decode(p, grid, Some(&b), None)
        }
        ViewValues::U8(v) => {
            match p {
                ProductId::Lithology => decode_cat::<Lithology>(v).map(|_| ())?,
                ProductId::CrustType => decode_cat::<CrustType>(v).map(|_| ())?,
                ProductId::Biome => decode_cat::<Biome>(v).map(|_| ())?,
                _ => {}
            }
            codec::decode(p, grid, Some(v), None)
        }
        ViewValues::U32(owner) => {
            if owner.len() != n {
                return Err(StoreError::Invalid("owner length mismatch".into()));
            }
            let mut cells: BTreeMap<u32, Vec<Cell>> = BTreeMap::new();
            for (c, &o) in owner.iter().enumerate() {
                if o != 0 {
                    cells.entry(o).or_default().push(c as Cell);
                }
            }
            let mut rows = serde_json::Map::new();
            let mut unfilled = Vec::new();
            for (&k, cs) in &cells {
                let key = Key(u64::from(k));
                let row = match view.table.get(&key) {
                    Some(Some(r)) => Some(r.clone()),
                    _ => fill(p, key, cs),
                };
                match row {
                    Some(r) => {
                        rows.insert(k.to_string(), r);
                    }
                    None => unfilled.push(k),
                }
            }
            if !unfilled.is_empty() {
                return Err(StoreError::Invalid(format!("{}: keys {unfilled:?} own cells but have no row in the table file", p.name())));
            }
            let b: Vec<u8> = owner.iter().flat_map(|x| x.to_le_bytes()).collect();
            codec::decode(p, grid, Some(&b), Some(serde_json::Value::Object(rows)))
        }
    }
}

// ---------------------------------------------------------------- patches

/// Cells changed by an edit, as runs `[start, len]` over cell indices.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Patch {
    pub product: ProductId,
    /// Content hash of the product the patch was made against.
    pub base: Hash,
    pub runs: Vec<[u32; 2]>,
    pub values: PatchValues,
    /// Partition rows added or changed by the edit.
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub rows: BTreeMap<Key, Option<serde_json::Value>>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "op", content = "v", rename_all = "snake_case")]
pub enum PatchValues {
    /// Added to the base value (continuous fields).
    Delta(Vec<f32>),
    /// Added to the base vector (vector fields).
    Delta2(Vec<[f32; 2]>),
    /// Replaces the base value (categories, partition owners).
    Set(Vec<u32>),
}

impl Patch {
    #[must_use]
    pub fn cell_count(&self) -> usize {
        self.runs.iter().map(|r| r[1] as usize).sum()
    }
    fn cells(&self) -> impl Iterator<Item = u32> + '_ {
        self.runs.iter().flat_map(|r| r[0]..r[0] + r[1])
    }
    /// Hash of the patch itself (part of downstream provenance).
    ///
    /// # Errors
    /// On serialization errors.
    pub fn hash(&self) -> Result<Hash> {
        Ok(hash_parts(&[b"patch", &serde_json::to_vec(self)?]))
    }
}

fn runs_of(cells: &[u32]) -> Vec<[u32; 2]> {
    let mut runs: Vec<[u32; 2]> = Vec::new();
    for &c in cells {
        match runs.last_mut() {
            Some(r) if r[0] + r[1] == c => r[1] += 1,
            _ => runs.push([c, 1]),
        }
    }
    runs
}

/// Diff an edited view against the product it was exported from.
///
/// # Errors
/// If the view and baseline are different products or shapes.
pub fn diff(baseline: &ProductData, view: &View) -> Result<Patch> {
    let p = baseline.id();
    if view.product != p {
        return Err(StoreError::Invalid("view and baseline are different products".into()));
    }
    let base = codec::content_hash(baseline)?;
    let mut cells = Vec::new();
    let values = match (baseline, &view.values) {
        (d, ViewValues::V2(v)) if d.vector().is_some() => {
            let f = d.vector().expect("vector");
            let (lo, hi) = gray_range(p).expect("vector field");
            let scale = (hi - lo) / 65535.0;
            let mut out = Vec::new();
            for (c, (b, e)) in f.values.iter().zip(v).enumerate() {
                // A component whose pixel value is unchanged keeps its exact base value.
                let changed = [0, 1].map(|k| to_u16(f64::from(b[k]), lo, scale) != to_u16(f64::from(e[k]), lo, scale));
                if changed[0] || changed[1] {
                    cells.push(c as u32);
                    out.push([0, 1].map(|k| if changed[k] { e[k] - b[k] } else { 0.0 }));
                }
            }
            PatchValues::Delta2(out)
        }
        (d, ViewValues::F32(v)) if d.scalar().is_some() => {
            let f = d.scalar().expect("scalar");
            let (lo, hi) = gray_range(p).expect("continuous field");
            let scale = (hi - lo) / 65535.0;
            let mut d = Vec::new();
            for (c, (&b, &e)) in f.values.iter().zip(v).enumerate() {
                // Compare in pixel space so the export's own rounding is never an "edit".
                if to_u16(f64::from(b), lo, scale) != to_u16(f64::from(e), lo, scale) {
                    cells.push(c as u32);
                    d.push(e - b);
                }
            }
            PatchValues::Delta(d)
        }
        (ProductData::CrustType(f), ViewValues::U8(v)) => PatchValues::Set(set_diff(f.values.iter().map(|c| c.to_u8()), v.iter().copied(), &mut cells)),
        (ProductData::Biome(f), ViewValues::U8(v)) => PatchValues::Set(set_diff(f.values.iter().map(|c| c.to_u8()), v.iter().copied(), &mut cells)),
        (ProductData::Lithology(f), ViewValues::U8(v)) => PatchValues::Set(set_diff(f.values.iter().map(|c| c.to_u8()), v.iter().copied(), &mut cells)),
        (ProductData::Plates(pt), ViewValues::U32(v)) => PatchValues::Set(set_diff(pt.owner.values.iter().copied(), v.iter().copied(), &mut cells)),
        (ProductData::Provinces(pt), ViewValues::U32(v)) => PatchValues::Set(set_diff(pt.owner.values.iter().copied(), v.iter().copied(), &mut cells)),
        (ProductData::Lakes(pt), ViewValues::U32(v)) => PatchValues::Set(set_diff(pt.owner.values.iter().copied(), v.iter().copied(), &mut cells)),
        (ProductData::Watersheds(pt), ViewValues::U32(v)) => PatchValues::Set(set_diff(pt.owner.values.iter().copied(), v.iter().copied(), &mut cells)),
        _ => return Err(StoreError::Invalid(format!("{}: only raster views can be imported as patches; import documents as authored", p.name()))),
    };
    // Rows the edit added or changed.
    let mut rows = BTreeMap::new();
    if p.shape() == Shape::Partition {
        let (_, doc) = codec::encode(baseline)?;
        let base_rows: BTreeMap<String, serde_json::Value> = serde_json::from_value(doc.unwrap_or_default())?;
        for (k, r) in &view.table {
            let before = base_rows.get(&k.0.to_string());
            if r.as_ref() != before {
                rows.insert(*k, r.clone());
            }
        }
    }
    Ok(Patch { product: p, base, runs: runs_of(&cells), values, rows })
}

fn set_diff(base: impl Iterator<Item = impl Into<u32>>, edited: impl Iterator<Item = impl Into<u32>>, cells: &mut Vec<u32>) -> Vec<u32> {
    let mut out = Vec::new();
    for (c, (b, e)) in base.zip(edited).enumerate() {
        let (b, e) = (b.into(), e.into());
        if b != e {
            cells.push(c as u32);
            out.push(e);
        }
    }
    out
}

/// Result of applying a patch to a (possibly regenerated) product.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Applied {
    /// The base changed since the patch was made (the patch was re-applied on top).
    pub rebased: bool,
    /// Edits that could not be applied (cells outside the grid, values with no meaning).
    pub orphaned: usize,
}

/// Apply a patch in place. `fill` supplies rows for keys the patch introduces without one.
///
/// # Errors
/// If the patch is for a different product or kind of value.
pub fn apply(data: &mut ProductData, patch: &Patch, fill: &dyn Fn(ProductId, Key, &[Cell]) -> Option<serde_json::Value>) -> Result<Applied> {
    let p = data.id();
    if patch.product != p {
        return Err(StoreError::Invalid(format!("patch is for {}, not {}", patch.product, p.name())));
    }
    let rebased = codec::content_hash(data)? != patch.base;
    let mut orphaned = 0;
    match (data, &patch.values) {
        (data, PatchValues::Delta2(d)) if data.vector().is_some() => {
            let f = data.vector_mut().expect("vector");
            for (c, dv) in patch.cells().zip(d) {
                match f.values.get_mut(c as usize) {
                    Some(v) => {
                        v[0] += dv[0];
                        v[1] += dv[1];
                    }
                    None => orphaned += 1,
                }
            }
        }
        (data, PatchValues::Delta(d)) if data.scalar().is_some() => {
            let f = data.scalar_mut().expect("scalar");
            for (c, &dv) in patch.cells().zip(d) {
                match f.values.get_mut(c as usize) {
                    Some(v) => *v += dv,
                    None => orphaned += 1,
                }
            }
        }
        (ProductData::CrustType(f), PatchValues::Set(s)) => orphaned += set_cat(&mut f.values, patch.cells(), s),
        (ProductData::Biome(f), PatchValues::Set(s)) => orphaned += set_cat(&mut f.values, patch.cells(), s),
        (ProductData::Lithology(f), PatchValues::Set(s)) => orphaned += set_cat(&mut f.values, patch.cells(), s),
        (ProductData::Plates(pt), PatchValues::Set(s)) => orphaned += set_part(p, &mut pt.owner, &mut pt.table.rows, patch, s, fill)?,
        (ProductData::Provinces(pt), PatchValues::Set(s)) => orphaned += set_part(p, &mut pt.owner, &mut pt.table.rows, patch, s, fill)?,
        (ProductData::Lakes(pt), PatchValues::Set(s)) => orphaned += set_part(p, &mut pt.owner, &mut pt.table.rows, patch, s, fill)?,
        (ProductData::Watersheds(pt), PatchValues::Set(s)) => orphaned += set_part(p, &mut pt.owner, &mut pt.table.rows, patch, s, fill)?,
        _ => return Err(StoreError::Invalid(format!("{}: patch values do not fit the product", p.name()))),
    }
    Ok(Applied { rebased, orphaned })
}

fn set_cat<C: Category>(values: &mut [C], cells: impl Iterator<Item = u32>, s: &[u32]) -> usize {
    let mut orphaned = 0;
    for (c, &v) in cells.zip(s) {
        match (values.get_mut(c as usize), u8::try_from(v).ok().and_then(C::from_u8)) {
            (Some(slot), Some(cat)) => *slot = cat,
            _ => orphaned += 1,
        }
    }
    orphaned
}

fn set_part<E: serde::de::DeserializeOwned>(
    p: ProductId,
    owner: &mut Field<u32>,
    rows: &mut BTreeMap<Key, E>,
    patch: &Patch,
    s: &[u32],
    fill: &dyn Fn(ProductId, Key, &[Cell]) -> Option<serde_json::Value>,
) -> Result<usize> {
    let mut orphaned = 0;
    let mut prev: BTreeMap<Cell, u32> = BTreeMap::new();
    for (c, &v) in patch.cells().zip(s) {
        match owner.values.get_mut(c as usize) {
            Some(o) => {
                prev.insert(c, *o);
                *o = v;
            }
            None => orphaned += 1,
        }
    }
    for (k, r) in &patch.rows {
        if let Some(r) = r {
            rows.insert(*k, serde_json::from_value(r.clone())?);
        }
    }
    // Keys that own cells but have no row: fill, else give the cells back to their base owner.
    let mut need: BTreeMap<u32, Vec<Cell>> = BTreeMap::new();
    for (c, &o) in owner.values.iter().enumerate() {
        if o != 0 && !rows.contains_key(&Key(u64::from(o))) {
            need.entry(o).or_default().push(c as Cell);
        }
    }
    for (k, cs) in need {
        match fill(p, Key(u64::from(k)), &cs).map(serde_json::from_value::<E>) {
            Some(Ok(r)) => {
                rows.insert(Key(u64::from(k)), r);
            }
            _ => {
                orphaned += cs.len();
                for c in cs {
                    let back = prev.get(&c).copied().unwrap_or(0);
                    owner.values[c as usize] = if rows.contains_key(&Key(u64::from(back))) { back } else { 0 };
                }
            }
        }
    }
    let used: BTreeSet<u32> = owner.values.iter().copied().collect();
    rows.retain(|k, _| used.contains(&(k.0 as u32)));
    Ok(orphaned)
}

/// A fill function that fills nothing.
#[must_use]
pub fn no_fill() -> impl Fn(ProductId, Key, &[Cell]) -> Option<serde_json::Value> {
    |_, _, _| None
}

/// Combine an existing patch with edits made on top of it. `newer` must have been diffed
/// against the product with `older` applied; the result keeps `older`'s base.
///
/// # Errors
/// If the patches are for different products or kinds of value.
pub fn merge(older: &Patch, newer: &Patch) -> Result<Patch> {
    if older.product != newer.product {
        return Err(StoreError::Invalid("patches are for different products".into()));
    }
    let collect = |p: &Patch| -> Vec<u32> { p.cells().collect() };
    let (oc, nc) = (collect(older), collect(newer));
    let values = match (&older.values, &newer.values) {
        (PatchValues::Delta(a), PatchValues::Delta(b)) => {
            let mut m: BTreeMap<u32, f32> = oc.iter().copied().zip(a.iter().copied()).collect();
            for (c, d) in nc.iter().zip(b) {
                *m.entry(*c).or_insert(0.0) += d;
            }
            m.retain(|_, d| d.abs() > 0.0);
            let cells: Vec<u32> = m.keys().copied().collect();
            return Ok(Patch { product: older.product, base: older.base, runs: runs_of(&cells), values: PatchValues::Delta(m.into_values().collect()), rows: BTreeMap::new() });
        }
        (PatchValues::Delta2(a), PatchValues::Delta2(b)) => {
            let mut m: BTreeMap<u32, [f32; 2]> = oc.iter().copied().zip(a.iter().copied()).collect();
            for (c, d) in nc.iter().zip(b) {
                let e = m.entry(*c).or_insert([0.0, 0.0]);
                e[0] += d[0];
                e[1] += d[1];
            }
            let cells: Vec<u32> = m.keys().copied().collect();
            return Ok(Patch { product: older.product, base: older.base, runs: runs_of(&cells), values: PatchValues::Delta2(m.into_values().collect()), rows: BTreeMap::new() });
        }
        (PatchValues::Set(a), PatchValues::Set(b)) => {
            let mut m: BTreeMap<u32, u32> = oc.iter().copied().zip(a.iter().copied()).collect();
            for (c, v) in nc.iter().zip(b) {
                m.insert(*c, *v);
            }
            m
        }
        _ => return Err(StoreError::Invalid("patches hold different kinds of value".into())),
    };
    let cells: Vec<u32> = values.keys().copied().collect();
    let mut rows = older.rows.clone();
    rows.extend(newer.rows.clone());
    Ok(Patch { product: older.product, base: older.base, runs: runs_of(&cells), values: PatchValues::Set(values.into_values().collect()), rows })
}
