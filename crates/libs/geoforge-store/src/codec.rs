//! Canonical (lossless) storage, written once per shape.
//!
//! Fields and partition owners are raw little-endian cell arrays (`<product>.bin`); every
//! product also has an envelope document (`<product>.json`) holding its schema version,
//! kind, provenance, content hash, grid, unit, and any table or document body.

use crate::{hash_parts, Result, StoreError};
use geoforge_core::grid::GridSpec;
use geoforge_core::product::{Envelope, Hash, ProductData, ProductId, Shape, ValueType};
use geoforge_core::shapes::{Field, Partition, Table};
use geoforge_core::surface::{Biome, Category, CrustType, Lithology};
use serde::{de::DeserializeOwned, Deserialize, Serialize};
use std::path::Path;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Stats {
    pub min: f64,
    pub max: f64,
    pub mean: f64,
}

/// The envelope body of a stored product.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StoredMeta {
    /// Hash of the product's own content.
    pub content: Hash,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub grid: Option<GridSpec>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub value_type: Option<ValueType>,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub unit: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub stats: Option<Stats>,
    /// Partition table or whole-document body.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub data: Option<serde_json::Value>,
}

fn f32_bytes(v: &[f32]) -> Vec<u8> {
    v.iter().flat_map(|x| x.to_le_bytes()).collect()
}
fn u32_bytes(v: &[u32]) -> Vec<u8> {
    v.iter().flat_map(|x| x.to_le_bytes()).collect()
}
fn cat_bytes<C: Category>(v: &[C]) -> Vec<u8> {
    v.iter().map(|c| c.to_u8()).collect()
}

fn v2_bytes(v: &[[f32; 2]]) -> Vec<u8> {
    v.iter().flat_map(|x| x[0].to_le_bytes().into_iter().chain(x[1].to_le_bytes())).collect()
}

/// Cell bytes and document body of a product (shared by storage and hashing).
pub(crate) fn encode(data: &ProductData) -> Result<(Option<Vec<u8>>, Option<serde_json::Value>)> {
    if let Some(f) = data.scalar() {
        return Ok((Some(f32_bytes(&f.values)), None));
    }
    if let Some(f) = data.vector() {
        return Ok((Some(v2_bytes(&f.values)), None));
    }
    Ok(match data {
        ProductData::CrustType(f) => (Some(cat_bytes(&f.values)), None),
        ProductData::Lithology(f) => (Some(cat_bytes(&f.values)), None),
        ProductData::Biome(f) => (Some(cat_bytes(&f.values)), None),
        ProductData::Plates(p) => (Some(u32_bytes(&p.owner.values)), Some(serde_json::to_value(&p.table)?)),
        ProductData::Provinces(p) => (Some(u32_bytes(&p.owner.values)), Some(serde_json::to_value(&p.table)?)),
        ProductData::Lakes(p) => (Some(u32_bytes(&p.owner.values)), Some(serde_json::to_value(&p.table)?)),
        ProductData::Watersheds(p) => (Some(u32_bytes(&p.owner.values)), Some(serde_json::to_value(&p.table)?)),
        ProductData::TectonicEvents(t) | ProductData::GeologyEvents(t) => (None, Some(serde_json::to_value(t)?)),
        ProductData::Rivers(n) => (None, Some(serde_json::to_value(n)?)),
        _ => unreachable!("fields handled above"),
    })
}

/// Hash of a product's content.
///
/// # Errors
/// If the product cannot be serialized.
pub fn content_hash(data: &ProductData) -> Result<Hash> {
    let (bin, doc) = encode(data)?;
    let doc_bytes = doc.map(|d| serde_json::to_vec(&d)).transpose()?.unwrap_or_default();
    Ok(hash_parts(&[data.id().name().as_bytes(), bin.as_deref().unwrap_or(&[]), &doc_bytes]))
}

fn stats_f32(v: &[f32]) -> Stats {
    let (mut mn, mut mx, mut s) = (f64::INFINITY, f64::NEG_INFINITY, 0.0);
    for &x in v {
        let x = f64::from(x);
        mn = mn.min(x);
        mx = mx.max(x);
        s += x;
    }
    Stats { min: mn, max: mx, mean: s / v.len().max(1) as f64 }
}

/// Write a product; returns its content hash.
///
/// # Errors
/// On IO or serialization errors.
pub fn write_canonical(dir: &Path, data: &ProductData, provenance: Hash) -> Result<Hash> {
    std::fs::create_dir_all(dir)?;
    let id = data.id();
    let (bin, doc) = encode(data)?;
    let content = content_hash(data)?;
    let grid = data.grid();
    let stats = data.scalar().map(|f| stats_f32(&f.values));
    let meta = StoredMeta { content, grid, value_type: id.value_type(), unit: id.unit().into(), stats, data: doc };
    let env = Envelope { schema_version: geoforge_core::SCHEMA_VERSION, kind: id, provenance, body: meta };
    if let Some(b) = bin {
        std::fs::write(dir.join(format!("{}.bin", id.name())), b)?;
    }
    std::fs::write(dir.join(format!("{}.json", id.name())), serde_json::to_string_pretty(&env)?)?;
    Ok(content)
}

fn read_env(dir: &Path, p: ProductId) -> Result<Envelope<StoredMeta>> {
    let env: Envelope<StoredMeta> = serde_json::from_str(&std::fs::read_to_string(dir.join(format!("{}.json", p.name())))?)?;
    if env.schema_version != geoforge_core::SCHEMA_VERSION || env.kind != p {
        return Err(StoreError::Invalid(format!("{}: wrong schema or kind", p.name())));
    }
    Ok(env)
}

/// Provenance recorded in a stored product, if present and readable.
#[must_use]
pub fn read_provenance(dir: &Path, p: ProductId) -> Option<Hash> {
    read_env(dir, p).ok().map(|e| e.provenance)
}

/// Content hash recorded in a stored product.
#[must_use]
pub fn read_content_hash(dir: &Path, p: ProductId) -> Option<Hash> {
    read_env(dir, p).ok().map(|e| e.body.content)
}

pub(crate) fn decode_f32(b: &[u8]) -> Vec<f32> {
    b.chunks_exact(4).map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect()
}
pub(crate) fn decode_u32(b: &[u8]) -> Vec<u32> {
    b.chunks_exact(4).map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect()
}
pub(crate) fn decode_cat<C: Category>(b: &[u8]) -> Result<Vec<C>> {
    b.iter().map(|&v| C::from_u8(v).ok_or_else(|| StoreError::Invalid(format!("unknown category value {v}")))).collect()
}

fn table<E: DeserializeOwned>(v: Option<serde_json::Value>) -> Result<Table<E>> {
    Ok(serde_json::from_value(v.unwrap_or(serde_json::json!({})))?)
}

/// Build a product from its cell bytes and document body (canonical or authored).
///
/// # Errors
/// If the bytes or body do not match the product's shape or the grid.
pub fn decode(p: ProductId, grid: GridSpec, bin: Option<&[u8]>, doc: Option<serde_json::Value>) -> Result<ProductData> {
    let n = grid.cell_count();
    let need_bin = |w: usize| -> Result<&[u8]> {
        let b = bin.ok_or_else(|| StoreError::Invalid(format!("{}: missing cell data", p.name())))?;
        if b.len() != n * w {
            return Err(StoreError::Invalid(format!("{}: {} bytes but the grid needs {}", p.name(), b.len(), n * w)));
        }
        Ok(b)
    };
    let part = |b: &[u8]| Field::from_values(grid, decode_u32(b));
    let data = match p {
        ProductId::CrustType => ProductData::CrustType(Field::from_values(grid, decode_cat::<CrustType>(need_bin(1)?)?)),
        ProductId::Lithology => ProductData::Lithology(Field::from_values(grid, decode_cat::<Lithology>(need_bin(1)?)?)),
        ProductId::Biome => ProductData::Biome(Field::from_values(grid, decode_cat::<Biome>(need_bin(1)?)?)),
        ProductId::Plates => ProductData::Plates(Partition { owner: part(need_bin(4)?), table: table(doc)? }),
        ProductId::Provinces => ProductData::Provinces(Partition { owner: part(need_bin(4)?), table: table(doc)? }),
        ProductId::Lakes => ProductData::Lakes(Partition { owner: part(need_bin(4)?), table: table(doc)? }),
        ProductId::Watersheds => ProductData::Watersheds(Partition { owner: part(need_bin(4)?), table: table(doc)? }),
        ProductId::TectonicEvents => ProductData::TectonicEvents(serde_json::from_value(doc.unwrap_or_default())?),
        ProductId::GeologyEvents => ProductData::GeologyEvents(serde_json::from_value(doc.unwrap_or_default())?),
        ProductId::Rivers => ProductData::Rivers(serde_json::from_value(doc.unwrap_or_default())?),
        _ if p.value_type() == Some(ValueType::V2F32) => {
            let v = decode_f32(need_bin(8)?).chunks_exact(2).map(|c| [c[0], c[1]]).collect();
            ProductData::from_vector(p, Field::from_values(grid, v)).expect("vector product")
        }
        _ => ProductData::from_scalar(p, Field::from_values(grid, decode_f32(need_bin(4)?))).expect("scalar product"),
    };
    validate(&data)?;
    Ok(data)
}

/// Shape invariants checked on every load and import.
///
/// # Errors
/// Describes the first broken invariant.
pub fn validate(data: &ProductData) -> Result<()> {
    if let Some(f) = data.scalar() {
        if f.values.iter().any(|v| !v.is_finite()) {
            return Err(StoreError::Invalid("field contains non-finite values".into()));
        }
    }
    if let Some(f) = data.vector() {
        if f.values.iter().any(|v| !v[0].is_finite() || !v[1].is_finite()) {
            return Err(StoreError::Invalid("vector field contains non-finite values".into()));
        }
    }
    match data {
        ProductData::Plates(p) => p.validate()?,
        ProductData::Provinces(p) => p.validate()?,
        ProductData::Lakes(p) => p.validate()?,
        ProductData::Watersheds(p) => p.validate()?,
        _ => {}
    }
    Ok(())
}

/// # Errors
/// On IO, format or validation errors.
pub fn read_canonical(dir: &Path, p: ProductId, grid: GridSpec) -> Result<ProductData> {
    let env = read_env(dir, p)?;
    if let Some(g) = env.body.grid {
        if g != grid {
            return Err(StoreError::Invalid(format!("{}: stored on {g:?}, world uses {grid:?}", p.name())));
        }
    }
    let bin = match p.shape() {
        Shape::Field | Shape::Partition => Some(std::fs::read(dir.join(format!("{}.bin", p.name())))?),
        _ => None,
    };
    decode(p, grid, bin.as_deref(), env.body.data)
}
