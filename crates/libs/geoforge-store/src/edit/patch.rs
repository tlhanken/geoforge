//! Patches: the cells an edit changed, kept apart from the generated product so they ride
//! on top when the layer is regenerated.

use super::*;

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
        (ProductData::Biome(pt), ViewValues::U32(v)) => PatchValues::Set(set_diff(pt.owner.values.iter().copied(), v.iter().copied(), &mut cells)),
        (ProductData::Lithology(pt), ViewValues::U32(v)) => PatchValues::Set(set_diff(pt.owner.values.iter().copied(), v.iter().copied(), &mut cells)),
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
pub fn apply(data: &mut ProductData, patch: &Patch, fill: &Fill) -> Result<Applied> {
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
        (ProductData::Biome(pt), PatchValues::Set(s)) => orphaned += set_part(p, &mut pt.owner, &mut pt.table.rows, patch, s, fill)?,
        (ProductData::Lithology(pt), PatchValues::Set(s)) => orphaned += set_part(p, &mut pt.owner, &mut pt.table.rows, patch, s, fill)?,
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

fn set_part<E: serde::de::DeserializeOwned>(p: ProductId, owner: &mut Field<u32>, rows: &mut BTreeMap<Key, E>, patch: &Patch, s: &[u32], fill: &Fill) -> Result<usize> {
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
    if !p.is_catalog() {
        let used: BTreeSet<u32> = owner.values.iter().copied().collect();
        rows.retain(|k, _| used.contains(&(k.0 as u32)));
    }
    Ok(orphaned)
}

/// A fill function that fills nothing.
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
