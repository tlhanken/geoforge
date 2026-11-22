//! Reading an edited view back: resampled to the grid, polar caps folded in, pixels decoded
//! to values; and a whole view turned into a product.

use super::*;

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
    /// The image's size, when it was not the grid's and was resampled.
    pub resampled_from: Option<(u32, u32)>,
    /// Cells taken from edited polar cap views.
    pub cap_cells: usize,
}

/// The encoding a view of `p` uses when no sidecar is present (a map drawn from scratch).
///
/// # Errors
/// For partitions, which need their table file.
pub fn default_encoding(p: ProductId) -> Result<Encoding> {
    match p {
        ProductId::CrustType => Ok(Encoding::Palette { entries: category_palette::<CrustType>() }),
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

/// Read an edited view from `dir`. An image of another 2:1 size is resampled to the grid;
/// edited polar caps beside it (which need the export's `baseline`) override it near the
/// poles.
///
/// # Errors
/// If files are missing, the image does not cover the globe, colours are unknown, or caps
/// are present without a baseline.
pub fn read_view(dir: &Path, p: ProductId, grid: GridSpec, baseline: Option<&ProductData>) -> Result<View> {
    let name = p.name();
    let blank = |values| View { product: p, values, table: BTreeMap::new(), snapped: 0, resampled_from: None, cap_cells: 0 };
    if matches!(p.shape(), Shape::Timeline | Shape::Network | Shape::Table) {
        let doc: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(dir.join(format!("{name}.json")))?)?;
        return Ok(blank(ViewValues::Doc(doc)));
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
    let depth = match encoding {
        Encoding::Gray16 { .. } => Depth::Gray16,
        Encoding::Vector16 { .. } => Depth::Rgb16,
        Encoding::Palette { .. } => Depth::Rgb8,
        Encoding::Json => return Err(StoreError::Invalid(format!("{name}: a raster product cannot use a JSON encoding"))),
    };
    let (mut pixels, iw, ih) = Pixels::load(&dir.join(format!("{name}.png")), depth)?;
    let mut view = blank(ViewValues::U8(Vec::new()));
    if (iw, ih) != (w, h) {
        if iw.abs_diff(2 * ih) > 1 {
            return Err(StoreError::Invalid(format!("{name}.png is {iw}x{ih}: a view must cover the globe, twice as wide as tall (the grid is {w}x{h})")));
        }
        pixels = pixels.resample((iw, ih), (w, h));
        view.resampled_from = Some((iw, ih));
    }
    for north in [true, false] {
        let path = dir.join(format!("{name}.{}.png", if north { "north" } else { "south" }));
        if !path.exists() {
            continue;
        }
        let base = baseline.ok_or_else(|| StoreError::Invalid(format!("{} needs the export's baseline: import it as a patch, or delete it", path.display())))?;
        let exported = raster(base)?.ok_or_else(|| StoreError::Invalid(format!("{name} has no raster view")))?.pixels;
        let (edited, cw, ch) = Pixels::load(&path, depth)?;
        if cw != ch {
            return Err(StoreError::Invalid(format!("{} is {cw}x{ch}; a cap view is square", path.display())));
        }
        let cap = Cap { north, side: cw };
        view.cap_cells += cap.merge(&edited, &cap.render(&exported, w, h), &mut pixels, w, h)?;
    }
    let (values, snapped) = match (encoding, pixels) {
        (Encoding::Gray16 { offset, scale, .. }, Pixels::Gray16(g)) => (ViewValues::F32(g.into_iter().map(|v| (offset + scale * f64::from(v)) as f32).collect()), 0),
        (Encoding::Vector16 { offset, scale, .. }, Pixels::Rgb16(g)) => {
            (ViewValues::V2(g.into_iter().map(|px| [(offset + scale * f64::from(px[0])) as f32, (offset + scale * f64::from(px[1])) as f32]).collect()), 0)
        }
        (Encoding::Palette { mut entries }, Pixels::Rgb8(rgb)) => {
            if p.shape() == Shape::Partition {
                let rows: Vec<TableRow> = serde_json::from_str(&std::fs::read_to_string(dir.join(format!("{name}.table.json")))?)?;
                // A catalog row's own colour is its palette colour.
                let color = |r: &TableRow| {
                    let own = r.row.as_ref().filter(|_| p.is_catalog()).and_then(|v| v.get("color")).and_then(|c| c.as_str());
                    own.unwrap_or(&r.color).to_string()
                };
                entries = rows.iter().map(|r| PaletteEntry { value: r.key.0 as u32, color: color(r), name: String::new() }).collect();
                for r in rows {
                    if r.key.0 == 0 || r.key.0 > u64::from(u32::MAX) {
                        return Err(StoreError::Invalid(format!("{name}: key {} is out of range (1..=4294967295)", r.key)));
                    }
                    view.table.insert(r.key, r.row);
                }
            }
            let mut pal: Vec<([u8; 3], u32)> = entries.iter().map(|e| Ok((parse_hex(&e.color)?, e.value))).collect::<Result<_>>()?;
            if p.shape() == Shape::Partition {
                pal.push(([0, 0, 0], 0));
            }
            let exact: BTreeMap<[u8; 3], u32> = pal.iter().copied().collect();
            let mut snapped = 0;
            let mut cache: BTreeMap<[u8; 3], u32> = BTreeMap::new();
            let vals: Vec<u32> = rgb
                .into_iter()
                .map(|c| {
                    if let Some(&v) = exact.get(&c) {
                        return v;
                    }
                    snapped += 1;
                    *cache.entry(c).or_insert_with(|| pal.iter().min_by_key(|(q, _)| (0..3).map(|i| (i32::from(c[i]) - i32::from(q[i])).pow(2)).sum::<i32>()).map_or(0, |(_, v)| *v))
                })
                .collect();
            if p.shape() == Shape::Partition {
                (ViewValues::U32(vals), snapped)
            } else {
                (ViewValues::U8(vals.into_iter().map(|v| v as u8).collect()), snapped)
            }
        }
        _ => unreachable!("pixels are loaded at the encoding's depth"),
    };
    view.values = values;
    view.snapped = snapped;
    Ok(view)
}

/// Attributes for a partition key an edit introduced without a table row: `(product, key,
/// cells) -> row`, or `None` to leave the row empty.
pub type Fill<'a> = dyn Fn(ProductId, Key, &[Cell]) -> Option<serde_json::Value> + 'a;

/// Turn a complete view into a product. `fill` supplies rows for partition keys that have
/// none (it gets the key and its cells); a key it cannot fill is an error.
///
/// # Errors
/// If the view does not form a valid product.
pub fn view_to_product(view: &View, grid: GridSpec, fill: &Fill) -> Result<ProductData> {
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
            if p == ProductId::CrustType {
                decode_cat::<CrustType>(v)?;
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
            if p.is_catalog() {
                for (k, r) in &view.table {
                    if let Some(r) = r {
                        rows.entry(k.0.to_string()).or_insert_with(|| r.clone());
                    }
                }
            }
            let b: Vec<u8> = owner.iter().flat_map(|x| x.to_le_bytes()).collect();
            codec::decode(p, grid, Some(&b), Some(serde_json::Value::Object(rows)))
        }
    }
}
