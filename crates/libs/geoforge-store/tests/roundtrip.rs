use geoforge_core::grid::GridSpec;
use geoforge_core::key::Key;
use geoforge_core::product::{Hash, ProductData, ProductId};
use geoforge_core::shapes::{Field, Partition, Table};
use geoforge_core::surface::{BiomeClass, Lake};
use geoforge_store::codec;
use geoforge_store::edit::{self, BodyRef};

const G: GridSpec = GridSpec::Equirect { width: 64, height: 32 };

fn tmp(name: &str) -> std::path::PathBuf {
    let d = std::env::temp_dir().join(format!("geoforge-store-test-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn elevation() -> ProductData {
    ProductData::Elevation(Field::from_values(G, (0..G.cell_count()).map(|i| (i as f32 * 7.3) % 9000.0 - 4500.0).collect()))
}

fn lakes() -> ProductData {
    let mut t = Table::new();
    t.insert(Key(1), Lake { surface_elevation_m: 120.0, endorheic: false });
    t.insert(Key(2), Lake { surface_elevation_m: -20.0, endorheic: true });
    let owner = (0..G.cell_count() as u32).map(|i| if i % 50 == 0 { 1 } else if i % 77 == 0 { 2 } else { 0 }).collect();
    ProductData::Lakes(Partition { owner: Field::from_values(G, owner), table: t })
}

/// A catalog: three kinds, the third unused.
fn biomes() -> ProductData {
    let class = |slug: &str, color| BiomeClass { slug: slug.into(), name: slug.into(), color, tags: Vec::new() };
    let mut t = Table::new();
    t.insert(Key(1), class("ocean", [28, 62, 120]));
    t.insert(Key(2), class("savanna", [200, 180, 80]));
    t.insert(Key(3), class("glowcap", [90, 58, 138]));
    let owner = (0..G.cell_count() as u32).map(|i| if i % 3 == 0 { 1 } else { 2 }).collect();
    ProductData::Biome(Partition { owner: Field::from_values(G, owner), table: t })
}

#[test]
fn canonical_roundtrip_is_lossless() {
    let d = tmp("canon");
    for data in [elevation(), lakes(), biomes()] {
        let h = codec::write_canonical(&d, &data, Hash([7; 32])).unwrap();
        let back = codec::read_canonical(&d, data.id(), G).unwrap();
        assert_eq!(back, data);
        assert_eq!(codec::content_hash(&back).unwrap(), h);
        assert_eq!(codec::read_provenance(&d, data.id()), Some(Hash([7; 32])));
    }
    // A different grid is refused.
    assert!(codec::read_canonical(&d, ProductId::Elevation, GridSpec::Equirect { width: 32, height: 16 }).is_err());
}

#[test]
fn untouched_export_imports_as_empty_patch() {
    let d = tmp("noop");
    let body = BodyRef { name: "Test".into(), radius_km: 6371.0 };
    for data in [elevation(), lakes(), biomes()] {
        edit::export_view(&d, &data, &body).unwrap();
        let v = edit::read_view(&d, data.id(), G).unwrap();
        let p = edit::diff(&data, &v).unwrap();
        assert_eq!(p.cell_count(), 0, "{}", data.id());
        assert!(p.rows.is_empty());
    }
    assert!(d.join("elevation.pgw").exists() && d.join("elevation.prj").exists());
}

#[test]
fn painted_edits_become_a_patch_that_rebases() {
    let d = tmp("paint");
    let body = BodyRef { name: "Test".into(), radius_km: 6371.0 };
    let base = elevation();
    edit::export_view(&d, &base, &body).unwrap();
    // "Paint" a 4x4 block 1000 m higher.
    let mut img = image::open(d.join("elevation.png")).unwrap().into_luma16();
    let step = 24000.0 / 65535.0;
    for y in 10..14 {
        for x in 20..24 {
            let p = img.get_pixel_mut(x, y);
            p.0[0] = (f64::from(p.0[0]) + 1000.0 / step).round() as u16;
        }
    }
    img.save(d.join("elevation.png")).unwrap();
    let v = edit::read_view(&d, ProductId::Elevation, G).unwrap();
    let patch = edit::diff(&base, &v).unwrap();
    assert_eq!(patch.cell_count(), 16);
    assert_eq!(patch.runs.len(), 4);

    // Same base: applies exactly.
    let mut same = base.clone();
    let a = edit::apply(&mut same, &patch, &edit::no_fill()).unwrap();
    assert!(!a.rebased && a.orphaned == 0);
    let (ProductData::Elevation(b), ProductData::Elevation(s)) = (&base, &same) else { unreachable!() };
    let c = 12 * 64 + 21;
    assert!((s.values[c] - b.values[c] - 1000.0).abs() < 0.5);
    assert_eq!(s.values[0], b.values[0]);

    // Regenerated base (everything shifted 50 m): the edit rides on top.
    let ProductData::Elevation(mut moved) = base.clone() else { unreachable!() };
    moved.values.iter_mut().for_each(|v| *v += 50.0);
    let mut moved = ProductData::Elevation(moved);
    let a = edit::apply(&mut moved, &patch, &edit::no_fill()).unwrap();
    assert!(a.rebased);
    let ProductData::Elevation(m) = &moved else { unreachable!() };
    assert!((m.values[c] - b.values[c] - 1050.0).abs() < 0.5);
}

#[test]
fn partition_paint_with_new_key() {
    let d = tmp("part");
    let body = BodyRef { name: "Test".into(), radius_km: 6371.0 };
    let base = lakes();
    edit::export_view(&d, &base, &body).unwrap();
    // Add a lake: new table row with a new colour, painted over a block.
    let mut rows: Vec<edit::TableRow> = serde_json::from_str(&std::fs::read_to_string(d.join("lakes.table.json")).unwrap()).unwrap();
    rows.push(edit::TableRow { key: Key(9), color: "#ff00ff".into(), row: Some(serde_json::json!({"surface_elevation_m": 300.0, "endorheic": false})) });
    std::fs::write(d.join("lakes.table.json"), serde_json::to_string(&rows).unwrap()).unwrap();
    let mut img = image::open(d.join("lakes.png")).unwrap().into_rgb8();
    for y in 2..5 {
        for x in 2..5 {
            img.put_pixel(x, y, image::Rgb([255, 0, 255]));
        }
    }
    // An anti-aliased edge pixel snaps to the nearest palette colour.
    img.put_pixel(5, 3, image::Rgb([250, 4, 250]));
    img.save(d.join("lakes.png")).unwrap();
    let v = edit::read_view(&d, ProductId::Lakes, G).unwrap();
    assert_eq!(v.snapped, 1);
    let patch = edit::diff(&base, &v).unwrap();
    assert_eq!(patch.rows.len(), 1);
    let mut out = base.clone();
    edit::apply(&mut out, &patch, &edit::no_fill()).unwrap();
    let ProductData::Lakes(l) = &out else { unreachable!() };
    l.validate().unwrap();
    assert_eq!(l.owner.values[3 * 64 + 3], 9);
    assert_eq!(l.owner.values[3 * 64 + 5], 9);
    assert!(l.table.get(Key(9)).is_some());

    // Authored from the same view.
    let authored = edit::view_to_product(&v, G, &edit::no_fill()).unwrap();
    assert_eq!(authored, out);
}

#[test]
fn vector_field_views_round_trip() {
    let d = tmp("vec");
    let body = BodyRef { name: "Test".into(), radius_km: 6371.0 };
    let v: Vec<[f32; 2]> = (0..G.cell_count()).map(|i| [((i % 17) as f32 - 8.0) * 1.5, ((i % 11) as f32 - 5.0) * 2.0]).collect();
    let data = ProductData::WindJul(Field::from_values(G, v));
    let h = codec::write_canonical(&d, &data, Hash([1; 32])).unwrap();
    assert_eq!(codec::read_canonical(&d, ProductId::WindJul, G).unwrap(), data);
    assert_eq!(codec::content_hash(&data).unwrap(), h);
    edit::export_view(&d, &data, &body).unwrap();
    let view = edit::read_view(&d, ProductId::WindJul, G).unwrap();
    let patch = edit::diff(&data, &view).unwrap();
    assert_eq!(patch.cell_count(), 0);
    // Turn the wind around in one cell: a two-component delta.
    let mut img = image::open(d.join("wind_jul.png")).unwrap().into_rgb16();
    let px = img.get_pixel_mut(5, 5);
    px.0[0] = 65535 - px.0[0];
    img.save(d.join("wind_jul.png")).unwrap();
    let view = edit::read_view(&d, ProductId::WindJul, G).unwrap();
    let patch = edit::diff(&data, &view).unwrap();
    assert_eq!(patch.cell_count(), 1);
    let mut out = data.clone();
    edit::apply(&mut out, &patch, &edit::no_fill()).unwrap();
    let (ProductData::WindJul(a), ProductData::WindJul(b)) = (&data, &out) else { unreachable!() };
    let c = 5 * 64 + 5;
    assert!((a.values[c][0] + b.values[c][0]).abs() < 0.01, "{:?} {:?}", a.values[c], b.values[c]);
    assert_eq!(a.values[c][1], b.values[c][1]);
}

/// Catalog views paint in each row's own colour, keep unused rows, and take new kinds from
/// added rows.
#[test]
fn catalog_views_keep_and_add_kinds() {
    let d = tmp("catalog");
    let body = BodyRef { name: "Test".into(), radius_km: 6371.0 };
    let base = biomes();
    edit::export_view(&d, &base, &body).unwrap();
    let side: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(d.join("biome.edit.json")).unwrap()).unwrap();
    let entries = side["encoding"]["entries"].as_array().unwrap();
    assert!(entries.iter().any(|e| e["color"] == "#5a3a8a" && e["name"] == "glowcap"), "row colour is the palette colour");

    // Add a kind and paint the first row with it.
    let table_path = d.join("biome.table.json");
    let mut rows: Vec<serde_json::Value> = serde_json::from_str(&std::fs::read_to_string(&table_path).unwrap()).unwrap();
    rows.push(serde_json::json!({ "key": 101, "color": "#ff00aa", "row": { "slug": "ashlands", "name": "Ashlands", "color": "#ff00aa" } }));
    std::fs::write(&table_path, serde_json::to_string(&rows).unwrap()).unwrap();
    let png = d.join("biome.png");
    let mut img = image::open(&png).unwrap().into_rgb8();
    for x in 0..64 {
        img.put_pixel(x, 0, image::Rgb([0xff, 0x00, 0xaa]));
    }
    img.save(&png).unwrap();

    let v = edit::read_view(&d, ProductId::Biome, G).unwrap();
    let edited = edit::view_to_product(&v, G, &edit::no_fill()).unwrap();
    let ProductData::Biome(pt) = &edited else { panic!() };
    assert_eq!(pt.table.rows.len(), 4, "unused glowcap kept, ashlands added");
    assert_eq!(pt.owner.values[0], 101);
    assert_eq!(pt.table.rows[&Key(101)].color, [0xff, 0x00, 0xaa]);

    // As a patch on the original, the unused row survives too.
    let patch = edit::diff(&base, &v).unwrap();
    assert_eq!(patch.cell_count(), 64);
    let mut again = base.clone();
    edit::apply(&mut again, &patch, &edit::no_fill()).unwrap();
    let ProductData::Biome(pt) = &again else { panic!() };
    assert!(pt.table.rows.contains_key(&Key(3)) && pt.table.rows.contains_key(&Key(101)));
}
