//! Known cases, tectonics → geology → elevation: Earth's real plates (Bird's PB2002 outlines
//! and Euler poles) and crust (continental where ETOPO1 is above −2000 m) go in; boundary
//! types, crust ages, provinces and the heightmap come out and are checked against PB2002's
//! boundary classes and ETOPO1. Half-degree grid: arcs and ridge offsets need ~50 km cells.
//!
//! `cargo test --test known_earth_tectonics -- --nocapture` prints the scorecard;
//! `cargo test --test known_earth_tectonics -- --ignored` asserts the known gaps.

mod known;

use geoforge_core::product::{LayerId, ProductData, Products};
use geoforge_core::surface::{BoundaryKind, CrustType, GeologicProvince as G, Plate};
use geoforge_core::{Archetype, Cell, Field, Grid, Key, Partition, Table, UnitVec3, Vec3};
use geoforge_surface::tectonics::boundaries;
use known::{at, cells_in, mean, Scorecard};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::OnceLock;

/// PB2002's poles are relative to the Pacific plate; adding the Pacific's own motion over
/// the hotspots (HS3-NUVEL1A, Gripp & Gordon 2002) gives every plate its motion over the
/// mantle, which hotspot tracks need.
const PACIFIC_OVER_HOTSPOTS: (f64, f64, f64) = (-61.467, 90.326, 1.0613);

fn rotation(lat: f64, lon: f64, deg_per_myr: f64) -> Vec3 {
    UnitVec3::from_lat_lon(lat, lon).v() * deg_per_myr
}

struct Earth {
    grid: geoforge_core::grid::AnyGrid,
    products: Products,
    /// PB2002 plate code per partition key.
    codes: BTreeMap<u32, String>,
    etopo: Vec<f32>,
}

fn earth() -> &'static Earth {
    static E: OnceLock<Earth> = OnceLock::new();
    E.get_or_init(|| {
        let g = known::grid_of(720, 360);
        let idx = known::plate_index(&g);
        let etopo = known::etopo(&g);
        let rows = known::csv("pb2002_plates.csv");
        let (plat, plon, prate) = PACIFIC_OVER_HOTSPOTS;
        let pacific = rotation(plat, plon, prate);
        let mut table = Table::new();
        let mut codes = BTreeMap::new();
        for (i, r) in rows.iter().enumerate() {
            let w = rotation(
                r[2].parse().unwrap(),
                r[3].parse().unwrap(),
                r[4].parse().unwrap(),
            ) + pacific;
            let cells: Vec<UnitVec3> = g
                .cells()
                .filter(|&c| usize::from(idx[c as usize]) == i)
                .map(|c| g.center(c))
                .collect();
            let Some(home) = UnitVec3::centroid(cells) else {
                continue;
            };
            let pole = w.normalized().unwrap();
            table.insert(
                Key(i as u64 + 1),
                Plate {
                    name: r[1].clone(),
                    euler_pole: pole,
                    rate_deg_per_myr: w.norm(),
                    seed_point: home,
                },
            );
            codes.insert(i as u32 + 1, r[0].clone());
        }
        let owner = Field::from_values(g.spec(), idx.iter().map(|&i| u32::from(i) + 1).collect());
        let crust = etopo
            .iter()
            .map(|&e| {
                if e > -2000.0 {
                    CrustType::Continental
                } else {
                    CrustType::Oceanic
                }
            })
            .collect();
        let mut inputs = Products::new();
        inputs.insert(ProductData::Plates(Partition { owner, table }));
        inputs.insert(ProductData::CrustType(Field::from_values(g.spec(), crust)));
        let body = Archetype::Earth.context(Key(1), "Earth");
        let products = known::run(
            &body,
            &g,
            inputs,
            &[LayerId::Tectonics, LayerId::Geology, LayerId::Elevation],
        );
        Earth {
            grid: g,
            products,
            codes,
            etopo,
        }
    })
}

/// A PB2002 boundary step: midpoint, class, the two plates, length.
struct Step {
    lat: f64,
    lon: f64,
    class: String,
    plates: (String, String),
    km: f64,
}

fn steps() -> Vec<Step> {
    known::csv("pb2002_steps.csv")
        .into_iter()
        .map(|r| {
            let (a, b) = r[3].split_at(2);
            Step {
                lat: r[0].parse().unwrap(),
                lon: r[1].parse().unwrap(),
                class: r[2].clone(),
                plates: (a.into(), b[1..].into()),
                km: r[4].parse().unwrap(),
            }
        })
        .collect()
}

/// What the model calls each PB2002 class.
fn expected_kind(class: &str) -> BoundaryKind {
    match class {
        "SUB" | "OCB" | "CCB" => BoundaryKind::Convergent,
        "OSR" | "CRB" => BoundaryKind::Divergent,
        _ => BoundaryKind::Transform,
    }
}

/// Cells within `km` of a place.
fn around(e: &Earth, lat: f64, lon: f64, km: f64) -> Vec<Cell> {
    let here = UnitVec3::from_lat_lon(lat, lon);
    let step = 0.25;
    let n = (km / 111.2 / step).ceil() as i32 + 1;
    let mut v = BTreeSet::new();
    for dy in -n..=n {
        for dx in -n..=n {
            let la = (lat + f64::from(dy) * step).clamp(-89.9, 89.9);
            let lo = lon + f64::from(dx) * step / la.to_radians().cos().max(0.05);
            let c = at(&e.grid, la, lo);
            if e.grid.center(c).angle_to(here) * 6371.0 <= km {
                v.insert(c);
            }
        }
    }
    v.insert(at(&e.grid, lat, lon));
    v.into_iter().collect()
}

fn province_kinds_near(e: &Earth, lat: f64, lon: f64, km: f64) -> BTreeSet<G> {
    let p = e.products.provinces().unwrap();
    around(e, lat, lon, km)
        .into_iter()
        .filter_map(|c| p.row_at(c).map(|row| row.kind))
        .collect()
}

fn province_at(e: &Earth, lat: f64, lon: f64) -> String {
    let p = e.products.provinces().unwrap();
    p.row_at(at(&e.grid, lat, lon))
        .map_or("none".into(), |row| format!("{:?}", row.kind))
}

/// The model's boundary type per (cell, plate code, plate code), both orders.
fn boundary_kinds(e: &Earth) -> BTreeMap<(Cell, String, String), BoundaryKind> {
    let plates = e.products.plates().unwrap();
    let segs = boundaries::analyze(
        &e.grid,
        plates,
        6371.0,
        &geoforge_surface::tectonics::TectonicsParams::default().boundaries,
    );
    let mut kind_at = BTreeMap::new();
    for seg in &segs {
        let (a, b) = (&e.codes[&(seg.a.0 as u32)], &e.codes[&(seg.b.0 as u32)]);
        for &c in &seg.cells {
            kind_at.insert((c, a.clone(), b.clone()), seg.kind);
            kind_at.insert((c, b.clone(), a.clone()), seg.kind);
        }
    }
    kind_at
}

/// Each PB2002 step against the model's boundary type there, by length.
fn boundary_types(
    s: &mut Scorecard,
    e: &Earth,
    steps: &[Step],
    kind_at: &BTreeMap<(Cell, String, String), BoundaryKind>,
) {
    let mut by_class: BTreeMap<&str, (f64, f64, f64)> = BTreeMap::new(); // (agree, found, total) km
    for st in steps {
        let found = around(e, st.lat, st.lon, 80.0).into_iter().find_map(|c| {
            kind_at
                .get(&(c, st.plates.0.clone(), st.plates.1.clone()))
                .copied()
        });
        let t = by_class.entry(&st.class).or_default();
        t.2 += st.km;
        if let Some(k) = found {
            t.1 += st.km;
            if k == expected_kind(&st.class) {
                t.0 += st.km;
            }
        }
    }
    let (agree, found, total) = by_class
        .values()
        .fold((0.0, 0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1, a.2 + b.2));
    s.range(
        "boundaries",
        "PB2002 boundary length found on the grid",
        found / total,
        0.8..=1.0,
        "",
    );
    s.range(
        "boundaries",
        "type agrees with PB2002 (of length found)",
        agree / found,
        0.7..=1.0,
        "",
    );
    for (class, what) in [
        ("SUB", "subduction → convergent"),
        ("OSR", "spreading ridge → divergent"),
        ("OTF", "oceanic transform → transform"),
        ("CCB", "continental collision → convergent"),
        ("CRB", "continental rift → divergent"),
    ] {
        let (a, f, _) = by_class[class];
        s.range(
            "boundaries",
            format!("{class}: {what}"),
            a / f,
            0.6..=1.0,
            "",
        );
        if class == "OTF" {
            s.gap("ridge–transform staircases read as oblique ridge: the boundary normal is fitted over 150 km and stretches under 100 km merge");
        }
    }
}

fn crust_age(s: &mut Scorecard, e: &Earth, steps: &[Step]) {
    let g = &e.grid;
    let age = known::scalar(&e.products, Products::crust_age);
    let crust = e.products.crust_type().unwrap();
    let ocean: Vec<Cell> = g
        .cells()
        .filter(|&c| *crust.get(c) == CrustType::Oceanic)
        .collect();
    s.range(
        "crust age",
        "mean ocean-floor age",
        mean(g, &ocean, |c| age[c as usize]),
        45.0..=85.0,
        "Ma",
    );
    s.range(
        "crust age",
        "oldest ocean floor",
        ocean.iter().map(|&c| age[c as usize]).fold(0.0, f64::max),
        150.0..=280.0,
        "Ma",
    );
    let mut ridge_ages: Vec<f64> = steps
        .iter()
        .filter(|st| st.class == "OSR")
        .map(|st| at(g, st.lat, st.lon))
        .filter(|&c| *crust.get(c) == CrustType::Oceanic)
        .map(|c| age[c as usize])
        .collect();
    ridge_ages.sort_by(f64::total_cmp);
    s.range(
        "crust age",
        "median age at PB2002 spreading ridges",
        ridge_ages[ridge_ages.len() / 2],
        0.0..=10.0,
        "Ma",
    );
}

/// Provinces at places with a known geology.
fn provinces(s: &mut Scorecard, e: &Earth, steps: &[Step]) {
    let g = &e.grid;
    let near_share = |class: &str, kind: &[G], km: f64| {
        let st: Vec<&Step> = steps.iter().filter(|st| st.class == class).collect();
        let hit = st
            .iter()
            .filter(|st| {
                kind.iter()
                    .any(|k| province_kinds_near(e, st.lat, st.lon, km).contains(k))
            })
            .map(|st| st.km)
            .sum::<f64>();
        hit / st.iter().map(|st| st.km).sum::<f64>()
    };
    s.range(
        "provinces",
        "spreading ridges with a mid-ocean ridge nearby",
        near_share("OSR", &[G::MidOceanRidge], 100.0),
        0.8..=1.0,
        "",
    );
    s.range(
        "provinces",
        "subduction zones with a trench nearby",
        near_share("SUB", &[G::OceanTrench], 150.0),
        0.7..=1.0,
        "",
    );
    s.range(
        "provinces",
        "subduction zones with a volcanic arc within 350 km",
        near_share("SUB", &[G::VolcanicArc], 350.0),
        0.6..=1.0,
        "",
    );
    s.range(
        "provinces",
        "continental collisions with an orogen nearby",
        near_share("CCB", &[G::CollisionOrogen], 200.0),
        0.6..=1.0,
        "",
    );
    s.range(
        "provinces",
        "continental rifts with a rift nearby",
        near_share("CRB", &[G::ContinentalRift], 200.0),
        0.5..=1.0,
        "",
    );
    // Arc–trench gap: from each subduction step to the nearest arc cell within 800 km.
    let prov = e.products.provinces().unwrap();
    let mut gaps_km: Vec<f64> = steps
        .iter()
        .filter(|st| st.class == "SUB")
        .filter_map(|st| {
            let here = UnitVec3::from_lat_lon(st.lat, st.lon);
            around(e, st.lat, st.lon, 800.0)
                .into_iter()
                .filter(|&c| prov.row_at(c).is_some_and(|r| r.kind == G::VolcanicArc))
                .map(|c| g.center(c).angle_to(here) * 6371.0)
                .reduce(f64::min)
        })
        .collect();
    gaps_km.sort_by(f64::total_cmp);
    s.range(
        "provinces",
        "median arc–trench distance (Earth 150–300)",
        gaps_km[gaps_km.len() / 2],
        120.0..=320.0,
        "km",
    );
    let stable = ["Craton", "Platform", "IntracratonicBasin", "PaleoOrogen"];
    for (name, lat, lon) in [
        ("Canadian Shield", 58.0, -95.0),
        ("Siberian platform", 63.0, 105.0),
        ("West Africa", 15.0, -5.0),
        ("Western Australia", -26.0, 122.0),
        ("Amazon craton", -8.0, -55.0),
    ] {
        s.one_of(
            "provinces",
            format!("{name} is a stable interior"),
            &province_at(e, lat, lon),
            &stable,
        );
    }
    let arc = [
        "VolcanicArc",
        "ForearcBasin",
        "BackarcBasin",
        "AccretionaryWedge",
        "CollisionOrogen",
        "OrogenicPlateau",
    ];
    s.one_of("provinces", "central Andes (20°S 68°W)", &province_at(e, -20.0, -68.0), &arc);
    s.one_of(
        "provinces",
        "Himalaya (29°N 84°E)",
        &province_at(e, 29.0, 84.0),
        &["CollisionOrogen"],
    );
    s.one_of(
        "provinces",
        "Zagros (32°N 50°E)",
        &province_at(e, 32.0, 50.0),
        &["CollisionOrogen"],
    );
    // India and Arabia go under Asia: the plateau grows on Asia, the foreland on them.
    for (name, lat, lon, want) in [
        ("Tibet (33°N 88°E)", 33.0, 88.0, "OrogenicPlateau"),
        ("Ganges plain (26°N 82°E)", 26.0, 82.0, "ForelandBasin"),
        ("Mesopotamia (32°N 45°E)", 32.0, 45.0, "ForelandBasin"),
    ] {
        s.one_of("provinces", name, &province_at(e, lat, lon), &[want]);
    }
}

/// The heightmap against ETOPO1.
fn elevation(s: &mut Scorecard, e: &Earth, steps: &[Step]) {
    let g = &e.grid;
    let z = known::scalar(&e.products, Products::elevation);
    let real: Vec<f64> = e.etopo.iter().map(|&v| f64::from(v)).collect();
    let all: Vec<Cell> = g.cells().collect();
    let land_share = |v: &[f64]| mean(g, &all, |c| f64::from(u8::from(v[c as usize] >= 0.0)));
    s.range(
        "elevation",
        format!(
            "land share of the surface (ETOPO1 {:.2})",
            land_share(&real)
        ),
        land_share(&z),
        0.25..=0.33,
        "",
    );
    let corr = |cells: &[Cell]| {
        let (mx, my) = (
            mean(g, cells, |c| z[c as usize]),
            mean(g, cells, |c| real[c as usize]),
        );
        let cov = mean(g, cells, |c| (z[c as usize] - mx) * (real[c as usize] - my));
        let sx = mean(g, cells, |c| (z[c as usize] - mx).powi(2)).sqrt();
        let sy = mean(g, cells, |c| (real[c as usize] - my).powi(2)).sqrt();
        cov / (sx * sy)
    };
    let mid: Vec<Cell> = all
        .iter()
        .copied()
        .filter(|&c| known::lat_lon(g, c).0.abs() < 60.0)
        .collect();
    let land: Vec<Cell> = mid
        .iter()
        .copied()
        .filter(|&c| real[c as usize] >= 0.0)
        .collect();
    let deep: Vec<Cell> = mid
        .iter()
        .copied()
        .filter(|&c| real[c as usize] < -2000.0)
        .collect();
    s.range(
        "elevation",
        "correlation with ETOPO1, 60°S–60°N",
        corr(&mid),
        0.85..=1.0,
        "",
    );
    s.range(
        "elevation",
        "correlation over land",
        corr(&land),
        0.4..=1.0,
        "",
    );
    s.range(
        "elevation",
        "correlation over the deep ocean",
        corr(&deep),
        0.5..=1.0,
        "",
    );

    // Each statistic for the model and for ETOPO1.
    let both = |f: &dyn Fn(&[f64]) -> f64| (f(&z), f(&real));
    let where_mean = |v: &[f64], cells: &[Cell], keep: &dyn Fn(f64) -> bool| {
        mean(
            g,
            &cells
                .iter()
                .copied()
                .filter(|&c| keep(v[c as usize]))
                .collect::<Vec<_>>(),
            |c| v[c as usize],
        )
    };
    let (m, r) = both(&|v| where_mean(v, &all, &|x| x < 0.0));
    s.range(
        "elevation",
        format!("mean ocean depth (ETOPO1 {r:.0})"),
        m,
        -4200.0..=-3200.0,
        "m",
    );
    let (m, r) = both(&|v| where_mean(v, &mid, &|x| x >= 0.0));
    s.range(
        "elevation",
        format!("mean land height, 60°S–60°N (ETOPO1 {r:.0})"),
        m,
        450.0..=950.0,
        "m",
    );
    let (m, r) = both(&|v| mean(g, &land, |c| f64::from(u8::from(v[c as usize] > 2000.0))));
    s.range(
        "elevation",
        format!("land above 2000 m, 60°S–60°N (ETOPO1 {r:.2})"),
        m,
        0.03..=0.10,
        "",
    );
    let ridge: Vec<Cell> = steps
        .iter()
        .filter(|st| st.class == "OSR")
        .map(|st| at(g, st.lat, st.lon))
        .collect();
    let (m, r) = both(&|v| mean(g, &ridge, |c| v[c as usize]));
    s.range(
        "elevation",
        format!("depth at spreading ridges (ETOPO1 {r:.0})"),
        m,
        -3300.0..=-2200.0,
        "m",
    );
    let subs: Vec<&Step> = steps
        .iter()
        .filter(|st| st.class == "SUB" && st.lat.abs() < 60.0)
        .collect();
    let trench = |v: &[f64]| {
        subs.iter()
            .map(|st| {
                around(e, st.lat, st.lon, 200.0)
                    .into_iter()
                    .map(|c| v[c as usize])
                    .fold(f64::INFINITY, f64::min)
            })
            .sum::<f64>()
            / subs.len() as f64
    };
    let (m, r) = both(&trench);
    s.range(
        "elevation",
        format!("deepest point near subduction zones (ETOPO1 {r:.0})"),
        m,
        -9000.0..=-4500.0,
        "m",
    );
    for (name, lat, lon0, lon1, want, gap) in [
        ("Tibet, 30–35°N 80–95°E", 30.0..=35.0, 80.0, 95.0, 3500.0..=5500.0, None),
        ("Himalaya, 27.5–29.5°N 80–90°E", 27.5..=29.5, 80.0, 90.0, 2500.0..=5500.0, None),
        ("Ganges plain, 24–27°N 78–88°E", 24.0..=27.0, 78.0, 88.0, 0.0..=500.0, None),
        ("Iranian plateau, 30–35°N 52–60°E", 30.0..=35.0, 52.0, 60.0, 800.0..=2200.0, None),
        ("Mesopotamia, 30–34°N 43–47°E", 30.0..=34.0, 43.0, 47.0, -100.0..=400.0, None),
        (
            "central Andes, 15–25°S 70–65°W",
            -25.0..=-15.0,
            -70.0,
            -65.0,
            2000.0..=4000.0,
            None,
        ),
        (
            "Peruvian Andes, 8–12°S 78–74°W",
            -12.0..=-8.0,
            -78.0,
            -74.0,
            1200.0..=3000.0,
            Some("over Peru's flat slab the cordillera is ~300 km of 4 km highlands; the model's young-slab cordillera is ~150 km wide and its forearc erodes to sea level"),
        ),
        (
            "Pampean ranges, 28–34°S 68–64°W",
            -34.0..=-28.0,
            -68.0,
            -64.0,
            200.0..=1500.0,
            Some("no flat-slab segments: the model's slab is as old here as under the Altiplano, so the plateau runs on to 37°S"),
        ),
        (
            "Cascades, 43–48°N 123–120°W",
            43.0..=48.0,
            -123.0,
            -120.0,
            400.0..=1500.0,
            Some("a slow margin's cordillera stands only one or two cells wide, and erodes"),
        ),
        ("Amazon basin, 5°S–0 65–55°W", -5.0..=0.0, -65.0, -55.0, 0.0..=300.0, Some("the retroarc foreland ends ~300 km past the Andes; the lowland beyond it sits at the continental base, ~300 m")),
        ("Rocky Mountains, 38–42°N 110–105°W", 38.0..=42.0, -110.0, -105.0, 1500.0..=3000.0, Some("no mountain building far inland of a plate boundary (the Laramide uplift, Tian Shan)")),
        ("Ethiopian highlands, 5–10°N 36–40°E", 5.0..=10.0, 36.0, 40.0, 1000.0..=2600.0, Some("no plume uplift: the Ethiopian and East African domes are 1.5–2 km over 1000 km")),
    ] {
        let cells = cells_in(g, lat, lon0, lon1);
        let (m, r) = both(&|v| mean(g, &cells, |c| v[c as usize]));
        s.range("elevation", format!("{name} (ETOPO1 {r:.0})"), m, want, "m");
        if let Some(why) = gap {
            s.gap(why);
        }
    }
}

/// With `KNOWN_CASES_DUMP` set, the fields for a look at the map.
fn dump(e: &Earth, kind_at: &BTreeMap<(Cell, String, String), BoundaryKind>) {
    let prov = e.products.provinces().unwrap();
    known::dump(
        "tectonics_elevation",
        known::scalar(&e.products, Products::elevation),
    );
    known::dump(
        "tectonics_crust_age",
        known::scalar(&e.products, Products::crust_age),
    );
    known::dump(
        "tectonics_provinces",
        e.grid
            .cells()
            .map(|c| prov.row_at(c).map_or(-1.0, |r| f64::from(r.kind as u8))),
    );
    let mut kinds = vec![-1.0; e.grid.len()];
    for ((c, _, _), k) in kind_at {
        kinds[*c as usize] = f64::from(*k as u8);
    }
    known::dump("tectonics_boundary_kind", kinds);
}

fn scorecard() -> &'static Scorecard {
    static S: OnceLock<Scorecard> = OnceLock::new();
    S.get_or_init(|| {
        let e = earth();
        let steps = steps();
        let kind_at = boundary_kinds(e);
        dump(e, &kind_at);
        let mut s = Scorecard::default();
        boundary_types(&mut s, e, &steps, &kind_at);
        crust_age(&mut s, e, &steps);
        provinces(&mut s, e, &steps);
        elevation(&mut s, e, &steps);
        s
    })
}

#[test]
fn tectonics_known_cases() {
    let s = scorecard();
    println!("{}", s.report());
    s.assert_cases("earth_tectonics");
}

#[test]
#[ignore = "known gaps: run with --ignored"]
fn tectonics_known_gaps() {
    scorecard().assert_gaps();
}
