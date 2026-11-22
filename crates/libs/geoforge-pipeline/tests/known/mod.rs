//! Shared pieces of the known-case tests: fixtures on a one-degree grid, layers run in
//! memory, and a scorecard.
//!
//! A known case feeds a layer real inputs and checks its output against what the real body
//! does. Cases the model is known to miss are recorded as gaps: they never fail the main
//! tests, and each file's `#[ignore]`d `…_known_gaps` test asserts them, so
//! `cargo test -- --ignored` lists what is still wrong.
#![allow(dead_code)]

use geoforge_core::grid::{AnyGrid, GridSpec};
use geoforge_core::product::{LayerId, ProductData, Products, RunContext};
use geoforge_core::{BodyContext, Cell, Field, Grid, Seed, UnitVec3};
use geoforge_pipeline::registry;
use std::fmt::Write;
use std::ops::RangeInclusive;
use std::path::PathBuf;

/// The one-degree Earth the climate cases run on.
pub fn grid() -> AnyGrid {
    grid_of(360, 180)
}

pub fn grid_of(width: u32, height: u32) -> AnyGrid {
    GridSpec::Equirect { width, height }.build().unwrap()
}

pub fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures").join(name)
}

/// An equirectangular image's value at (lat, lon), bilinear between pixel centres.
fn sample(px: &[f32], w: usize, h: usize, lat: f64, lon: f64) -> f32 {
    let fy = ((90.0 - lat) * h as f64 / 180.0 - 0.5).clamp(0.0, h as f64 - 1.0);
    let fx = (lon + 180.0) * w as f64 / 360.0 - 0.5;
    let (y0, x0) = (fy.floor() as usize, fx.floor());
    let (ty, tx) = (fy - y0 as f64, fx - x0);
    let y1 = (y0 + 1).min(h - 1);
    let col = |x: f64| (x as i64).rem_euclid(w as i64) as usize;
    let (xa, xb) = (col(x0), col(x0 + 1.0));
    let v = |y: usize, x: usize| f64::from(px[y * w + x]);
    let top = v(y0, xa) * (1.0 - tx) + v(y0, xb) * tx;
    let bot = v(y1, xa) * (1.0 - tx) + v(y1, xb) * tx;
    (top * (1.0 - ty) + bot * ty) as f32
}

/// ETOPO1 one-degree block means (see `scripts/known-cases/make_fixtures.py`), m, sampled
/// at each cell centre of `g`.
pub fn etopo(g: &AnyGrid) -> Vec<f32> {
    let img = image::open(fixture("earth_etopo1_360x180.png")).unwrap().to_luma16();
    let px: Vec<f32> = img.pixels().map(|p| f32::from(p.0[0]) - 12000.0).collect();
    g.cells()
        .map(|c| {
            let (lat, lon) = lat_lon(g, c);
            sample(&px, 360, 180, lat, lon)
        })
        .collect()
}

/// PB2002 plate index per cell of a 720×360 grid.
pub fn plate_index(g: &AnyGrid) -> Vec<u8> {
    let img = image::open(fixture("pb2002_plates_720x360.png")).unwrap().to_luma8();
    let (w, h) = img.dimensions();
    assert_eq!(g.len(), (w * h) as usize, "plate fixture is {w}×{h}");
    let mut out = vec![0; g.len()];
    for (x, y, p) in img.enumerate_pixels() {
        let lat = 90.0 - (f64::from(y) + 0.5) * 180.0 / f64::from(h);
        let lon = -180.0 + (f64::from(x) + 0.5) * 360.0 / f64::from(w);
        out[g.locate(UnitVec3::from_lat_lon(lat, lon)) as usize] = p.0[0];
    }
    out
}

pub fn csv(name: &str) -> Vec<Vec<String>> {
    std::fs::read_to_string(fixture(name)).unwrap().lines().skip(1).map(|l| l.split(',').map(str::to_string).collect()).collect()
}

pub fn at(g: &AnyGrid, lat: f64, lon: f64) -> Cell {
    g.locate(UnitVec3::from_lat_lon(lat, lon))
}

pub fn lat_lon(g: &AnyGrid, c: Cell) -> (f64, f64) {
    let ll = g.center(c).to_lat_lon();
    (ll.lat_deg, ll.lon_deg)
}

/// Cells in a lat/lon box (longitudes may wrap: `lon0 > lon1` crosses the dateline).
pub fn cells_in(g: &AnyGrid, lat: RangeInclusive<f64>, lon0: f64, lon1: f64) -> Vec<Cell> {
    g.cells()
        .filter(|&c| {
            let (la, lo) = lat_lon(g, c);
            lat.contains(&la) && if lon0 <= lon1 { (lon0..=lon1).contains(&lo) } else { lo >= lon0 || lo <= lon1 }
        })
        .collect()
}

/// Area-weighted mean of `f` over `cells`.
pub fn mean(g: &AnyGrid, cells: &[Cell], f: impl Fn(Cell) -> f64) -> f64 {
    let (s, w) = cells.iter().fold((0.0, 0.0), |(s, w), &c| (s + f(c) * g.area_sr(c), w + g.area_sr(c)));
    s / w
}

/// The cell at a place, or the nearest cell within `r` cells that passes `ok` (a coastal
/// station whose cell is sea takes its nearest land cell).
pub fn near(g: &AnyGrid, lat: f64, lon: f64, r: i32, ok: impl Fn(Cell) -> bool) -> Cell {
    let c = at(g, lat, lon);
    if ok(c) {
        return c;
    }
    let mut best: Option<(f64, Cell)> = None;
    for dy in -r..=r {
        for dx in -r..=r {
            let la = (lat + f64::from(dy)).clamp(-89.5, 89.5);
            let n = at(g, la, lon + f64::from(dx));
            let d = f64::from(dx * dx + dy * dy);
            if ok(n) && best.is_none_or(|(bd, _)| d < bd) {
                best = Some((d, n));
            }
        }
    }
    best.map_or(c, |(_, n)| n)
}

/// Run `layers` in order on `inputs`, each by the generator its body calls for, with
/// default parameters. Returns inputs and outputs together.
pub fn run(body: &BodyContext, g: &AnyGrid, products: Products, layers: &[LayerId]) -> Products {
    run_seeded(1, body, g, products, layers)
}

/// [`run`] with another world seed.
pub fn run_seeded(seed: u64, body: &BodyContext, g: &AnyGrid, mut products: Products, layers: &[LayerId]) -> Products {
    for &l in layers {
        let i = registry::select(l, body).unwrap_or_else(|e| panic!("{}: {e:?}", l.name()));
        let ctx = RunContext { root_seed: Seed(seed), body, grid: g };
        let out = i.run(&ctx, &products, &i.default_params()).unwrap();
        products.extend(out);
    }
    products
}

pub fn scalar(p: &Products, d: fn(&Products) -> Option<&Field<f32>>) -> Vec<f64> {
    d(p).unwrap().values.iter().map(|&v| f64::from(v)).collect()
}

/// With `KNOWN_CASES_DUMP` set, write a field as little-endian f32 in grid-cell order to
/// `target/known-cases/<name>.f32`, for plotting.
pub fn dump(name: &str, values: impl IntoIterator<Item = f64>) {
    if std::env::var_os("KNOWN_CASES_DUMP").is_none() {
        return;
    }
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../../target/known-cases");
    let _ = std::fs::create_dir_all(&dir);
    let bytes: Vec<u8> = values.into_iter().flat_map(|v| (v as f32).to_le_bytes()).collect();
    let _ = std::fs::write(dir.join(format!("{name}.f32")), bytes);
}

pub fn elevation_product(g: &AnyGrid, e: Vec<f32>) -> ProductData {
    ProductData::Elevation(Field::from_values(g.spec(), e))
}

/// A real heightmap given to the pipeline as its rough heightmap too, so the climate sees it
/// and erosion leaves it as it is.
pub fn rough_elevation_product(g: &AnyGrid, e: Vec<f32>) -> ProductData {
    ProductData::RoughElevation(Field::from_values(g.spec(), e))
}

/// One known case.
pub struct Check {
    pub layer: &'static str,
    pub name: String,
    pub value: String,
    pub expect: String,
    pub pass: bool,
    /// Why the model misses it, for a known gap.
    pub gap: Option<&'static str>,
}

#[derive(Default)]
pub struct Scorecard {
    pub checks: Vec<Check>,
}

impl Scorecard {
    /// A number that should fall in `range`.
    pub fn range(&mut self, layer: &'static str, name: impl Into<String>, value: f64, range: RangeInclusive<f64>, unit: &str) -> &mut Self {
        let fmt = |v: f64| {
            if v.abs() >= 100.0 {
                format!("{v:.0}")
            } else {
                format!("{v:.2}")
            }
        };
        self.checks.push(Check {
            layer,
            name: name.into(),
            value: format!("{} {unit}", fmt(value)),
            expect: format!("{}…{} {unit}", fmt(*range.start()), fmt(*range.end())),
            pass: range.contains(&value),
            gap: None,
        });
        self
    }

    /// A category that should be one of `want`.
    pub fn one_of(&mut self, layer: &'static str, name: impl Into<String>, value: &str, want: &[&str]) -> &mut Self {
        self.checks.push(Check { layer, name: name.into(), value: value.into(), expect: want.join(" or "), pass: want.contains(&value), gap: None });
        self
    }

    /// Mark the last case as a known gap.
    pub fn gap(&mut self, why: &'static str) -> &mut Self {
        self.checks.last_mut().unwrap().gap = Some(why);
        self
    }

    /// Mark the last case a known gap if `key` is listed in `gaps`.
    pub fn gap_if_listed(&mut self, key: &str, gaps: &[(&str, &'static str)]) -> &mut Self {
        if let Some(&(_, why)) = gaps.iter().find(|(k, _)| *k == key) {
            self.gap(why);
        }
        self
    }

    pub fn report(&self) -> String {
        let mut s = String::new();
        for c in &self.checks {
            let mark = match (c.pass, c.gap.is_some()) {
                (true, false) => "pass",
                (false, false) => "FAIL",
                (false, true) => "gap ",
                (true, true) => "gap passes now",
            };
            let _ = writeln!(s, "{mark:5} {:11} {:52} {:>18}  expect {}", c.layer, c.name, c.value, c.expect);
            if let (Some(why), false) = (c.gap, c.pass) {
                let _ = writeln!(s, "{:18}{why}", "");
            }
        }
        let pass = self.checks.iter().filter(|c| c.pass).count();
        let _ = writeln!(s, "{pass} of {} known cases pass; {} are known gaps", self.checks.len(), self.checks.iter().filter(|c| c.gap.is_some()).count());
        s
    }

    /// Every case that is not a known gap passes. Writes the scorecard to
    /// `target/known-cases/<name>.txt` (and `.json`) either way.
    pub fn assert_cases(&self, name: &str) {
        let report = self.report();
        let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../../target/known-cases");
        let _ = std::fs::create_dir_all(&dir);
        let _ = std::fs::write(dir.join(format!("{name}.txt")), &report);
        let rows: Vec<serde_json::Value> =
            self.checks.iter().map(|c| serde_json::json!({"layer": c.layer, "name": c.name, "value": c.value, "expect": c.expect, "pass": c.pass, "gap": c.gap})).collect();
        let _ = std::fs::write(dir.join(format!("{name}.json")), serde_json::to_string_pretty(&rows).unwrap_or_default());
        let failed: Vec<&Check> = self.checks.iter().filter(|c| c.gap.is_none() && !c.pass).collect();
        assert!(failed.is_empty(), "{} known cases fail:\n{report}", failed.len());
    }

    /// Every known gap passes (expected to fail until the model improves).
    pub fn assert_gaps(&self) {
        let open: Vec<&Check> = self.checks.iter().filter(|c| c.gap.is_some() && !c.pass).collect();
        assert!(open.is_empty(), "{} known gaps are open:\n{}", open.len(), self.report());
    }
}
