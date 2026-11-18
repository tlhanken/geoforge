//! Preview renders (RGB PNG) of products, for looking at a world rather than editing it.
//! Renders use other products for context when present: hillshade from elevation under
//! biomes, rivers and lakes drawn over the land, plate outlines over plates.

use crate::edit::{hsv, key_color};
use crate::{Result, StoreError};
use geoforge_core::grid::{AnyGrid, Grid, GridSpec};
use geoforge_core::product::{ProductData, ProductId, Products};
use geoforge_core::surface::Category;
use image::{Rgb, RgbImage};

/// Hypsometric tint (sea level = 0 m).
#[must_use]
pub fn elevation_color(e: f64) -> [u8; 3] {
    const SEA: [(f64, [u8; 3]); 6] = [
        (-11000.0, [8, 14, 48]),
        (-6000.0, [16, 32, 92]),
        (-4000.0, [28, 62, 132]),
        (-2000.0, [48, 98, 170]),
        (-200.0, [92, 150, 210]),
        (0.0, [150, 200, 235]),
    ];
    const LAND: [(f64, [u8; 3]); 8] = [
        (0.0, [70, 120, 60]),
        (200.0, [110, 150, 80]),
        (500.0, [170, 180, 105]),
        (1000.0, [200, 180, 120]),
        (2000.0, [170, 130, 90]),
        (3000.0, [140, 105, 85]),
        (4500.0, [190, 185, 180]),
        (6500.0, [250, 250, 250]),
    ];
    ramp(if e < 0.0 { &SEA } else { &LAND }, e)
}

fn ramp(stops: &[(f64, [u8; 3])], x: f64) -> [u8; 3] {
    if x <= stops[0].0 {
        return stops[0].1;
    }
    for w in stops.windows(2) {
        let ((a, ca), (b, cb)) = (w[0], w[1]);
        if x <= b {
            let t = (x - a) / (b - a);
            return [0, 1, 2].map(|i| (f64::from(ca[i]) + t * (f64::from(cb[i]) - f64::from(ca[i]))).round() as u8);
        }
    }
    stops[stops.len() - 1].1
}

fn temperature_color(t: f64) -> [u8; 3] {
    ramp(
        &[
            (-45.0, [60, 0, 110]),
            (-25.0, [40, 70, 200]),
            (-10.0, [90, 160, 230]),
            (0.0, [210, 235, 250]),
            (10.0, [250, 240, 160]),
            (20.0, [250, 170, 60]),
            (28.0, [220, 70, 30]),
            (35.0, [130, 0, 20]),
        ],
        t,
    )
}

fn precip_color(p: f64) -> [u8; 3] {
    ramp(
        &[(0.0, [150, 90, 40]), (200.0, [220, 190, 120]), (500.0, [235, 235, 170]), (1000.0, [120, 200, 110]), (2000.0, [30, 140, 90]), (3500.0, [20, 90, 160]), (6000.0, [40, 20, 120])],
        p,
    )
}

fn age_color(a: f64) -> [u8; 3] {
    ramp(&[(0.0, [220, 40, 40]), (20.0, [250, 160, 40]), (60.0, [240, 230, 80]), (120.0, [80, 190, 110]), (180.0, [40, 110, 190]), (400.0, [70, 40, 130]), (2500.0, [40, 30, 60])], a)
}

fn shade(c: [u8; 3], k: f64) -> [u8; 3] {
    c.map(|v| (f64::from(v) * k).clamp(0.0, 255.0) as u8)
}

fn blend(a: [u8; 3], b: [u8; 3], t: f64) -> [u8; 3] {
    [0, 1, 2].map(|i| (f64::from(a[i]) * (1.0 - t) + f64::from(b[i]) * t).round() as u8)
}

/// Hillshade factor per cell (≈0.6–1.25), light from the north-west.
fn hillshade(grid: &AnyGrid, elev: &[f32], radius_km: f64, exaggeration: f64) -> Vec<f64> {
    let Some(eq) = grid.as_equirect() else { return vec![1.0; elev.len()] };
    let (w, h) = (eq.width(), eq.height());
    let dy_m = eq.lat_res_deg().to_radians() * radius_km * 1000.0;
    let mut out = vec![1.0; elev.len()];
    for y in 0..h {
        let (lat, _) = eq.lat_lon_xy(0, y);
        let dx_m = (eq.lon_res_deg().to_radians() * radius_km * 1000.0 * lat.to_radians().cos()).max(dy_m * 0.05);
        for x in 0..w {
            let at = |xx: i64, yy: i64| {
                let yy = yy.clamp(0, i64::from(h) - 1) as u32;
                let xx = xx.rem_euclid(i64::from(w)) as u32;
                f64::from(elev[eq.cell(xx, yy) as usize].max(0.0))
            };
            let (xi, yi) = (i64::from(x), i64::from(y));
            let gx = (at(xi + 1, yi) - at(xi - 1, yi)) / (2.0 * dx_m) * exaggeration;
            let gy = (at(xi, yi + 1) - at(xi, yi - 1)) / (2.0 * dy_m) * exaggeration;
            // Normal (-gx, gy, 1) with y pointing south; light from NW, 45° up.
            let n = (gx * gx + gy * gy + 1.0).sqrt();
            let (lx, ly, lz) = (-0.5f64, -0.5f64, std::f64::consts::FRAC_1_SQRT_2);
            let d = (-gx * lx + gy * ly + lz) / n;
            out[eq.cell(x, y) as usize] = (0.35 + 0.9 * d).clamp(0.55, 1.25);
        }
    }
    out
}

fn image_of(grid: GridSpec, px: impl Fn(usize) -> [u8; 3]) -> Result<RgbImage> {
    let GridSpec::Equirect { width, height } = grid else {
        return Err(StoreError::Invalid("renders need an equirectangular grid".into()));
    };
    let mut img = RgbImage::new(width, height);
    for (i, p) in img.pixels_mut().enumerate() {
        *p = Rgb(px(i));
    }
    Ok(img)
}

/// Render one product, using `ctx` for context layers. `None` for products with no raster.
///
/// # Errors
/// For grids without a raster form.
pub fn render(p: ProductId, ctx: &Products, grid: &AnyGrid, radius_km: f64) -> Result<Option<RgbImage>> {
    let spec = grid.spec();
    let elev = ctx.elevation().map(|f| f.values.as_slice());
    let shade_f = elev.map(|e| hillshade(grid, e, radius_km, 12.0));
    let lit = |i: usize, c: [u8; 3]| match (&shade_f, elev) {
        (Some(s), Some(e)) if e[i] >= 0.0 => shade(c, s[i]),
        _ => c,
    };
    let Some(data) = ctx.get(p) else { return Ok(None) };
    let mut img = match data {
        ProductData::Elevation(f) => image_of(spec, |i| lit(i, elevation_color(f64::from(f.values[i]))))?,
        ProductData::CrustAge(f) => {
            let ct = ctx.crust_type();
            image_of(spec, |i| {
                let c = age_color(f64::from(f.values[i]));
                match ct {
                    Some(ct) if ct.values[i].to_u8() == 1 => blend(c, [128, 128, 128], 0.55),
                    _ => c,
                }
            })?
        }
        ProductData::Temperature(f) => image_of(spec, |i| lit(i, temperature_color(f64::from(f.values[i]))))?,
        ProductData::TemperatureRange(f) => image_of(spec, |i| lit(i, ramp(&[(0.0, [240, 240, 250]), (10.0, [150, 200, 240]), (25.0, [240, 200, 90]), (45.0, [200, 60, 40]), (70.0, [90, 0, 40])], f64::from(f.values[i]))))?,
        ProductData::Precipitation(f) => image_of(spec, |i| {
            let c = precip_color(f64::from(f.values[i]));
            match elev {
                Some(e) if e[i] < 0.0 => blend(c, [30, 40, 70], 0.55),
                _ => lit(i, c),
            }
        })?,
        ProductData::CrustType(f) => image_of(spec, |i| f.values[i].color())?,
        ProductData::Lithology(f) => image_of(spec, |i| lit(i, f.values[i].color()))?,
        ProductData::RockAge(f) => image_of(spec, |i| {
            let a = f64::from(f.values[i]);
            ramp(&[(0.0, [230, 60, 50]), (66.0, [245, 170, 60]), (252.0, [240, 230, 90]), (541.0, [110, 190, 120]), (1000.0, [60, 150, 200]), (2500.0, [70, 70, 170]), (4000.0, [50, 30, 80])], a)
        })?,
        ProductData::SstAnomaly(f) => image_of(spec, |i| match elev {
            Some(e) if e[i] >= 0.0 => [215, 210, 200],
            _ => ramp(&[(-8.0, [30, 60, 160]), (-3.0, [120, 170, 225]), (0.0, [245, 245, 245]), (3.0, [245, 170, 120]), (8.0, [170, 30, 30])], f64::from(f.values[i])),
        })?,
        ProductData::PrecipitationWarmShare(f) => image_of(spec, |i| match elev {
            Some(e) if e[i] < 0.0 => [40, 50, 75],
            _ => lit(i, ramp(&[(0.15, [180, 90, 40]), (0.35, [235, 200, 120]), (0.5, [240, 240, 230]), (0.65, [140, 200, 150]), (0.85, [40, 110, 160])], f64::from(f.values[i]))),
        })?,
        ProductData::WindJan(f) | ProductData::WindJul(f) | ProductData::OceanCurrents(f) => {
            let ocean = p == ProductId::OceanCurrents;
            let top = if ocean { 1.0 } else { 15.0 };
            let mut img = image_of(spec, |i| {
                let v = f.values[i];
                let sp = f64::from(v[0]).hypot(f64::from(v[1]));
                let base = ramp(&[(0.0, [235, 238, 245]), (top * 0.3, [170, 200, 230]), (top * 0.7, [70, 120, 200]), (top, [40, 40, 140])], sp);
                match elev {
                    Some(e) if e[i] >= 0.0 => if ocean { [200, 195, 180] } else { blend(base, [180, 170, 150], 0.35) },
                    _ => base,
                }
            })?;
            arrows(&mut img, grid, &f.values, top);
            img
        }
        ProductData::Biome(f) => image_of(spec, |i| {
            let c = f.values[i].color();
            lit(i, c)
        })?,
        ProductData::Plates(pt) => image_of(spec, |i| key_color(pt.owner.values[i]))?,
        ProductData::Provinces(pt) => image_of(spec, |i| match pt.owner_of(i as u32).and_then(|k| pt.table.get(k)) {
            Some(prov) => prov.kind.color(),
            None => [30, 30, 30],
        })?,
        ProductData::Watersheds(pt) => image_of(spec, |i| match pt.owner.values[i] {
            0 => [20, 30, 60],
            k => lit(i, key_color(k)),
        })?,
        ProductData::Lakes(_) | ProductData::Rivers(_) => {
            // Drawn as an overlay on the elevation map.
            match elev {
                Some(e) => image_of(spec, |i| {
                    let v = f64::from(e[i]);
                    let c = if v < 0.0 { [190, 205, 220] } else { blend(elevation_color(v), [235, 235, 235], 0.6) };
                    lit(i, c)
                })?,
                None => image_of(spec, |_| [235, 235, 235])?,
            }
        }
        ProductData::TectonicEvents(_) | ProductData::GeologyEvents(_) => return Ok(None),
    };
    // Overlays.
    let plates_outline = matches!(p, ProductId::Plates | ProductId::CrustAge | ProductId::Provinces);
    if plates_outline {
        if let Some(pt) = ctx.plates() {
            for c in grid.cells() {
                let o = pt.owner.values[c as usize];
                if grid.neighbors(c).iter().any(|nb| pt.owner.values[nb as usize] != o && nb > c) {
                    put(&mut img, grid, c, [20, 20, 20]);
                }
            }
        }
    }
    if matches!(p, ProductId::Biome | ProductId::Rivers | ProductId::Lakes | ProductId::Elevation | ProductId::Watersheds | ProductId::Precipitation) {
        if let Some(lk) = ctx.lakes() {
            for (c, &o) in lk.owner.values.iter().enumerate() {
                if o != 0 {
                    let col = match lk.table.get(geoforge_core::key::Key(u64::from(o))) {
                        Some(l) if l.endorheic => [120, 170, 190],
                        _ => [50, 110, 200],
                    };
                    put(&mut img, grid, c as u32, col);
                }
            }
        }
        if let Some(rv) = ctx.rivers() {
            let min_order = if p == ProductId::Rivers { 1 } else { 3 };
            for e in &rv.edges {
                if e.attrs.strahler_order < min_order {
                    continue;
                }
                let t = (f64::from(e.attrs.strahler_order) - f64::from(min_order)) / 6.0;
                let col = blend([110, 160, 230], [10, 40, 150], t.clamp(0.0, 1.0));
                for k in [e.from, e.to] {
                    if let Some(n) = rv.nodes.get(k) {
                        put(&mut img, grid, n.cell, col);
                    }
                }
            }
        }
    }
    Ok(Some(img))
}

/// Sparse direction arrows over a vector field.
fn arrows(img: &mut RgbImage, grid: &AnyGrid, v: &[[f32; 2]], top: f64) {
    let Some(eq) = grid.as_equirect() else { return };
    let (w, h) = (eq.width(), eq.height());
    let step = (w / 48).max(10);
    let mut y = step / 2;
    while y < h {
        let mut x = step / 2;
        while x < w {
            let c = eq.cell(x, y) as usize;
            let (u, n) = (f64::from(v[c][0]), f64::from(v[c][1]));
            let sp = u.hypot(n);
            if sp > top * 0.03 {
                let len = (f64::from(step) * 0.8) * (sp / top).clamp(0.45, 1.0);
                let (dx, dy) = (u / sp * len, -n / sp * len);
                let (x0, y0) = (f64::from(x) - dx / 2.0, f64::from(y) - dy / 2.0);
                let steps = len.ceil() as i32 * 2;
                for k in 0..=steps {
                    let t = f64::from(k) / f64::from(steps);
                    let (px, py) = (x0 + dx * t, y0 + dy * t);
                    if px >= 0.0 && py >= 0.0 && (px as u32) < w && (py as u32) < h {
                        img.put_pixel(px as u32, py as u32, Rgb([20, 20, 30]));
                    }
                }
                // Head: a small dot at the downstream end.
                let (hx, hy) = (x0 + dx, y0 + dy);
                for (ox, oy) in [(0, 0), (1, 0), (0, 1), (-1, 0), (0, -1)] {
                    let (px, py) = (hx as i64 + ox, hy as i64 + oy);
                    if px >= 0 && py >= 0 && (px as u32) < w && (py as u32) < h {
                        img.put_pixel(px as u32, py as u32, Rgb([200, 30, 30]));
                    }
                }
            }
            x += step;
        }
        y += step;
    }
}

fn put(img: &mut RgbImage, grid: &AnyGrid, c: u32, col: [u8; 3]) {
    if let Some(eq) = grid.as_equirect() {
        let (x, y) = eq.xy(c);
        img.put_pixel(x, y, Rgb(col));
    }
}

/// A legend-free colour for an arbitrary index (used by tools).
#[must_use]
pub fn index_color(i: usize) -> [u8; 3] {
    hsv((i as f64 * 0.618_033_988_749_895).fract(), 0.6, 0.9)
}

/// Downscale an image by an integer factor (box filter), for quick looks.
#[must_use]
pub fn downscale(img: &RgbImage, f: u32) -> RgbImage {
    if f <= 1 {
        return img.clone();
    }
    let (w, h) = (img.width() / f, img.height() / f);
    let mut out = RgbImage::new(w, h);
    for y in 0..h {
        for x in 0..w {
            let mut s = [0u32; 3];
            for dy in 0..f {
                for dx in 0..f {
                    let p = img.get_pixel(x * f + dx, y * f + dy).0;
                    for i in 0..3 {
                        s[i] += u32::from(p[i]);
                    }
                }
            }
            out.put_pixel(x, y, Rgb(s.map(|v| (v / (f * f)) as u8)));
        }
    }
    out
}
