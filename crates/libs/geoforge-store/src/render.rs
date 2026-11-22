//! Preview renders (RGB PNG) of products, for looking at a world rather than editing it.
//! Renders use other products for context when present: hillshade from elevation under
//! biomes, rivers and lakes drawn over the land, plate outlines over plates.

use crate::edit::{hsv, key_color};
use crate::{Result, StoreError};
use geoforge_core::body::BodyClass;
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

/// Elevation tint for a world without seas (relative to the mean radius), after the
/// colouring of Mars laser-altimetry maps: blue lows to white highs.
#[must_use]
pub fn dry_elevation_color(e: f64) -> [u8; 3] {
    const STOPS: [(f64, [u8; 3]); 8] = [
        (-9000.0, [30, 20, 90]),
        (-5000.0, [50, 90, 190]),
        (-2000.0, [70, 170, 200]),
        (0.0, [110, 190, 110]),
        (2000.0, [220, 200, 90]),
        (5000.0, [200, 120, 60]),
        (9000.0, [150, 80, 60]),
        (16000.0, [250, 245, 240]),
    ];
    ramp(&STOPS, e)
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
/// Render a product. `sea`: the world has surface water, so ground below 0 m is sea; on a
/// dry world elevation is relative to the mean radius and every cell is land.
/// `class` picks the cloud palette of a giant.
pub fn render(p: ProductId, ctx: &Products, grid: &AnyGrid, radius_km: f64, sea: bool, class: BodyClass) -> Result<Option<RgbImage>> {
    let spec = grid.spec();
    // The rough heightmap is shaded by itself; everything else by the final heightmap.
    let elev = if p == ProductId::RoughElevation { ctx.rough_elevation() } else { ctx.elevation().or(ctx.rough_elevation()) }.map(|f| f.values.as_slice());
    let shade_f = elev.map(|e| hillshade(grid, e, radius_km, 12.0));
    let is_sea = |i: usize| sea && elev.is_some_and(|e| e[i] < 0.0);
    let lit = |i: usize, c: [u8; 3]| match &shade_f {
        Some(s) if !is_sea(i) => shade(c, s[i]),
        _ => c,
    };
    let Some(data) = ctx.get(p) else { return Ok(None) };
    let mut img = match data {
        ProductData::Elevation(f) | ProductData::RoughElevation(f) if sea => image_of(spec, |i| lit(i, elevation_color(f64::from(f.values[i]))))?,
        ProductData::Elevation(f) | ProductData::RoughElevation(f) => image_of(spec, |i| lit(i, dry_elevation_color(f64::from(f.values[i]))))?,
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
            if is_sea(i) { blend(c, [30, 40, 70], 0.55) } else { lit(i, c) }
        })?,
        ProductData::CrustType(f) => image_of(spec, |i| f.values[i].color())?,
        ProductData::Lithology(pt) => image_of(spec, |i| lit(i, pt.row_at(i as u32).map_or([30, 30, 30], |r| r.color)))?,
        ProductData::RockAge(f) => image_of(spec, |i| {
            let a = f64::from(f.values[i]);
            ramp(&[(0.0, [230, 60, 50]), (66.0, [245, 170, 60]), (252.0, [240, 230, 90]), (541.0, [110, 190, 120]), (1000.0, [60, 150, 200]), (2500.0, [70, 70, 170]), (4000.0, [50, 30, 80])], a)
        })?,
        ProductData::SstAnomaly(f) => image_of(spec, |i| match is_sea(i) {
            false => [215, 210, 200],
            true => ramp(&[(-8.0, [30, 60, 160]), (-3.0, [120, 170, 225]), (0.0, [245, 245, 245]), (3.0, [245, 170, 120]), (8.0, [170, 30, 30])], f64::from(f.values[i])),
        })?,
        ProductData::PrecipitationWarmShare(f) => image_of(spec, |i| match is_sea(i) {
            true => [40, 50, 75],
            false => lit(i, ramp(&[(0.15, [180, 90, 40]), (0.35, [235, 200, 120]), (0.5, [240, 240, 230]), (0.65, [140, 200, 150]), (0.85, [40, 110, 160])], f64::from(f.values[i]))),
        })?,
        ProductData::WindJan(f) | ProductData::WindJul(f) | ProductData::OceanCurrents(f) => {
            let ocean = p == ProductId::OceanCurrents;
            let top = if ocean { 1.0 } else { 15.0 };
            let mut img = image_of(spec, |i| {
                let v = f.values[i];
                let sp = f64::from(v[0]).hypot(f64::from(v[1]));
                let base = ramp(&[(0.0, [235, 238, 245]), (top * 0.3, [170, 200, 230]), (top * 0.7, [70, 120, 200]), (top, [40, 40, 140])], sp);
                match is_sea(i) {
                    false => if ocean { [200, 195, 180] } else { blend(base, [180, 170, 150], 0.35) },
                    true => base,
                }
            })?;
            arrows(&mut img, grid, &f.values, top);
            img
        }
        ProductData::Biome(pt) => image_of(spec, |i| lit(i, pt.row_at(i as u32).map_or([30, 30, 30], |r| r.color)))?,
        ProductData::Plates(pt) => image_of(spec, |i| key_color(pt.owner.values[i]))?,
        ProductData::Provinces(pt) => image_of(spec, |i| match pt.row_at(i as u32) {
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
        ProductData::CloudBands(f) => {
            let stops: &[(f64, [u8; 3])] = if class == BodyClass::IceGiant {
                &[(0.0, [60, 110, 175]), (0.5, [130, 190, 225]), (1.0, [215, 240, 250])]
            } else {
                &[(0.0, [120, 70, 45]), (0.35, [175, 120, 80]), (0.6, [220, 195, 160]), (1.0, [250, 245, 232])]
            };
            image_of(spec, |i| ramp(stops, f64::from(f.values[i])))?
        }
        ProductData::ZonalWind(f) => {
            let top = f.values.iter().fold(1.0_f32, |m, v| m.max(v.abs()));
            image_of(spec, |i| ramp(&[(-1.0, [40, 70, 170]), (0.0, [240, 240, 240]), (1.0, [180, 40, 40])], f64::from(f.values[i] / top)))?
        }
        ProductData::Deposits(t) => {
            // Markers over the land, coloured by what they yield.
            let mut img = match elev {
                Some(e) => image_of(spec, |i| {
                    let v = f64::from(e[i]);
                    let c = if is_sea(i) { [200, 212, 225] } else { blend(dry_elevation_color(v.max(0.0)), [235, 235, 235], 0.5) };
                    lit(i, c)
                })?,
                None => image_of(spec, |_| [235, 235, 235])?,
            };
            let (w, h) = (img.width() as i64, img.height() as i64);
            let r = (w / 400).max(2);
            for d in t.rows.values() {
                let x0 = ((d.lon_deg + 180.0) / 360.0 * w as f64) as i64;
                let y0 = ((90.0 - d.lat_deg) / 180.0 * h as f64) as i64;
                for dy in -r - 1..=r + 1 {
                    for dx in -r - 1..=r + 1 {
                        let d2 = dx * dx + dy * dy;
                        let y = y0 + dy;
                        if d2 > (r + 1) * (r + 1) || !(0..h).contains(&y) {
                            continue;
                        }
                        let col = if d2 > r * r { [30, 30, 30] } else { deposit_color(&d.kind) };
                        img.put_pixel((x0 + dx).rem_euclid(w) as u32, y as u32, Rgb(col));
                    }
                }
            }
            img
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

/// Sea in class maps.
pub const CLASS_SEA: [u8; 3] = [225, 232, 240];

/// An image of per-cell class colours, sea (`None`) in [`CLASS_SEA`]. With `relief` (the
/// products and the body's radius) land is hillshaded by the heightmap, for previews; without
/// it every pixel is its class's exact colour, for exports.
///
/// # Errors
/// For grids without a raster form.
pub fn render_classes(colors: &[Option<[u8; 3]>], grid: &AnyGrid, relief: Option<(&Products, f64)>) -> Result<RgbImage> {
    let shading = relief.and_then(|(p, r)| p.elevation().or(p.rough_elevation()).map(|e| hillshade(grid, &e.values, r, 12.0)));
    image_of(grid.spec(), |i| match (colors[i], &shading) {
        (None, _) => CLASS_SEA,
        (Some(c), Some(s)) => shade(c, s[i]),
        (Some(c), None) => c,
    })
}

/// Preview ramp for chance maps: white for 0 through gold and orange to dark red for 1.
const CHANCE_RAMP: [(f64, [u8; 3]); 5] = [(0.0, [250, 250, 245]), (0.1, [250, 225, 140]), (0.35, [240, 150, 60]), (0.7, [200, 50, 40]), (1.0, [110, 0, 30])];

/// A preview of a chance map (values 0–1), hillshaded with `relief`; values at sea (cells
/// where `sea` is set) are tinted blue, and cells without a value drawn as plain sea or land.
///
/// # Errors
/// For grids without a raster form.
pub fn render_chances(values: &[Option<f64>], sea: &[bool], grid: &AnyGrid, relief: Option<(&Products, f64)>) -> Result<RgbImage> {
    let colors: Vec<Option<[u8; 3]>> = values.iter().map(|v| v.map(|x| ramp(&CHANCE_RAMP, x.clamp(0.0, 1.0)))).collect();
    let mut img = render_classes(&colors, grid, relief)?;
    for (i, p) in img.pixels_mut().enumerate() {
        match (sea[i], values[i]) {
            (true, Some(_)) => *p = Rgb(blend(p.0, [150, 180, 215], 0.45)),
            // A sea map's land.
            (false, None) => *p = Rgb([205, 200, 190]),
            _ => {}
        }
    }
    Ok(img)
}

/// An exact 8-bit grayscale image of a chance map: 0–255 for 0–1, 0 where `None`.
///
/// # Errors
/// For grids without a raster form.
pub fn gray_chances(values: &[Option<f64>], grid: &AnyGrid) -> Result<image::GrayImage> {
    let GridSpec::Equirect { width, height } = grid.spec() else {
        return Err(StoreError::Invalid("exports need an equirectangular grid".into()));
    };
    Ok(image::GrayImage::from_fn(width, height, |x, y| {
        let v = values[(y * width + x) as usize].unwrap_or(0.0);
        image::Luma([(v.clamp(0.0, 1.0) * 255.0).round() as u8])
    }))
}

/// Marker colour of a deposit kind.
fn deposit_color(kind: &str) -> [u8; 3] {
    match kind {
        "iron" => [170, 60, 30],
        "copper" => [20, 150, 140],
        "tin" => [110, 110, 120],
        "lead_zinc_silver" => [120, 140, 200],
        "gold" => [230, 180, 0],
        "coal" => [20, 20, 20],
        "salt" => [250, 250, 250],
        "peat" => [100, 70, 40],
        "gemstones" => [200, 40, 200],
        _ => [255, 0, 0],
    }
}

/// A preview of a quantity map: blue below its mean through white to red above it (±`span`),
/// with contours every `contour` units, hillshaded with `relief` and coasts drawn.
///
/// # Errors
/// For grids without a raster form.
pub fn render_quantity(values: &[f64], sea: &[bool], contour: f64, span: f64, grid: &AnyGrid, relief: Option<(&Products, f64)>) -> Result<RgbImage> {
    let mean = values.iter().sum::<f64>() / values.len().max(1) as f64;
    let diverging = [(0.0, [40, 70, 160]), (0.5, [245, 245, 245]), (1.0, [170, 30, 30])];
    let colors: Vec<Option<[u8; 3]>> = values.iter().map(|v| Some(ramp(&diverging, (0.5 + (v - mean) / (2.0 * span)).clamp(0.0, 1.0)))).collect();
    let mut img = render_classes(&colors, grid, relief)?;
    let GridSpec::Equirect { width, height } = grid.spec() else { return Ok(img) };
    let band = |i: usize| (values[i] / contour).floor();
    for y in 0..height as usize {
        for x in 0..width as usize {
            let i = y * width as usize + x;
            let right = y * width as usize + (x + 1) % width as usize;
            let below = (y + 1).min(height as usize - 1) * width as usize + x;
            let coast = sea[i] != sea[right] || sea[i] != sea[below];
            if band(i) != band(right) || band(i) != band(below) {
                img.put_pixel(x as u32, y as u32, Rgb([60, 60, 60]));
            } else if coast {
                img.put_pixel(x as u32, y as u32, Rgb([110, 110, 110]));
            }
        }
    }
    Ok(img)
}

/// A 16-bit grayscale image.
pub type Gray16 = image::ImageBuffer<image::Luma<u16>, Vec<u16>>;

/// An exact 16-bit grayscale image of a quantity map and the range it spans: 0 is `min`,
/// 65535 is `max`.
///
/// # Errors
/// For grids without a raster form.
pub fn gray_quantity(values: &[f64], grid: &AnyGrid) -> Result<(Gray16, f64, f64)> {
    let GridSpec::Equirect { width, height } = grid.spec() else {
        return Err(StoreError::Invalid("exports need an equirectangular grid".into()));
    };
    let (min, max) = values.iter().fold((f64::MAX, f64::MIN), |(lo, hi), &v| (lo.min(v), hi.max(v)));
    let span = (max - min).max(1e-9);
    let img = image::ImageBuffer::from_fn(width, height, |x, y| image::Luma([((values[(y * width + x) as usize] - min) / span * 65535.0).round() as u16]));
    Ok((img, min, max))
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

/// Something drawn on a sky chart.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SkyPoint {
    pub ra_deg: f64,
    pub dec_deg: f64,
    /// Visual magnitude.
    pub v: f64,
    pub temperature_k: f64,
    /// Circled: a planet or moon rather than a star.
    pub ring: bool,
}

/// An equirectangular sky chart `width` wide (half as tall): RA 360° at the left edge to 0°
/// at the right, as the sky looks from below, Dec +90° at the top. `glow` is a grid
/// `glow_width` wide (0–1), drawn as the galaxy's haze; `path` (RA, Dec) a dotted line (the
/// star's yearly path). Lines every 30°, the equator brighter.
#[must_use]
pub fn render_sky(width: u32, glow: &[f32], glow_width: usize, points: &[SkyPoint], path: &[(f64, f64)]) -> RgbImage {
    let height = width / 2;
    let glow_height = glow.len() / glow_width.max(1);
    let mut acc = vec![[0f32; 3]; (width * height) as usize];
    for (i, a) in acc.iter_mut().enumerate() {
        let (x, y) = (i as u32 % width, i as u32 / width);
        let gx = (f64::from(x) + 0.5) / f64::from(width) * glow_width as f64 - 0.5;
        let gy = (f64::from(y) + 0.5) / f64::from(height) * glow_height as f64 - 0.5;
        let g = bilinear(glow, glow_width, glow_height, gx, gy).powf(0.8) as f32;
        *a = [6.0 + 70.0 * g, 8.0 + 74.0 * g, 16.0 + 92.0 * g];
        let (ra, dec) = (360.0 * (1.0 - (f64::from(x) + 0.5) / f64::from(width)), 90.0 - 180.0 * (f64::from(y) + 0.5) / f64::from(height));
        let step = 360.0 / f64::from(width);
        let line = |v: f64| (v / 30.0 - (v / 30.0).round()).abs() * 30.0 < step * 0.5;
        if dec.abs() < step * 0.5 {
            *a = [a[0] + 40.0, a[1] + 48.0, a[2] + 70.0];
        } else if line(ra) || line(dec) {
            *a = [a[0] + 14.0, a[1] + 18.0, a[2] + 30.0];
        }
    }
    let to_px = |ra: f64, dec: f64| ((1.0 - ra.rem_euclid(360.0) / 360.0) * f64::from(width), (90.0 - dec) / 180.0 * f64::from(height));
    for &(ra, dec) in path {
        let (x, y) = to_px(ra, dec);
        splat(&mut acc, width, height, x, y, 1.0, [110.0, 90.0, 40.0]);
    }
    let mut sorted: Vec<&SkyPoint> = points.iter().collect();
    sorted.sort_by(|a, b| b.v.total_cmp(&a.v));
    for p in sorted {
        let (x, y) = to_px(p.ra_deg, p.dec_deg);
        let flux = 10f64.powf(-0.4 * (p.v - 6.5)).max(1.0);
        let radius = (0.6 + 0.85 * flux.log10()).min(12.0);
        let peak = (0.35 + 0.3 * flux.log10()).min(1.0);
        let c = star_rgb(p.temperature_k).map(|v| f32::from(v) * peak as f32);
        splat(&mut acc, width, height, x, y, radius, c);
        if p.ring {
            let r = radius + 4.0;
            for k in 0..64 {
                let a = std::f64::consts::TAU * f64::from(k) / 64.0;
                splat(&mut acc, width, height, x + r * a.cos(), y + r * a.sin(), 0.5, [120.0, 200.0, 140.0]);
            }
        }
    }
    let mut img = RgbImage::new(width, height);
    for (i, a) in acc.iter().enumerate() {
        img.put_pixel(i as u32 % width, i as u32 / width, Rgb(a.map(|v| v.clamp(0.0, 255.0) as u8)));
    }
    img
}

/// Add a soft spot of colour `c` (peak) and size `radius` px at (x, y), wrapping in RA.
fn splat(acc: &mut [[f32; 3]], width: u32, height: u32, x: f64, y: f64, radius: f64, c: [f32; 3]) {
    let sigma = (radius * 0.5).max(0.35);
    let reach = (radius + 1.5).ceil() as i64;
    for dy in -reach..=reach {
        let py = y.floor() as i64 + dy;
        if py < 0 || py >= i64::from(height) {
            continue;
        }
        for dx in -reach..=reach {
            let px = (x.floor() as i64 + dx).rem_euclid(i64::from(width));
            let (ex, ey) = (x.floor() + dx as f64 + 0.5 - x, py as f64 + 0.5 - y);
            let w = (-(ex * ex + ey * ey) / (2.0 * sigma * sigma)).exp() as f32;
            if w > 0.01 {
                let a = &mut acc[(py as u32 * width + px as u32) as usize];
                for k in 0..3 {
                    a[k] += c[k] * w;
                }
            }
        }
    }
}

fn bilinear(g: &[f32], w: usize, h: usize, x: f64, y: f64) -> f64 {
    if g.is_empty() {
        return 0.0;
    }
    let y = y.clamp(0.0, (h - 1) as f64);
    let (x0, y0) = (x.floor(), y.floor());
    let (fx, fy) = (x - x0, y - y0);
    let at = |i: i64, j: i64| f64::from(g[(j.clamp(0, h as i64 - 1) as usize) * w + i.rem_euclid(w as i64) as usize]);
    let (i, j) = (x0 as i64, y0 as i64);
    let top = at(i, j) * (1.0 - fx) + at(i + 1, j) * fx;
    let bottom = at(i, j + 1) * (1.0 - fx) + at(i + 1, j + 1) * fx;
    top * (1.0 - fy) + bottom * fy
}

/// The colour of a star of temperature `t` K, as the eye sees it: blackbody hues (Helland's
/// fit to Charity's table) halfway to white.
#[must_use]
pub fn star_rgb(t: f64) -> [u8; 3] {
    let t = t.clamp(1000.0, 40_000.0) / 100.0;
    let r = if t <= 66.0 { 255.0 } else { 329.7 * (t - 60.0).powf(-0.1332) };
    let g = if t <= 66.0 { 99.47 * t.ln() - 161.12 } else { 288.12 * (t - 60.0).powf(-0.0755) };
    let b = if t >= 66.0 {
        255.0
    } else if t <= 19.0 {
        0.0
    } else {
        138.52 * (t - 10.0).ln() - 305.04
    };
    [r, g, b].map(|v: f64| (0.5 * v.clamp(0.0, 255.0) + 127.5) as u8)
}

/// A weather map: `base` colours (row by row, `width` wide) under isobars of `pressure`
/// (hPa) every 4 hPa, rain shaded by `rain` (0–1), fronts as (lat, lon) lines (warm red,
/// cold blue) and storm centres (lows red, tropical cyclones magenta).
#[must_use]
pub fn render_weather(width: u32, base: &[[u8; 3]], pressure: &[f32], rain: &[f32], fronts: &[(Vec<(f64, f64)>, bool)], centres: &[(f64, f64, bool)]) -> RgbImage {
    let height = base.len() as u32 / width;
    let (w, h) = (width as usize, height as usize);
    let band = |v: f32| (v / 4.0).floor() as i32;
    let mut img = RgbImage::new(width, height);
    for y in 0..h {
        for x in 0..w {
            let i = y * w + x;
            let mut c = base[i];
            let r = rain[i].clamp(0.0, 1.0);
            if r > 0.15 {
                c = blend(c, [40, 170, 90], f64::from(0.25 + 0.5 * (r - 0.15) / 0.85));
            }
            let right = y * w + (x + 1) % w;
            let below = (y + 1).min(h - 1) * w + x;
            if band(pressure[i]) != band(pressure[right]) || band(pressure[i]) != band(pressure[below]) {
                c = blend(c, [255, 255, 255], 0.7);
            }
            img.put_pixel(x as u32, y as u32, Rgb(c));
        }
    }
    let to_px = |lat: f64, lon: f64| ((lon + 180.0) / 360.0 * f64::from(width), (90.0 - lat) / 180.0 * f64::from(height));
    let dot = |img: &mut RgbImage, x: f64, y: f64, rad: f64, col: [u8; 3]| {
        let r = rad.ceil() as i64;
        for dy in -r..=r {
            for dx in -r..=r {
                if ((dx * dx + dy * dy) as f64) <= rad * rad {
                    let px = (x as i64 + dx).rem_euclid(i64::from(width));
                    let py = y as i64 + dy;
                    if (0..i64::from(height)).contains(&py) {
                        img.put_pixel(px as u32, py as u32, Rgb(col));
                    }
                }
            }
        }
    };
    for (line, warm) in fronts {
        let col = if *warm { [220, 40, 40] } else { [40, 80, 230] };
        for pair in line.windows(2) {
            let ((x0, y0), (x1, y1)) = (to_px(pair[0].0, pair[0].1), to_px(pair[1].0, pair[1].1));
            let dx = (x1 - x0 + f64::from(width) / 2.0).rem_euclid(f64::from(width)) - f64::from(width) / 2.0;
            let n = (dx.hypot(y1 - y0).ceil() as usize).max(1);
            for k in 0..=n {
                let t = k as f64 / n as f64;
                dot(&mut img, x0 + t * dx, y0 + t * (y1 - y0), 1.5, col);
            }
        }
    }
    for &(lat, lon, tropical) in centres {
        let (x, y) = to_px(lat, lon);
        dot(&mut img, x, y, 5.0, [255, 255, 255]);
        dot(&mut img, x, y, 3.5, if tropical { [210, 40, 210] } else { [220, 30, 30] });
    }
    img
}
