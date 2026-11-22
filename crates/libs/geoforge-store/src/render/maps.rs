//! Derived-view maps: classes in their colours, chances on a white-to-red ramp, quantities
//! on a diverging ramp with contours; and their grayscale exports.

use super::{blend, hillshade, image_of, ramp, shade};
use crate::{Result, StoreError};
use geoforge_core::grid::{AnyGrid, Grid, GridSpec};
use geoforge_core::product::Products;
use image::{Rgb, RgbImage};

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

/// The key to [`render_chances`]: chance 0–1.
#[must_use]
pub fn chance_legend() -> super::Legend {
    super::Legend::Ramp { stops: CHANCE_RAMP.to_vec(), units: "chance" }
}

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

/// A preview of a quantity map: blue below its mean through white to red above it (±`span`),
/// with contours every `contour` units, hillshaded with `relief` and coasts drawn.
///
/// # Errors
/// For grids without a raster form.
pub fn render_quantity(values: &[f64], sea: &[bool], contour: f64, span: f64, grid: &AnyGrid, relief: Option<(&Products, f64)>) -> Result<RgbImage> {
    let mean = mean(values);
    let colors: Vec<Option<[u8; 3]>> = values.iter().map(|v| Some(ramp(&DIVERGING, (0.5 + (v - mean) / (2.0 * span)).clamp(0.0, 1.0)))).collect();
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

const DIVERGING: [(f64, [u8; 3]); 3] = [(0.0, [40, 70, 160]), (0.5, [245, 245, 245]), (1.0, [170, 30, 30])];

fn mean(values: &[f64]) -> f64 {
    values.iter().sum::<f64>() / values.len().max(1) as f64
}

/// The key to [`render_quantity`] of `values` with `span`, in `units`.
#[must_use]
pub fn quantity_legend(values: &[f64], span: f64, units: &'static str) -> super::Legend {
    let m = mean(values);
    super::Legend::Ramp { stops: DIVERGING.iter().map(|&(t, c)| (m + (t - 0.5) * 2.0 * span, c)).collect(), units }
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
