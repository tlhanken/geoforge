//! Weather maps: isobars, rain, fronts and storm centres over the land.

use super::blend;
use image::{Rgb, RgbImage};

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
