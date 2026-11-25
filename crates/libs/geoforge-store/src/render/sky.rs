//! Star charts: an equirectangular sky with the galaxy's glow, nebulae, stars by magnitude
//! and colour, planets circled, the star's yearly path and the world's constellations.

use image::{Rgb, RgbImage};

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

/// A nebula on a sky chart: a soft patch of `colour`, `strength` (0–1, its surface
/// brightness against the night's) and `radius_deg`; a `shell` is brightest at its edge (a
/// supernova remnant's).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SkyPatch {
    pub ra_deg: f64,
    pub dec_deg: f64,
    pub radius_deg: f64,
    pub colour: [u8; 3],
    pub strength: f32,
    pub shell: bool,
}

/// What a sky chart shows: `glow` a grid `glow_width` wide (0–1) of the galaxy's haze,
/// `zodiacal` the same grid of the zodiacal light (in the glow's units), drawn warm;
/// `patches` the nebulae, under the stars; `points` the stars and bodies; `path` (RA, Dec) a
/// dotted line (the star's yearly path); `lines` constellation lines between (RA, Dec)
/// ends, along great circles.
#[derive(Debug, Clone, Copy, Default)]
pub struct SkyChart<'a> {
    pub glow: &'a [f32],
    pub zodiacal: &'a [f32],
    pub glow_width: usize,
    pub patches: &'a [SkyPatch],
    pub points: &'a [SkyPoint],
    pub path: &'a [(f64, f64)],
    pub lines: &'a [[(f64, f64); 2]],
}

/// An equirectangular sky chart of `chart`, `width` wide (half as tall): RA 360° at the left
/// edge to 0° at the right, as the sky looks from below, Dec +90° at the top. Lines every
/// 30°, the equator brighter.
#[must_use]
pub fn render_sky(width: u32, chart: &SkyChart) -> RgbImage {
    let SkyChart { glow, zodiacal, glow_width, patches, points, path, lines } = *chart;
    let height = width / 2;
    let glow_height = glow.len() / glow_width.max(1);
    let mut acc = vec![[0f32; 3]; (width * height) as usize];
    for (i, a) in acc.iter_mut().enumerate() {
        let (x, y) = (i as u32 % width, i as u32 / width);
        let gx = (f64::from(x) + 0.5) / f64::from(width) * glow_width as f64 - 0.5;
        let gy = (f64::from(y) + 0.5) / f64::from(height) * glow_height as f64 - 0.5;
        let g = bilinear(glow, glow_width, glow_height, gx, gy).powf(0.8) as f32;
        let z = if zodiacal.len() == glow.len() { bilinear(zodiacal, glow_width, glow_height, gx, gy).min(1.5).powf(0.8) as f32 } else { 0.0 };
        *a = [6.0 + 70.0 * g + 72.0 * z, 8.0 + 74.0 * g + 62.0 * z, 16.0 + 92.0 * g + 40.0 * z];
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
    for p in patches {
        patch(&mut acc, width, height, p);
    }
    for &(ra, dec) in path {
        let (x, y) = to_px(ra, dec);
        splat(&mut acc, width, height, x, y, 1.0, [110.0, 90.0, 40.0]);
    }
    let unit = |(ra, dec): (f64, f64)| {
        let (r, d) = (f64::to_radians(ra), f64::to_radians(dec));
        [d.cos() * r.cos(), d.cos() * r.sin(), d.sin()]
    };
    for [a, b] in lines {
        let (u, v) = (unit(*a), unit(*b));
        let angle = (u[0] * v[0] + u[1] * v[1] + u[2] * v[2]).clamp(-1.0, 1.0).acos();
        let steps = (angle.to_degrees() / 360.0 * f64::from(width) * 1.5).ceil().max(1.0) as usize;
        for i in 0..=steps {
            // Along the great circle from u to v.
            let t = i as f64 / steps as f64;
            let (wa, wb) = if angle < 1e-9 { (1.0 - t, t) } else { (((1.0 - t) * angle).sin() / angle.sin(), (t * angle).sin() / angle.sin()) };
            let p = [0, 1, 2].map(|k| wa * u[k] + wb * v[k]);
            let (x, y) = to_px(p[1].atan2(p[0]).to_degrees(), p[2].clamp(-1.0, 1.0).asin().to_degrees());
            splat(&mut acc, width, height, x, y, 0.8, [36.0, 52.0, 80.0]);
        }
    }
    let mut sorted: Vec<&SkyPoint> = points.iter().collect();
    sorted.sort_by(|a, b| b.v.total_cmp(&a.v));
    for p in sorted {
        let (x, y) = to_px(p.ra_deg, p.dec_deg);
        draw_point(&mut acc, width, height, x, y, p, 1.0);
    }
    let mut img = RgbImage::new(width, height);
    for (i, a) in acc.iter().enumerate() {
        img.put_pixel(i as u32 % width, i as u32 / width, Rgb(a.map(|v| v.clamp(0.0, 255.0) as u8)));
    }
    img
}

/// A key strip for [`render_sky`]: each point drawn as on the chart but `scale` times the
/// size, centred in its own `cell` px square on the chart's night, left to right (their RA
/// and Dec are ignored).
#[must_use]
pub fn render_sky_key(cell: u32, scale: f64, points: &[SkyPoint]) -> RgbImage {
    let (width, height) = (cell * points.len().max(1) as u32, cell);
    let mut acc = vec![[6f32, 8.0, 16.0]; (width * height) as usize];
    for (i, p) in points.iter().enumerate() {
        let (x, y) = ((i as f64 + 0.5) * f64::from(cell), 0.5 * f64::from(cell));
        draw_point(&mut acc, width, height, x, y, p, scale);
    }
    let mut img = RgbImage::new(width, height);
    for (i, a) in acc.iter().enumerate() {
        img.put_pixel(i as u32 % width, i as u32 / width, Rgb(a.map(|v| v.clamp(0.0, 255.0) as u8)));
    }
    img
}

/// A star by its magnitude and colour, circled if it is a planet or moon; `scale` times the
/// chart's size.
fn draw_point(acc: &mut [[f32; 3]], width: u32, height: u32, x: f64, y: f64, p: &SkyPoint, scale: f64) {
    let flux = 10f64.powf(-0.4 * (p.v - 6.5)).max(1.0);
    let radius = scale * (0.6 + 0.85 * flux.log10()).min(12.0);
    let peak = (0.35 + 0.3 * flux.log10()).min(1.0);
    let c = star_rgb(p.temperature_k).map(|v| f32::from(v) * peak as f32);
    splat(acc, width, height, x, y, radius, c);
    if p.ring {
        let r = radius + 4.0 * scale;
        for k in 0..64 {
            let a = std::f64::consts::TAU * f64::from(k) / 64.0;
            splat(acc, width, height, x + r * a.cos(), y + r * a.sin(), 0.5, [120.0, 200.0, 140.0]);
        }
    }
}

/// Add nebula `p`: a soft disc, or a shell, at least a pixel and a half across.
fn patch(acc: &mut [[f32; 3]], width: u32, height: u32, p: &SkyPatch) {
    let px_per_deg = f64::from(width) / 360.0;
    let r_px = (p.radius_deg * px_per_deg).max(1.5);
    let (cx, cy) = ((1.0 - p.ra_deg.rem_euclid(360.0) / 360.0) * f64::from(width), (90.0 - p.dec_deg) / 180.0 * f64::from(height));
    // RA squeezes toward the poles: a patch spans 1/cos(Dec) more columns.
    let stretch = 1.0 / p.dec_deg.to_radians().cos().max(0.05);
    let reach_y = (1.6 * r_px).ceil() as i64;
    let reach_x = (1.6 * r_px * stretch).ceil().min(f64::from(width) / 2.0) as i64;
    let c = p.colour.map(|v| f32::from(v) * 0.5 * p.strength);
    for dy in -reach_y..=reach_y {
        let py = cy.floor() as i64 + dy;
        if py < 0 || py >= i64::from(height) {
            continue;
        }
        for dx in -reach_x..=reach_x {
            let px = (cx.floor() as i64 + dx).rem_euclid(i64::from(width));
            let (ex, ey) = ((cx.floor() + dx as f64 + 0.5 - cx) / stretch, py as f64 + 0.5 - cy);
            let r = (ex * ex + ey * ey).sqrt() / r_px;
            let w = if p.shell { (-((r - 0.9) / 0.15).powi(2)).exp() } else { (-2.0 * r * r).exp() } as f32;
            if w > 0.01 {
                let a = &mut acc[(py as u32 * width + px as u32) as usize];
                for k in 0..3 {
                    a[k] += c[k] * w;
                }
            }
        }
    }
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
