//! The atmosphere's own coarse latitude–longitude grid, and the world's surface on it.

use geoforge_core::grid::Grid;
use geoforge_core::UnitVec3;

/// A regular latitude–longitude grid for the atmosphere.
#[derive(Debug, Clone, Copy)]
pub struct AtmoGrid {
    pub w: usize,
    pub h: usize,
    pub res: f64,
}

impl AtmoGrid {
    #[must_use]
    pub fn new(res_deg: f64) -> Self {
        let res = 180.0 / (180.0 / res_deg.clamp(0.25, 5.0)).round();
        Self { w: (360.0 / res).round() as usize, h: (180.0 / res).round() as usize, res }
    }
    #[must_use]
    pub fn len(&self) -> usize {
        self.w * self.h
    }
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
    #[must_use]
    pub fn lat(&self, y: usize) -> f64 {
        90.0 - (y as f64 + 0.5) * self.res
    }
    #[must_use]
    pub fn lon(&self, x: usize) -> f64 {
        -180.0 + (x as f64 + 0.5) * self.res
    }
    #[must_use]
    pub fn index(&self, lat: f64, lon: f64) -> usize {
        let y = (((90.0 - lat) / self.res).floor().max(0.0) as usize).min(self.h - 1);
        let x = (((lon + 180.0) / self.res).floor().rem_euclid(self.w as f64) as usize).min(self.w - 1);
        y * self.w + x
    }
    /// Bilinear sample with longitude wrap and latitude clamp.
    #[must_use]
    pub fn sample(&self, f: &[f64], lat: f64, lon: f64) -> f64 {
        let fy = ((90.0 - lat) / self.res - 0.5).clamp(0.0, (self.h - 1) as f64);
        let fx = (lon + 180.0) / self.res - 0.5;
        let (y0, ty) = (fy.floor() as usize, fy.fract());
        let y1 = (y0 + 1).min(self.h - 1);
        let x0f = fx.floor();
        let tx = fx - x0f;
        let x0 = (x0f as i64).rem_euclid(self.w as i64) as usize;
        let x1 = (x0 + 1) % self.w;
        let a = f[y0 * self.w + x0] * (1.0 - tx) + f[y0 * self.w + x1] * tx;
        let b = f[y1 * self.w + x0] * (1.0 - tx) + f[y1 * self.w + x1] * tx;
        a * (1.0 - ty) + b * ty
    }
    /// A smooth blur of about the spread of a box blur of `radius` cells: three box blurs,
    /// close to a Gaussian, so the field's slopes have no straight edges.
    #[must_use]
    pub fn smooth(&self, f: &[f64], radius: usize) -> Vec<f64> {
        // Three boxes of radius r spread as far as one of radius R when r(r + 1) = R(R + 1)/3.
        let r = ((((radius * (radius + 1)) as f64 / 3.0 + 0.25).sqrt() - 0.5).round() as usize).max(1);
        self.blur(&self.blur(&self.blur(f, r), r), r)
    }
    /// Separable box blur, radius in cells (longitude radius widened toward the poles).
    #[must_use]
    pub fn blur(&self, f: &[f64], radius: usize) -> Vec<f64> {
        let mut tmp = vec![0.0; f.len()];
        for y in 0..self.h {
            let c = self.lat(y).to_radians().cos().max(0.05);
            let r = ((radius as f64 / c).round() as usize).min(self.w / 2);
            for x in 0..self.w {
                let mut s = 0.0;
                for dx in 0..=2 * r {
                    let xx = (x + self.w + dx - r) % self.w;
                    s += f[y * self.w + xx];
                }
                tmp[y * self.w + x] = s / (2 * r + 1) as f64;
            }
        }
        let mut out = vec![0.0; f.len()];
        for y in 0..self.h {
            for x in 0..self.w {
                let (lo, hi) = (y.saturating_sub(radius), (y + radius).min(self.h - 1));
                let s: f64 = (lo..=hi).map(|yy| tmp[yy * self.w + x]).sum();
                out[y * self.w + x] = s / (hi - lo + 1) as f64;
            }
        }
        out
    }
}

/// Sea: below datum on a world with surface water. A dry world has none.
#[must_use]
pub fn is_sea(wet: bool, elevation_m: f64) -> bool {
    wet && elevation_m < 0.0
}

/// Height of what the air sits on: sea level over the sea, the ground elsewhere, m.
#[must_use]
pub fn ground_m(wet: bool, elevation_m: f64) -> f64 {
    if is_sea(wet, elevation_m) {
        0.0
    } else {
        elevation_m
    }
}

/// The world's surface on the atmosphere grid.
pub struct Atmo {
    pub g: AtmoGrid,
    /// Surface water exists, so cells below datum are sea.
    pub wet: bool,
    /// Land fraction per atmosphere cell.
    pub land: Vec<f64>,
    /// Mean height of what the air sits on (see [`ground_m`]), m.
    pub elev: Vec<f64>,
}

impl Atmo {
    /// Aggregate a world-grid elevation onto an atmosphere grid.
    #[must_use]
    pub fn from_world(grid: &dyn Grid, elev: &[f32], res_deg: f64, wet: bool) -> Self {
        let g = AtmoGrid::new(res_deg);
        let mut sum_e = vec![0.0; g.len()];
        let mut sum_l = vec![0.0; g.len()];
        let mut sum_a = vec![0.0; g.len()];
        for c in grid.cells() {
            let ll = grid.center(c).to_lat_lon();
            let i = g.index(ll.lat_deg, ll.lon_deg);
            let a = grid.area_sr(c);
            let e = f64::from(elev[c as usize]);
            sum_e[i] += ground_m(wet, e) * a;
            sum_l[i] += if is_sea(wet, e) { 0.0 } else { a };
            sum_a[i] += a;
        }
        let mut land = vec![0.0; g.len()];
        let mut elev_c = vec![0.0; g.len()];
        for i in 0..g.len() {
            if sum_a[i] > 0.0 {
                land[i] = sum_l[i] / sum_a[i];
                elev_c[i] = sum_e[i] / sum_a[i];
            } else {
                let (x, y) = (i % g.w, i / g.w);
                let c = grid.locate(UnitVec3::from_lat_lon(g.lat(y), g.lon(x)));
                let e = f64::from(elev[c as usize]);
                land[i] = if is_sea(wet, e) { 0.0 } else { 1.0 };
                elev_c[i] = ground_m(wet, e);
            }
        }
        Self { g, wet, land, elev: elev_c }
    }

    /// Area-weighted mean of a world-grid vector field on the atmosphere grid.
    #[must_use]
    pub fn vectors_from_world(&self, grid: &dyn Grid, v: &[[f32; 2]], scale: f64) -> Vec<(f64, f64)> {
        let g = &self.g;
        let mut su = vec![0.0; g.len()];
        let mut sv = vec![0.0; g.len()];
        let mut sa = vec![0.0; g.len()];
        for c in grid.cells() {
            let ll = grid.center(c).to_lat_lon();
            let i = g.index(ll.lat_deg, ll.lon_deg);
            let a = grid.area_sr(c);
            su[i] += f64::from(v[c as usize][0]) * a;
            sv[i] += f64::from(v[c as usize][1]) * a;
            sa[i] += a;
        }
        (0..g.len())
            .map(|i| {
                if sa[i] > 0.0 {
                    (su[i] / sa[i] / scale, sv[i] / sa[i] / scale)
                } else {
                    let (x, y) = (i % g.w, i / g.w);
                    let c = grid.locate(UnitVec3::from_lat_lon(g.lat(y), g.lon(x))) as usize;
                    (f64::from(v[c][0]) / scale, f64::from(v[c][1]) / scale)
                }
            })
            .collect()
    }

    /// Area-weighted mean of a world-grid scalar over sea cells only (0 where no sea).
    #[must_use]
    pub fn sea_mean_from_world(&self, grid: &dyn Grid, f: &[f32], elev: &[f32]) -> Vec<f64> {
        let g = &self.g;
        let mut s = vec![0.0; g.len()];
        let mut a = vec![0.0; g.len()];
        for c in grid.cells() {
            if !is_sea(self.wet, f64::from(elev[c as usize])) {
                continue;
            }
            let ll = grid.center(c).to_lat_lon();
            let i = g.index(ll.lat_deg, ll.lon_deg);
            let w = grid.area_sr(c);
            s[i] += f64::from(f[c as usize]) * w;
            a[i] += w;
        }
        s.iter().zip(&a).map(|(s, a)| if *a > 0.0 { s / a } else { 0.0 }).collect()
    }
}
