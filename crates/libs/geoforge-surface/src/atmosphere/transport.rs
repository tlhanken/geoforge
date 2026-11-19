//! Moving quantities with the wind on the atmosphere grid.

use super::grid::{Atmo, AtmoGrid};
use rayon::prelude::*;

/// Upwind sample: each cell takes the value one grid step upstream. Slow winds are sped up
/// to `speed_floor` so stagnant air still mixes.
pub fn semi_lagrangian(g: &AtmoGrid, wind: &[(f64, f64)], f: &[f64], speed_floor: f64) -> Vec<f64> {
    (0..g.len())
        .into_par_iter()
        .map(|i| {
            let (x, y) = (i % g.w, i / g.w);
            let (lat, lon) = (g.lat(y), g.lon(x));
            let (u, v) = wind[i];
            let sp = (u * u + v * v).sqrt();
            let k = if sp < speed_floor { speed_floor / sp.max(1e-9) } else { 1.0 };
            let (u, v) = (u * k, v * k);
            let ulat = (lat - v * g.res).clamp(-89.99, 89.99);
            let ulon = lon - u * g.res / lat.to_radians().cos().max(0.15);
            g.sample(f, ulat, ulon)
        })
        .collect()
}

pub fn diffuse(g: &AtmoGrid, f: &mut [f64], k: f64) {
    let src = f.to_vec();
    for y in 0..g.h {
        for x in 0..g.w {
            let i = y * g.w + x;
            let l = src[y * g.w + (x + g.w - 1) % g.w];
            let r = src[y * g.w + (x + 1) % g.w];
            let u = src[y.saturating_sub(1) * g.w + x];
            let d = src[(y + 1).min(g.h - 1) * g.w + x];
            f[i] = (1.0 - k) * src[i] + k * 0.25 * (l + r + u + d);
        }
    }
}

/// Diffusion with a coefficient per grid row (eddies mix more in the storm tracks), in as
/// many stable sub-steps as the largest coefficient needs.
pub fn diffuse_rows(g: &AtmoGrid, f: &mut [f64], k_row: &[f64]) {
    let k_max = k_row.iter().copied().fold(0.0, f64::max);
    let n = (k_max / 0.2).ceil().max(1.0) as usize;
    for _ in 0..n {
        let src = f.to_vec();
        for y in 0..g.h {
            let k = k_row[y] / n as f64;
            for x in 0..g.w {
                let i = y * g.w + x;
                let l = src[y * g.w + (x + g.w - 1) % g.w];
                let r = src[y * g.w + (x + 1) % g.w];
                let u = src[y.saturating_sub(1) * g.w + x];
                let d = src[(y + 1).min(g.h - 1) * g.w + x];
                f[i] = (1.0 - k) * src[i] + k * 0.25 * (l + r + u + d);
            }
        }
    }
}

/// A quantity carried by the winds: over sea it relaxes toward `sea` (25% per step), over
/// land it decays with distance (`reach_km`). With `sea` = 1 this is the maritime-air
/// fraction; with the sea-surface temperature anomaly it carries current-warmed or
/// upwelling-cooled air inland.
pub fn carry(a: &Atmo, wind: &[(f64, f64)], reach_km: f64, radius_km: f64, sea: &[f64]) -> Vec<f64> {
    let g = &a.g;
    let step_km = g.res.to_radians() * radius_km;
    let decay = (-step_km / reach_km).exp();
    let mut m: Vec<f64> = a.land.iter().zip(sea).map(|(l, s)| (1.0 - l) * s).collect();
    for _ in 0..(3.0 * reach_km / step_km).ceil() as usize + 20 {
        let up = semi_lagrangian(g, wind, &m, 0.35);
        m = up
            .iter()
            .zip(&a.land)
            .zip(sea)
            .map(|((&mu, &l), &s)| {
                let over_land = mu * decay;
                let over_sea = mu + (s - mu) * 0.25;
                l * over_land + (1.0 - l) * over_sea
            })
            .collect();
        diffuse(g, &mut m, 0.15);
    }
    m
}
