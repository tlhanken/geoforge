//! Water vapour carried by the winds, raining out where air rises or cools.

use crate::atmosphere::transport::{diffuse, semi_lagrangian};
use crate::atmosphere::{Atmo, Cells};
use rayon::prelude::*;

/// Saturation column water (precipitable water at saturation), mm, from temperature.
#[must_use]
pub fn w_sat(t: f64) -> f64 {
    (55.0 * (0.064 * (t - 27.0)).exp()).clamp(0.2, 90.0)
}

/// One season of the moisture model; returns raw precipitation (mm per step).
pub(super) fn moisture(a: &Atmo, wind: &[(f64, f64)], circ: &Cells, temp: &[f64], temp_mean: &[f64], rainout: f64) -> Vec<f64> {
    let g = &a.g;
    let ws: Vec<f64> = temp.iter().map(|&t| w_sat(t)).collect();
    // Hot continental summers do not stop convective rain: judge humidity against a damped
    // seasonal temperature.
    let s_rain: Vec<f64> = temp.iter().zip(temp_mean).map(|(&t, &m)| w_sat(m + 0.5 * (t - m))).collect();
    let conv: Vec<f64> = (0..g.h).map(|y| circ.convergence(g.lat(y))).collect();
    let mut w: Vec<f64> = ws.iter().zip(&a.land).map(|(s, l)| s * (0.75 - 0.4 * l)).collect();
    let mut p = vec![0.0; g.len()];
    let elev_up_cache: Vec<f64> = semi_lagrangian(g, wind, &a.elev, 0.35);
    for _ in 0..260 {
        let up = semi_lagrangian(g, wind, &w, 0.35);
        let res: Vec<(f64, f64)> = (0..g.len())
            .into_par_iter()
            .map(|i| {
                let y = i / g.w;
                let l = a.land[i];
                let s = ws[i];
                let mut wi = up[i];
                // Lifting over terrain rains out; sinking air dries.
                let dz = (a.elev[i] - elev_up_cache[i]) / 1000.0;
                let oro = if dz > 0.0 { 1.0 + 3.0 * dz } else { (1.0 + 2.0 * dz).max(0.35) };
                let rate = (rainout * conv[y] * oro).min(0.85);
                // Rain needs humid air: only vapour above ~25–35% of saturation can rain out.
                let mut pi = rate * (wi - (0.35 - 0.1 * l) * s_rain[i]).max(0.0);
                wi -= pi;
                if wi > s {
                    pi += wi - s;
                    wi = s;
                }
                // Evaporation: the sea relaxes air toward 85% humidity; land recycles rain.
                let t = temp[i];
                let recycle = 0.65 * ((t + 5.0) / 15.0).clamp(0.0, 1.0) + 0.15 * ((t - 18.0) / 8.0).clamp(0.0, 1.0);
                let e_sea = 0.28 * (0.85 * s - wi).max(0.0);
                let e_land = recycle * p[i];
                wi += (1.0 - l) * e_sea + l * e_land;
                (wi.max(0.0), pi)
            })
            .collect();
        for (i, (wi, pi)) in res.into_iter().enumerate() {
            w[i] = wi;
            p[i] = pi;
        }
        diffuse(g, &mut w, 0.12);
    }
    p
}
