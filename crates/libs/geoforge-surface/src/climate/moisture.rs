//! Water vapour carried by the winds, raining out where air rises or cools.

use crate::atmosphere::transport::{carry, diffuse_rows, semi_lagrangian};
use crate::atmosphere::{monsoon_flow, Atmo, Season};
use rayon::prelude::*;

/// Saturation column water (precipitable water at saturation), mm, from temperature.
#[must_use]
pub fn w_sat(t: f64) -> f64 {
    (55.0 * (0.064 * (t - 27.0)).exp()).clamp(0.2, 90.0)
}

/// Rain-out multiplier from large-scale vertical motion, per atmosphere cell:
/// - the season's cells and ITCZ;
/// - the subtropical highs: in the subtropics air turned poleward rises and air turned
///   equatorward sinks (the Sverdrup balance). The sinking on the eastern side of each ocean's
///   high travels inland with the wind and dries the west coasts; the rising on its western
///   side wets the east coasts;
/// - the northerlies on the poleward side of a summer heat low sink (the monsoon inflow on its
///   equatorward side is shallow and adds nothing);
/// - the stable air over cold water or fresh from it.
fn ascent(a: &Atmo, wind: &[(f64, f64)], season: &Season, monsoon: f64, sst_air: &[f64], radius_km: f64) -> Vec<f64> {
    let g = &a.g;
    let m = monsoon_flow(a, season, monsoon);
    let subtropics = |lat: f64| ((lat.abs() - 8.0) / 6.0).clamp(0.0, 1.0) * ((50.0 - lat.abs()) / 10.0).clamp(0.0, 1.0);
    let mut highs = vec![0.0; g.len()];
    let mut heat_low = vec![0.0; g.len()];
    for i in 0..g.len() {
        let (x, lat) = (i % g.w, g.lat(i / g.w));
        let k = 1.2 * subtropics(lat) * lat.signum();
        highs[i] = k * (wind[i].1 - season.wind(x, lat).1 - m[i].1);
        heat_low[i] = (k * m[i].1).min(0.0);
    }
    let highs = carry(a, wind, 800.0, radius_km, &highs);
    (0..g.len())
        .map(|i| {
            let stationary = (highs[i] + heat_low[i]).clamp(-1.5, 0.7).exp();
            let stable = (sst_air[i].min(0.0) / 2.0).exp();
            season.convergence(i % g.w, g.lat(i / g.w), a.land[i]) * stationary * stable
        })
        .collect()
}

/// One season of the moisture model; returns raw precipitation (mm per step).
#[allow(clippy::too_many_arguments)]
pub(super) fn moisture(a: &Atmo, wind: &[(f64, f64)], season: &Season, monsoon: f64, temp: &[f64], temp_mean: &[f64], sst_air: &[f64], rainout: f64, radius_km: f64) -> Vec<f64> {
    let g = &a.g;
    let ws: Vec<f64> = temp.iter().map(|&t| w_sat(t)).collect();
    // Humidity for rain is judged against a slightly damped seasonal temperature (hot
    // continental summers still have convective rain where moist air reaches them).
    let s_rain: Vec<f64> = temp.iter().zip(temp_mean).map(|(&t, &m)| w_sat(m + 0.8 * (t - m))).collect();
    let conv = ascent(a, wind, season, monsoon, sst_air, radius_km);
    // Eddies carry vapour across the storm tracks, far into the continents.
    let eddies: Vec<f64> = (0..g.h).map(|y| 0.12 + 0.35 * season.cells.storminess(g.lat(y))).collect();
    let mut w: Vec<f64> = ws.iter().zip(&a.land).map(|(s, l)| s * (0.75 - 0.4 * l)).collect();
    let mut p = vec![0.0; g.len()];
    let elev_up_cache: Vec<f64> = semi_lagrangian(g, wind, &a.elev, 0.35);
    for _ in 0..260 {
        let up = semi_lagrangian(g, wind, &w, 0.35);
        let res: Vec<(f64, f64)> = (0..g.len())
            .into_par_iter()
            .map(|i| {
                let l = a.land[i];
                let s = ws[i];
                let mut wi = up[i];
                // Lifting over terrain rains out; sinking air dries.
                let dz = (a.elev[i] - elev_up_cache[i]) / 1000.0;
                let oro = if dz > 0.0 { 1.0 + 2.0 * dz } else { (1.0 + 2.0 * dz).max(0.35) };
                let rate = (rainout * conv[i] * oro).min(0.6);
                // Rain needs humid air: rising air rains from ~30% of saturation, sinking air
                // only above ~70%.
                let humid = (0.55 - 0.22 * (conv[i] - 1.0)).clamp(0.28, 0.72) - 0.05 * l;
                let mut pi = rate * (wi - humid * s_rain[i]).max(0.0);
                wi -= pi;
                // Air arriving over colder (or higher) ground than it left condenses its excess
                // over a few steps, not at once: the ground warms the air as it cools it.
                if wi > s {
                    let d = 0.2 * (wi - s);
                    pi += d;
                    wi -= d;
                }
                // Evaporation: the sea relaxes air toward 85% humidity; land recycles rain, more of it
                // the warmer the season (evaporation is limited by the sun's energy).
                let t = temp[i];
                let recycle = 0.65 * ((t + 5.0) / 25.0).clamp(0.0, 1.0) + 0.15 * ((t - 18.0) / 8.0).clamp(0.0, 1.0);
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
        diffuse_rows(g, &mut w, &eddies);
    }
    p
}
