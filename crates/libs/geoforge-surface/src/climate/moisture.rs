//! Water vapour carried by the winds, raining out where air rises or cools.
//!
//! Each step moves the column's vapour one grid step along the season's wind, rains out what
//! the column's humidity and the air's rising allow, and adds what the surface evaporates:
//! the sea by the bulk formula, land what its soil water gives (see [`soil_evaporation`]).
//! Eddies mix the vapour across the storm tracks.

use crate::atmosphere::transport::{carry, diffuse_rows, semi_lagrangian};
use crate::atmosphere::{monsoon_flow, Atmo, Season};
use crate::circulation::WIND_REF_M_S;
use rayon::prelude::*;

/// Saturation column water (precipitable water at saturation), mm, from temperature.
#[must_use]
pub fn w_sat(t: f64) -> f64 {
    (55.0 * (0.064 * (t - 27.0)).exp()).clamp(0.2, 90.0)
}

/// Share of the ground's seasonal swing the air column follows when its capacity for vapour
/// is reckoned: a winter inversion keeps the air aloft well above the frozen ground, and a hot
/// summer surface heats only the lowest kilometre (Yakutsk swings 57 K at the ground, about
/// 35 K at 850 hPa).
const COLUMN_SWING: f64 = 0.6;
/// Bulk transfer coefficient for evaporation over the sea (Large & Yeager: ~1.2–1.4 × 10⁻³).
const BULK_TRANSFER: f64 = 1.3e-3;
/// Scale height of water vapour, m: the column holds the surface air's humidity over about this
/// depth.
const VAPOUR_SCALE_HEIGHT_M: f64 = 2200.0;
/// Gusts and convection keep evaporating where the mean wind is calm, m/s.
const GUST_M_S: f64 = 4.0;
/// Eddy diffusivity of vapour, m²/s: weak in the tropics, ~10⁶ in the storm tracks (moist
/// energy-balance models take 1 × 10⁶ for the mid-latitudes).
const EDDY_TROPICS_M2_S: f64 = 3.0e4;
const EDDY_STORM_M2_S: f64 = 8.0e5;
/// Drying of the column under sinking air over land, s⁻¹ per unit of ln(1/ascent multiplier):
/// the divergent outflow beneath subsidence carries the low-level vapour away, ~0.6 × 10⁻⁶ s⁻¹
/// under a desert's strong sinking (a multiplier of 0.4), as observed in the subtropics. Over
/// the sea the evaporation it drives replaces it.
const SUBSIDENCE_DRYING_S: f64 = 7.0e-7;
/// How much more a windward slope wrings out at most: lifting rains out the vapour within a
/// scale height or so, so a 4 km rise rains little more than a 1.5 km one.
const OROGRAPHIC_GAIN: f64 = 1.0;

/// Seconds of one moisture step: a wind of [`WIND_REF_M_S`] crosses one grid cell.
#[must_use]
pub fn step_seconds(res_deg: f64, radius_km: f64) -> f64 {
    res_deg.to_radians() * radius_km * 1000.0 / WIND_REF_M_S
}

/// Steps in a year of 365.25 days.
#[must_use]
pub fn steps_per_year(res_deg: f64, radius_km: f64) -> f64 {
    365.25 * 86400.0 / step_seconds(res_deg, radius_km)
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
fn ascent(a: &Atmo, wind: &[(f64, f64)], flat_wind: &[(f64, f64)], season: &Season, monsoon: f64, sst_air: &[f64], radius_km: f64) -> Vec<f64> {
    let g = &a.g;
    let m = monsoon_flow(a, season, monsoon);
    let subtropics = |lat: f64| ((lat.abs() - 8.0) / 6.0).clamp(0.0, 1.0) * ((50.0 - lat.abs()) / 10.0).clamp(0.0, 1.0);
    let mut highs = vec![0.0; g.len()];
    let mut heat_low = vec![0.0; g.len()];
    for i in 0..g.len() {
        let (x, lat) = (i % g.w, g.lat(i / g.w));
        let k = 1.2 * subtropics(lat) * lat.signum();
        highs[i] = k * (flat_wind[i].1 - season.wind(x, lat).1 - m[i].1);
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

/// What one season's moisture run reads.
pub(super) struct SeasonAir<'a> {
    pub atmo: &'a Atmo,
    pub wind: &'a [(f64, f64)],
    /// The season's winds over flat ground (for the subtropical highs' rising and sinking).
    pub flat_wind: &'a [(f64, f64)],
    pub season: &'a Season,
    pub monsoon: f64,
    /// The season's temperature, °C, and the year's mean, per cell.
    pub temp: &'a [f64],
    pub temp_mean: &'a [f64],
    /// Air warmed or chilled by the sea's currents, K.
    pub sst_air: &'a [f64],
    /// Rain-out rate per step.
    pub rainout: f64,
    pub radius_km: f64,
}

/// A season's moisture after a run: rain per step (mm) and the vapour left in each column (mm).
pub(super) struct Run {
    pub rain: Vec<f64>,
    pub vapour: Vec<f64>,
}

/// Run one season for `steps` steps, from `start`'s vapour or a guess, with land evaporating
/// `soil` (mm per step).
pub(super) fn moisture(air: &SeasonAir<'_>, soil: &[f64], start: Option<&[f64]>, steps: usize) -> Run {
    let a = air.atmo;
    let g = &a.g;
    let wind = air.wind;
    let column = |t: f64, m: f64, share: f64| w_sat(m + share * COLUMN_SWING * (t - m));
    let ws: Vec<f64> = air.temp.iter().zip(air.temp_mean).map(|(&t, &m)| column(t, m, 1.0)).collect();
    // Humidity for rain is judged against a slightly damped seasonal temperature (hot
    // continental summers still have convective rain where moist air reaches them).
    let s_rain: Vec<f64> = air.temp.iter().zip(air.temp_mean).map(|(&t, &m)| column(t, m, 0.8)).collect();
    let conv = ascent(a, wind, air.flat_wind, air.season, air.monsoon, air.sst_air, air.radius_km);
    let dt = step_seconds(g.res, air.radius_km);
    let dx = g.res.to_radians() * air.radius_km * 1000.0;
    // Eddies carry vapour across the storm tracks, far into the continents.
    let eddies: Vec<f64> = (0..g.h).map(|y| 4.0 * (EDDY_TROPICS_M2_S + EDDY_STORM_M2_S * air.season.cells.storminess(g.lat(y))) * dt / (dx * dx)).collect();
    // The sea evaporates by the bulk formula, E = C_E |U| (q_sat − q), here on the column.
    let bulk = BULK_TRANSFER * dt / VAPOUR_SCALE_HEIGHT_M;
    let mut w: Vec<f64> = match start {
        Some(w0) => w0.to_vec(),
        None => ws.iter().zip(&a.land).map(|(s, l)| s * (0.75 - 0.4 * l)).collect(),
    };
    let mut p = vec![0.0; g.len()];
    let elev_up: Vec<f64> = semi_lagrangian(g, wind, &a.elev, 0.35);
    // Sinking air over land spreads out near the ground and carries its vapour away.
    let export: Vec<f64> = conv.iter().zip(&a.land).map(|(c, l)| l * SUBSIDENCE_DRYING_S * dt * (1.0 / c.max(0.05)).ln().max(0.0)).collect();
    for _ in 0..steps {
        let up = semi_lagrangian(g, wind, &w, 0.35);
        let res: Vec<(f64, f64)> = (0..g.len())
            .into_par_iter()
            .map(|i| {
                let l = a.land[i];
                let s = ws[i];
                let mut wi = up[i] * (-export[i]).exp();
                // Lifting over terrain rains out; sinking air dries.
                let dz = (a.elev[i] - elev_up[i]) / 1000.0;
                let oro = if dz > 0.0 { 1.0 + OROGRAPHIC_GAIN * (1.0 - (-dz / 1.5).exp()) } else { (1.0 + 2.0 * dz).max(0.35) };
                let rate = (air.rainout * conv[i] * oro).min(0.6);
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
                let speed = wind[i].0.hypot(wind[i].1) * WIND_REF_M_S;
                let e_sea = bulk * speed.max(GUST_M_S) * (s - wi).max(0.0);
                let e_land = soil[i].min((s - wi).max(0.0));
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
    Run { rain: p, vapour: w }
}

/// Fu's form of the Budyko curve's shape parameter (2.6 is the global fit).
const BUDYKO_OMEGA: f64 = 2.6;

/// A year's evaporation from land, mm, from its rain and evaporative demand, mm: the Budyko
/// curve (Fu 1981). Dry land returns nearly all its rain, wet land nearly its demand.
#[must_use]
pub fn budyko_evaporation(rain_mm: f64, pet_mm: f64) -> f64 {
    if rain_mm <= 0.0 || pet_mm <= 0.0 {
        return 0.0;
    }
    let x = pet_mm / rain_mm;
    rain_mm * (1.0 + x - (1.0 + x.powf(BUDYKO_OMEGA)).powf(1.0 / BUDYKO_OMEGA))
}

/// Evaporation from land per moisture step in each season (northern summer, northern winter),
/// per cell: the year's Budyko evaporation from `annual_rain_mm`, spread over the seasons by
/// their evaporative demand `pet_mm` (mm/yr at each season's temperature). Soil water carries
/// the year's rain, winter snow included, into the summer, so a continental interior keeps
/// evaporating through a summer whose own rain could not feed it.
#[must_use]
pub fn soil_evaporation(annual_rain_mm: &[f64], pet_mm: &[[f64; 2]], steps_per_year: f64) -> [Vec<f64>; 2] {
    let season = |k: usize| -> Vec<f64> {
        annual_rain_mm
            .iter()
            .zip(pet_mm)
            .map(|(&pr, pet)| {
                let mean = 0.5 * (pet[0] + pet[1]);
                if mean <= 0.0 {
                    0.0
                } else {
                    budyko_evaporation(pr, mean) * pet[k] / mean / steps_per_year
                }
            })
            .collect()
    };
    [season(0), season(1)]
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Dry land gives back nearly all its rain, wet land nearly its demand, and neither more.
    #[test]
    fn budyko_limits() {
        let dry = budyko_evaporation(50.0, 2000.0);
        assert!(dry > 48.0 && dry <= 50.0, "{dry}");
        let wet = budyko_evaporation(3000.0, 600.0);
        assert!(wet > 570.0 && wet <= 600.0, "{wet}");
        for (p, e) in [(500.0, 500.0), (800.0, 1200.0), (1200.0, 800.0)] {
            let v = budyko_evaporation(p, e);
            assert!(v < p.min(e) && v > 0.5 * p.min(e), "{p} {e} {v}");
        }
        assert_eq!(budyko_evaporation(0.0, 1000.0), 0.0);
        assert_eq!(budyko_evaporation(800.0, 0.0), 0.0);
    }

    /// A season with all the demand gets all the evaporation, and the two seasons add up to
    /// the year's.
    #[test]
    fn soil_water_follows_demand() {
        let spy = 1000.0;
        let [summer, winter] = soil_evaporation(&[600.0, 600.0], &[[1200.0, 0.0], [600.0, 600.0]], spy);
        assert_eq!(winter[0], 0.0);
        for i in 0..2 {
            let year = 0.5 * (summer[i] + winter[i]) * spy;
            assert!((year - budyko_evaporation(600.0, 600.0)).abs() < 1e-9, "{year}");
        }
    }
}
