//! Climate layer: mean annual temperature, annual temperature range (warmest minus coldest
//! month) and annual precipitation.
//!
//! The model is deliberately physical in structure and cheap to run:
//!
//! 1. **Zonal energy balance** (Budyko–Sellers): absorbed starlight by latitude against
//!    outgoing longwave `A + B·T` and poleward heat transport `C·(T − T̄)`, with an ice
//!    albedo feedback. Calibrated so an Earth-like body gives Earth's zonal means.
//! 2. **Seasons** from the body's obliquity and eccentricity: the annual harmonic of
//!    insolation sets the seasonal swing, damped by the ocean's heat capacity.
//! 3. **Circulation** from rotation: Hadley, Ferrel and polar cells (the Hadley cell widens
//!    on slowly rotating bodies), shifted seasonally with the ITCZ, plus a monsoon flow
//!    toward heated continents in summer.
//! 4. **Continentality**: a maritime-air tracer carried by the winds, so west coasts in the
//!    westerlies are mild and continental interiors have large seasonal swings.
//! 5. **Moisture**: water vapour carried along the winds on a coarse atmosphere grid;
//!    evaporation over the sea relaxes it toward saturation, it rains out where air rises
//!    (convergence zones, windward slopes) or cools, and not where it sinks or warms (the
//!    subtropical highs, rain shadows, interiors). Land recycles part of its rain.
//!    The global total is set by a water budget scaled with temperature and ocean area.
//! 6. **Downscaling** to the world grid: lapse rate on the cell's own elevation, coast
//!    distance on the fine grid, local orographic enhancement.

use geoforge_core::derived::{daily_insolation, insolation_by_latitude};
use geoforge_core::grid::{distance_field, Cell};
use geoforge_core::product::{Generator, LayerId, Precondition, ProductData, ProductId, Products, RunContext};
use geoforge_core::shapes::Field;
use geoforge_core::{BodyContext, Error};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct ClimateParams {
    /// Atmosphere grid spacing, degrees.
    pub resolution_deg: f64,
    pub lapse_rate_c_per_km: f64,
    /// Extra greenhouse warming (or cooling, negative) on top of the atmosphere's own, °C.
    pub greenhouse_offset_c: f64,
    /// Poleward heat transport coefficient C, W/m²/°C (Earth ≈ 2.6).
    pub heat_transport: f64,
    /// Global mean precipitation, mm/yr. `None` = from temperature and ocean area.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub global_precip_mm: Option<f64>,
    /// How far maritime air reaches inland, km.
    pub maritime_reach_km: f64,
    /// Strength of the summer monsoon flow toward continents (0 = none).
    pub monsoon: f64,
    /// Rain-out rate per atmosphere step (sets how fast air dries inland).
    pub rainout: f64,
    /// Carbonate–silicate thermostat: the fraction of a radiative-forcing difference from
    /// Earth (fainter or brighter star, thinner or thicker air) that weathering has
    /// compensated by adjusting CO₂ over geological time. 0 = no regulation.
    pub weathering_feedback: f64,
}

impl Default for ClimateParams {
    fn default() -> Self {
        Self {
            resolution_deg: 1.0,
            lapse_rate_c_per_km: 6.5,
            greenhouse_offset_c: 0.0,
            heat_transport: 3.2,
            global_precip_mm: None,
            maritime_reach_km: 1200.0,
            monsoon: 0.7,
            rainout: 0.06,
            weathering_feedback: 0.75,
        }
    }
}

// Earth calibration of the energy balance (see `tests::earth_zonal_means`).
const OLR_A: f64 = 200.0;
const OLR_B: f64 = 2.09;

/// Ice-albedo feedback, smooth: 0.30 open, 0.62 ice-covered.
fn albedo(t: f64) -> f64 {
    0.46 - 0.16 * ((t + 8.0) / 10.0).tanh()
}

/// Zonal sea-level mean temperatures (°C) at the given latitudes, and their global mean.
#[must_use]
pub fn zonal_temperature(lats_deg: &[f64], flux: f64, obliquity: f64, e: f64, olr_a: f64, c: f64) -> (Vec<f64>, f64) {
    let s: Vec<f64> = lats_deg.iter().map(|&l| insolation_by_latitude(flux, l, obliquity, e).annual_w_m2).collect();
    let w: Vec<f64> = lats_deg.iter().map(|l| l.to_radians().cos()).collect();
    let wsum: f64 = w.iter().sum();
    let mut t: Vec<f64> = vec![10.0; lats_deg.len()];
    let mut mean = 10.0;
    for _ in 0..200 {
        mean = t.iter().zip(&w).map(|(a, b)| a * b).sum::<f64>() / wsum;
        let next: Vec<f64> = s.iter().zip(&t).map(|(&si, &ti)| (si * (1.0 - albedo(ti)) - olr_a + c * mean) / (OLR_B + c)).collect();
        let d = next.iter().zip(&t).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
        t = next.iter().zip(&t).map(|(a, b)| 0.5 * a + 0.5 * b).collect();
        if d < 1e-4 {
            break;
        }
    }
    (t, mean)
}

/// Peak-to-peak annual harmonic of daily insolation at a latitude (W/m²), positive when the
/// peak is in northern summer.
#[must_use]
pub fn seasonal_swing(flux: f64, lat_deg: f64, obliquity: f64, e: f64) -> f64 {
    let lat = lat_deg.to_radians().clamp(-PI / 2.0 + 1e-6, PI / 2.0 - 1e-6);
    let eps = obliquity.to_radians();
    let n = 48;
    let (mut re, mut im) = (0.0, 0.0);
    for i in 0..n {
        let ls = 2.0 * PI * (f64::from(i) + 0.5) / f64::from(n);
        let decl = (eps.sin() * ls.sin()).asin();
        // Distance factor: perihelion near northern winter, as on Earth.
        let s = flux * (1.0 + 2.0 * e * (ls + PI / 2.0 + 0.2).cos());
        let q = daily_insolation(s, lat, decl);
        re += q * ls.sin();
        im += q * ls.cos();
    }
    let amp = 2.0 * (re * re + im * im).sqrt() / f64::from(n);
    2.0 * amp * re.signum()
}

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

/// Circulation of one season.
#[derive(Debug, Clone, Copy)]
pub struct Circulation {
    /// Latitude of the ITCZ, degrees.
    pub itcz: f64,
    /// Width of each Hadley cell, degrees from the ITCZ.
    pub hadley: f64,
    /// Poleward edge of the Ferrel cells, degrees.
    pub ferrel: f64,
}

impl Circulation {
    #[must_use]
    pub fn for_body(rotation_h: f64, obliquity: f64, season: f64) -> Self {
        // Hadley width ~ Ω^(-2/3) (GCM scaling); Earth: 30°.
        let hadley = (30.0 * (rotation_h.abs().max(1.0) / 24.0).powf(2.0 / 3.0)).clamp(12.0, 80.0);
        let ferrel = (hadley * 2.0).min(85.0);
        Self { itcz: season * 0.3 * obliquity.min(45.0), hadley, ferrel }
    }

    /// Surface wind (east, north) and speed factor at a latitude.
    #[must_use]
    pub fn wind(&self, lat: f64) -> (f64, f64) {
        let rel = lat - self.itcz;
        let side = if rel >= 0.0 { 1.0 } else { -1.0 };
        let hemi = if lat >= 0.0 { 1.0 } else { -1.0 };
        let x = rel.abs() / self.hadley;
        if x < 1.0 {
            // Trades: easterly, blowing toward the ITCZ; calm at the ITCZ and at the high.
            let s = (PI * x).sin().max(0.0);
            return (-0.9 * s, -0.45 * s * side);
        }
        let edge = self.itcz + side * self.hadley;
        let ferrel_w = (self.ferrel - edge.abs()).max(5.0);
        let x2 = (lat.abs() - edge.abs()) / ferrel_w;
        if x2 < 1.0 {
            let s = (PI * x2.max(0.0)).sin();
            return (1.0 * s, 0.3 * s * hemi);
        }
        let x3 = ((lat.abs() - self.ferrel) / (90.0 - self.ferrel).max(1.0)).clamp(0.0, 1.0);
        let s = (PI * x3).sin();
        (-0.6 * s, -0.35 * s * hemi)
    }

    /// Rain-out multiplier from large-scale vertical motion.
    #[must_use]
    pub fn convergence(&self, lat: f64) -> f64 {
        let rel = (lat - self.itcz).abs();
        let itcz = 0.9 * (-(rel / 7.0).powi(2)).exp();
        let high = 0.65 * (-((rel - self.hadley) / 8.0).powi(2)).exp();
        // Mid-latitude storm tracks, and frequent frontal precipitation poleward of them.
        let storm_lat = (self.hadley + 18.0).min(self.ferrel);
        let storm = 0.8 * (-((lat.abs() - storm_lat) / 14.0).powi(2)).exp();
        let subpolar = if lat.abs() > storm_lat { 0.5 } else { 0.0 };
        (1.0 + itcz - high + storm + subpolar).max(0.2)
    }
}

/// The `A` of outgoing longwave `A + B·T` for a body: pressure broadening (≈5 °C per
/// e-fold of surface pressure), the user's offset, and the weathering thermostat's partial
/// compensation of the forcing difference from Earth.
#[must_use]
pub fn outgoing_a(flux: f64, pressure_bar: f64, p: &ClimateParams) -> f64 {
    let pressure = pressure_bar.max(0.001);
    let forcing_starlight = (flux - 1361.0) / 4.0 * (1.0 - 0.3);
    let forcing_air = OLR_B * 5.0 * (pressure / 1.013).ln();
    let regulated = (1.0 - p.weathering_feedback.clamp(0.0, 1.0)) * (forcing_starlight + forcing_air) - forcing_starlight;
    OLR_A - regulated - OLR_B * p.greenhouse_offset_c
}

/// Saturation column water (precipitable water at saturation), mm, from temperature.
#[must_use]
pub fn w_sat(t: f64) -> f64 {
    (55.0 * (0.064 * (t - 27.0)).exp()).clamp(0.2, 90.0)
}

pub struct Climate;

/// Coarse fields the moisture model works on.
struct Atmo {
    g: AtmoGrid,
    land: Vec<f64>,
    elev: Vec<f64>,
}

fn semi_lagrangian(g: &AtmoGrid, wind: &[(f64, f64)], f: &[f64], speed_floor: f64) -> Vec<f64> {
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

fn diffuse(g: &AtmoGrid, f: &mut [f64], k: f64) {
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

/// Winds for one season: cell pattern plus monsoon flow toward (summer) or away from
/// (winter) large heated land masses.
fn winds(a: &Atmo, circ: &Circulation, season: f64, monsoon: f64) -> Vec<(f64, f64)> {
    let g = &a.g;
    let rad = (12.0 / g.res).round().max(1.0) as usize;
    let ls = g.blur(&a.land, rad);
    (0..g.len())
        .map(|i| {
            let (x, y) = (i % g.w, i / g.w);
            let lat = g.lat(y);
            let (mut u, mut v) = circ.wind(lat);
            if monsoon > 0.0 {
                let c = lat.to_radians().cos().max(0.1);
                let xe = (x + 1) % g.w;
                let xw = (x + g.w - 1) % g.w;
                let dx = (ls[y * g.w + xe] - ls[y * g.w + xw]) / (2.0 * g.res * c);
                let yn = y.saturating_sub(1);
                let ys = (y + 1).min(g.h - 1);
                let dy = (ls[yn * g.w + x] - ls[ys * g.w + x]) / ((ys - yn).max(1) as f64 * g.res);
                let mag = (dx * dx + dy * dy).sqrt();
                if mag > 1e-9 {
                    let summer = lat * season > 0.0;
                    let taper = 1.0 - ((lat.abs() - 25.0) / 20.0).clamp(0.0, 1.0);
                    let strength = monsoon * taper * (mag / 0.02).min(1.0) * if summer { 1.0 } else { -0.4 };
                    u += strength * dx / mag;
                    v += strength * dy / mag;
                }
            }
            (u, v)
        })
        .collect()
}

/// Maritime-air fraction carried by the winds (1 over open sea, decaying inland).
fn maritime(a: &Atmo, wind: &[(f64, f64)], reach_km: f64, radius_km: f64) -> Vec<f64> {
    let g = &a.g;
    let step_km = g.res.to_radians() * radius_km;
    let decay = (-step_km / reach_km).exp();
    let mut m: Vec<f64> = a.land.iter().map(|l| 1.0 - l).collect();
    for _ in 0..(3.0 * reach_km / step_km).ceil() as usize + 20 {
        let up = semi_lagrangian(g, wind, &m, 0.35);
        m = up
            .iter()
            .zip(&a.land)
            .map(|(&mu, &l)| {
                let over_land = mu * decay;
                let over_sea = mu + (1.0 - mu) * 0.25;
                l * over_land + (1.0 - l) * over_sea
            })
            .collect();
        diffuse(g, &mut m, 0.15);
    }
    m
}

/// One season of the moisture model; returns raw precipitation (mm per step).
#[allow(clippy::too_many_arguments)]
fn moisture(a: &Atmo, wind: &[(f64, f64)], circ: &Circulation, temp: &[f64], temp_mean: &[f64], rainout: f64) -> Vec<f64> {
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

/// Everything the climate layer computes on the atmosphere grid.
pub struct Coarse {
    pub grid: AtmoGrid,
    /// Sea-level-equivalent mean temperature, °C.
    pub t_sl: Vec<f64>,
    /// Seasonal swing coefficient: range = swing × (0.13 + 0.87 × continentality).
    pub swing: Vec<f64>,
    pub maritime: Vec<f64>,
    pub precip: Vec<f64>,
    pub elev: Vec<f64>,
    pub global_mean_t: f64,
}

fn continentality(is_land: bool, maritime: f64, coast_km: f64) -> f64 {
    if is_land {
        let coastal = 0.55 * (-coast_km / 250.0).exp();
        (1.0 - maritime.max(coastal)).clamp(0.0, 1.0)
    } else {
        (0.5 * (1.0 - maritime)).clamp(0.0, 0.5)
    }
}

impl Climate {
    /// Run the atmosphere on the coarse grid.
    fn coarse(ctx: &RunContext<'_>, elev: &[f32], p: &ClimateParams) -> Coarse {
        let grid = ctx.grid;
        let body = ctx.body;
        let g = AtmoGrid::new(p.resolution_deg);
        let r = ctx.radius_km();
        // Aggregate the world grid onto the atmosphere grid.
        let mut sum_e = vec![0.0; g.len()];
        let mut sum_l = vec![0.0; g.len()];
        let mut sum_a = vec![0.0; g.len()];
        for c in grid.cells() {
            let ll = grid.center(c).to_lat_lon();
            let i = g.index(ll.lat_deg, ll.lon_deg);
            let a = grid.area_sr(c);
            let e = f64::from(elev[c as usize]);
            sum_e[i] += e.max(0.0) * a;
            sum_l[i] += if e >= 0.0 { a } else { 0.0 };
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
                let c = grid.locate(geoforge_core::UnitVec3::from_lat_lon(g.lat(y), g.lon(x)));
                let e = f64::from(elev[c as usize]);
                land[i] = if e >= 0.0 { 1.0 } else { 0.0 };
                elev_c[i] = e.max(0.0);
            }
        }
        let atmo = Atmo { g, land, elev: elev_c };

        // 1. Zonal energy balance.
        let lats: Vec<f64> = (0..g.h).map(|y| g.lat(y)).collect();
        let flux = body.irradiance.flux_w_m2;
        let obl = body.exterior.obliquity_deg;
        let ecc = body.orbit.e;
        let olr_a = outgoing_a(flux, body.interior.surface_pressure_bar(), p);
        let (tz, tmean) = zonal_temperature(&lats, flux, obl, ecc, olr_a, p.heat_transport);
        let swing_z: Vec<f64> = lats.iter().map(|&l| seasonal_swing(flux, l, obl, ecc)).collect();

        // 2–4. Circulation and continentality from annual-mean winds.
        let annual = Circulation::for_body(body.exterior.rotation_h, obl, 0.0);
        let wind_annual = winds(&atmo, &annual, 0.0, 0.0);
        let mar = maritime(&atmo, &wind_annual, p.maritime_reach_km, r);
        let mut t_sl = vec![0.0; g.len()];
        let mut swing = vec![0.0; g.len()];
        let mut t_season = [vec![0.0; g.len()], vec![0.0; g.len()]];
        let mut t_mean_here = vec![0.0; g.len()];
        for i in 0..g.len() {
            let y = i / g.w;
            let lat = lats[y];
            let is_land = atmo.land[i] >= 0.5;
            let cont = continentality(is_land, mar[i], if is_land { 400.0 } else { 0.0 });
            // Interiors at higher latitude are colder than the zonal mean; mild west coasts warmer.
            let hi = (lat.to_radians().sin().abs() - 0.3).max(0.0);
            let t = tz[y] - 7.0 * (cont - 0.35) * hi * atmo.land[i].max(0.3);
            t_sl[i] = t;
            swing[i] = swing_z[y] * 0.10;
            let range = swing[i].abs() * (0.13 + 0.87 * cont);
            let half = 0.5 * range * swing[i].signum();
            let t_here = t - p.lapse_rate_c_per_km * atmo.elev[i] / 1000.0;
            t_mean_here[i] = t_here;
            t_season[0][i] = t_here + half; // northern summer
            t_season[1][i] = t_here - half; // northern winter
        }

        // 5. Moisture, two solstice seasons.
        let mut precip = vec![0.0; g.len()];
        for (k, season) in [1.0, -1.0].into_iter().enumerate() {
            let circ = Circulation::for_body(body.exterior.rotation_h, obl, season);
            let w = winds(&atmo, &circ, season, p.monsoon);
            let ps = moisture(&atmo, &w, &circ, &t_season[k], &t_mean_here, p.rainout);
            for i in 0..g.len() {
                precip[i] += 0.5 * ps[i];
            }
        }
        Coarse { grid: g, t_sl, swing, maritime: mar, precip, elev: atmo.elev, global_mean_t: tmean }
    }
}

impl Generator for Climate {
    type Params = ClimateParams;
    const LAYER: LayerId = LayerId::Climate;

    fn needs() -> &'static [ProductId] {
        &[ProductId::Elevation]
    }
    fn applies(body: &BodyContext) -> Result<(), Precondition> {
        if !body.interior.has_solid_surface() {
            return Err(Precondition("climate (surface) needs a solid surface".into()));
        }
        if !body.interior.has_atmosphere() {
            return Err(Precondition("climate needs an atmosphere".into()));
        }
        Ok(())
    }

    fn run(ctx: &RunContext<'_>, inputs: &Products, p: &ClimateParams) -> Result<Products, Error> {
        let mut out = Products::new();
        if inputs.temperature().is_some() && inputs.temperature_range().is_some() && inputs.precipitation().is_some() {
            return Ok(out);
        }
        let elev = &inputs.elevation().ok_or(Error::MissingInput(ProductId::Elevation))?.values;
        let grid = ctx.grid;
        let r = ctx.radius_km();
        let co = Self::coarse(ctx, elev, p);
        let g = co.grid;

        // Fine coast distance (km), land side and sea side.
        let coast: Vec<Cell> = grid
            .cells()
            .filter(|&c| {
                let l = elev[c as usize] >= 0.0;
                grid.neighbors(c).iter().any(|n| (elev[n as usize] >= 0.0) != l)
            })
            .collect();
        let coast_km: Vec<f64> = distance_field(grid, coast, 1500.0 / r, |_| true).into_iter().map(|d| (d * r).min(1500.0)).collect();
        let e_smooth = g.blur(&co.elev, 1);

        let n = grid.len();
        // Stationary waves and regional weather regimes: smooth anomalies so climate zones
        // meander instead of following parallels (±~1 °C, ±~25% precipitation).
        let wave_t = crate::noise::TerrainNoise::new(ctx.seed(ProductId::Temperature).value());
        let wave_p = crate::noise::TerrainNoise::new(ctx.seed(ProductId::Precipitation).value());
        let fine: Vec<(f32, f32, f32)> = (0..n)
            .into_par_iter()
            .map(|i| {
                let pt = grid.center(i as Cell);
                let ll = pt.to_lat_lon();
                let e = f64::from(elev[i]);
                let is_land = e >= 0.0;
                let mar = g.sample(&co.maritime, ll.lat_deg, ll.lon_deg);
                let cont = continentality(is_land, mar, coast_km[i]);
                let t_sl = g.sample(&co.t_sl, ll.lat_deg, ll.lon_deg) + 3.5 * wave_t.fbm(pt, 3.0, 3);
                let t = t_sl - p.lapse_rate_c_per_km * e.max(0.0) / 1000.0;
                let swing = g.sample(&co.swing, ll.lat_deg, ll.lon_deg).abs();
                let range = swing * (0.13 + 0.87 * cont);
                let mut pr = g.sample(&co.precip, ll.lat_deg, ll.lon_deg).max(0.0) * (0.8 * wave_p.fbm(pt, 4.0, 3)).exp();
                if is_land {
                    // Local relief above the atmosphere grid's smooth surface catches more rain.
                    let rel = (e - g.sample(&e_smooth, ll.lat_deg, ll.lon_deg)) / 1000.0;
                    pr *= (1.0 + 0.45 * rel).clamp(0.55, 1.9);
                }
                (t as f32, range as f32, pr as f32)
            })
            .collect();
        // Dry land heats more (little evaporative cooling, clear skies): warmer means and
        // larger swings in deserts. Uses the normalised precipitation below.
        let arid = |pr_mm: f64| (1.0 - pr_mm / 900.0).clamp(0.0, 1.0);

        // Water budget: global mean precipitation from temperature and ocean area.
        let area: Vec<f64> = (0..n).map(|i| grid.area_sr(i as Cell)).collect();
        let total_a: f64 = area.iter().sum();
        let ocean_frac = (0..n).filter(|&i| elev[i] < 0.0).map(|i| area[i]).sum::<f64>() / total_a;
        let raw_mean = fine.iter().zip(&area).map(|(f, a)| f64::from(f.2) * a).sum::<f64>() / total_a;
        let target = p.global_precip_mm.unwrap_or_else(|| 1000.0 * 1.03f64.powf(co.global_mean_t - 14.5) * (ocean_frac / 0.71).max(0.0).powf(0.6));
        let k = if raw_mean > 0.0 { target / raw_mean } else { 0.0 };

        let spec = grid.spec();
        let dry_warm = |i: usize| -> (f32, f32) {
            if elev[i] < 0.0 {
                return (0.0, 0.0);
            }
            let a = arid(f64::from(fine[i].2) * k);
            let lat = grid.center(i as Cell).to_lat_lon().lat_deg.abs();
            let sub = (1.0 - ((lat - 25.0) / 30.0).powi(2)).clamp(0.0, 1.0);
            ((3.5 * a * sub) as f32, (5.0 * a * sub) as f32)
        };
        if inputs.temperature().is_none() {
            out.insert(ProductData::Temperature(Field::from_values(spec, (0..n).map(|i| fine[i].0 + dry_warm(i).0).collect())));
        }
        if inputs.temperature_range().is_none() {
            out.insert(ProductData::TemperatureRange(Field::from_values(spec, (0..n).map(|i| fine[i].1 + dry_warm(i).1).collect())));
        }
        if inputs.precipitation().is_none() {
            out.insert(ProductData::Precipitation(Field::from_values(spec, fine.iter().map(|f| (f64::from(f.2) * k) as f32).collect())));
        }
        Ok(out)
    }
}

/// Potential evapotranspiration, mm/yr, from mean annual temperature (Thornthwaite-like,
/// calibrated to ~1500 mm at 25 °C and ~170 mm at 0 °C).
#[must_use]
pub fn pet_mm(t: f64) -> f64 {
    if t <= -5.0 { 0.0 } else { 25.0 * (t + 5.0).powf(1.2) }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn earth_zonal_means() {
        let lats: Vec<f64> = (0..180).map(|i| 89.5 - f64::from(i)).collect();
        let (t, mean) = zonal_temperature(&lats, 1361.0, 23.44, 0.0167, OLR_A, 3.0);
        let at = |l: f64| t[(89.5 - l) as usize];
        // Earth (sea-level, zonal, annual): ~14–15 °C global, ~26 °C equator, ~0 °C at 60°,
        // ~-20 °C near the poles.
        assert!((13.5..17.0).contains(&mean), "global {mean}");
        assert!((25.0..28.5).contains(&at(0.5)), "equator {}", at(0.5));
        assert!((8.0..13.0).contains(&at(45.5)), "45N {}", at(45.5));
        assert!((-5.0..3.0).contains(&at(60.5)), "60N {}", at(60.5));
        assert!((-22.0..-10.0).contains(&at(89.5)), "pole {}", at(89.5));
    }

    #[test]
    fn seasons_follow_latitude() {
        let s0 = seasonal_swing(1361.0, 0.0, 23.44, 0.0).abs();
        let s60 = seasonal_swing(1361.0, 60.0, 23.44, 0.0);
        let s60s = seasonal_swing(1361.0, -60.0, 23.44, 0.0);
        assert!(s0 < 60.0, "{s0}");
        assert!(s60 > 380.0 && s60s < -380.0, "{s60} {s60s}");
        // No obliquity, no seasons.
        assert!(seasonal_swing(1361.0, 45.0, 0.0, 0.0).abs() < 1.0);
    }

    #[test]
    fn hadley_cell_scales_with_rotation() {
        let earth = Circulation::for_body(24.0, 23.0, 0.0);
        let slow = Circulation::for_body(240.0, 23.0, 0.0);
        assert!((earth.hadley - 30.0).abs() < 1e-9);
        assert!(slow.hadley > 60.0);
        // Trades in the tropics, westerlies in mid-latitudes.
        assert!(earth.wind(15.0).0 < -0.5 && earth.wind(-15.0).0 < -0.5);
        assert!(earth.wind(45.0).0 > 0.5 && earth.wind(-45.0).0 > 0.5);
    }
}
