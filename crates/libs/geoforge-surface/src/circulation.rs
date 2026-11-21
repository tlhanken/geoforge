//! Circulation layer: seasonal surface winds, ocean surface currents and the sea-surface
//! temperature anomaly they cause.
//!
//! - **Winds** (January and July): Hadley, Ferrel and polar cells scaled by rotation and
//!   shifted with the ITCZ (further over large heated land, where the trades cross the
//!   equator as monsoon westerlies); monsoon flow into the summer heat lows (and out of the
//!   winter highs), turned by the Coriolis force and blocked by high ranges; and the
//!   subtropical highs over each ocean basin, which turn the flow poleward
//!   along the basin's western side (moist air onto east coasts) and equatorward along its
//!   eastern side (dry west coasts).
//! - **Currents**: wind-driven gyres from the Sverdrup balance (transport set by the curl of
//!   the zonal wind stress, integrated westward from each basin's eastern shore) closed by a
//!   narrow western boundary current (Stommel), which crosses the basin where the gyres
//!   meet and runs up the far shore (the North Atlantic Drift); plus equatorial currents and
//!   a circumpolar current where an ocean band has no shores.
//! - **SST anomaly**: currents carry water across the temperature gradient (poleward
//!   western boundary currents run warm, equatorward eastern ones cold), and that heat is
//!   carried downstream along the current; coastal upwelling under equatorward winds cools
//!   eastern margins further.
//!
//! Climate reads these products, so editing a wind or current map changes rainfall and
//! temperature downstream.

use crate::atmosphere::transport::semi_lagrangian;
use crate::atmosphere::grid::is_sea;
use crate::atmosphere::{winds, Atmo, AtmoGrid, Radiation, Season};
use crate::climate::ClimateParams;
use geoforge_core::grid::Cell;
use geoforge_core::product::{Generator, LayerId, Precondition, ProductData, ProductId, Products, RunContext};
use geoforge_core::shapes::Field;
use geoforge_core::{BodyContext, Error};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// Widest a single gyre gets (about the Pacific's width), km.
const MAX_GYRE_KM: f64 = 12_000.0;

/// Wind speed that the moisture model treats as one atmosphere cell per step.
pub const WIND_REF_M_S: f64 = 8.0;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct CirculationParams {
    /// Atmosphere/ocean grid spacing, degrees.
    pub resolution_deg: f64,
    /// Strength of the summer monsoon flow toward continents (0 = none).
    pub monsoon: f64,
    /// Strength of the subtropical highs' meridional flow, relative to the trades.
    pub subtropical_highs: f64,
    /// Peak speed of western boundary currents, m/s.
    pub boundary_current_m_s: f64,
    /// Width of western boundary currents, km.
    pub boundary_current_km: f64,
    /// Speed of the drift that carries a western boundary current across its basin and up
    /// the far shore (the North Atlantic Drift), m/s.
    pub drift_m_s: f64,
    /// How long current-carried water keeps its temperature, days (sets anomaly size).
    pub heat_memory_days: f64,
    /// Cooling where winds drive coastal upwelling, °C.
    pub upwelling_c: f64,
}

impl Default for CirculationParams {
    fn default() -> Self {
        Self {
            resolution_deg: 1.0,
            monsoon: 0.7,
            subtropical_highs: 0.5,
            boundary_current_m_s: 1.5,
            boundary_current_km: 150.0,
            drift_m_s: 0.3,
            heat_memory_days: 60.0,
            upwelling_c: 12.0,
        }
    }
}

/// Ocean basins along each latitude row: for every sea cell, distance (km) to the shore to
/// its west and to its east, or `None` for a row of open sea all the way round.
struct Basins {
    west_km: Vec<f64>,
    east_km: Vec<f64>,
    open_row: Vec<bool>,
}

fn basins(a: &Atmo, radius_km: f64) -> Basins {
    let g = &a.g;
    let mut west_km = vec![0.0; g.len()];
    let mut east_km = vec![0.0; g.len()];
    let mut open_row = vec![false; g.h];
    for y in 0..g.h {
        let dx = g.res.to_radians() * radius_km * g.lat(y).to_radians().cos().max(0.02);
        let sea = |x: usize| a.land[y * g.w + x] < 0.5;
        if (0..g.w).all(sea) {
            // No shore anywhere round the circle.
            open_row[y] = true;
            for x in 0..g.w {
                west_km[y * g.w + x] = f64::INFINITY;
                east_km[y * g.w + x] = f64::INFINITY;
            }
            continue;
        }
        // Walk east from a land cell, then west from a land cell, counting sea cells.
        let start = (0..g.w).find(|&x| !sea(x)).expect("row has land");
        let mut run = 0usize;
        for k in 1..=g.w {
            let x = (start + k) % g.w;
            if sea(x) {
                run += 1;
                west_km[y * g.w + x] = (run as f64 - 0.5) * dx;
            } else {
                run = 0;
            }
        }
        run = 0;
        for k in 1..=g.w {
            let x = (start + g.w - k) % g.w;
            if sea(x) {
                run += 1;
                east_km[y * g.w + x] = (run as f64 - 0.5) * dx;
            } else {
                run = 0;
            }
        }
    }
    Basins { west_km, east_km, open_row }
}

/// Surface winds of one season on the atmosphere grid, in units of [`WIND_REF_M_S`].
#[must_use]
pub fn season_winds(a: &Atmo, season: &Season, p: &CirculationParams) -> Vec<(f64, f64)> {
    let g = &a.g;
    let cells = &season.cells;
    let mut w = winds(a, season, p.monsoon);
    if p.subtropical_highs > 0.0 {
        // Anticyclones over each basin: poleward on the western side, equatorward on the east.
        let b = basins(a, 6371.0);
        let mut v_high = vec![0.0; g.len()];
        for i in 0..g.len() {
            let y = i / g.w;
            if a.land[i] >= 0.5 || b.open_row[y] {
                continue;
            }
            let lat = g.lat(y);
            let hemi = lat.signum();
            let (we, ea) = (b.west_km[i], b.east_km[i]);
            let x_rel = we / (we + ea).max(1.0); // 0 at the western shore, 1 at the eastern
            let edge = (cells.itcz + hemi * cells.hadley).abs();
            let band = (-((lat.abs() - edge) / 11.0).powi(2)).exp();
            v_high[i] = hemi * p.subtropical_highs * band * (std::f64::consts::PI * x_rel).cos();
        }
        // Let the flow reach onto the coasts, and round off the basins' row-by-row edges.
        let rad = (5.0 / g.res).round().max(1.0) as usize;
        let v_high = g.blur(&v_high, rad);
        for (wi, vh) in w.iter_mut().zip(&v_high) {
            wi.1 += vh;
        }
    }
    w
}

/// Ocean surface currents (m/s, east and north) on the atmosphere grid.
fn currents(a: &Atmo, annual: &[(f64, f64)], radius_km: f64, p: &CirculationParams) -> Vec<(f64, f64)> {
    let g = &a.g;
    let b = basins(a, radius_km);
    // Zonal wind stress (quadratic drag) averaged over each row's sea.
    let tau: Vec<f64> = (0..g.h)
        .map(|y| {
            let (mut s, mut n) = (0.0, 0.0);
            for x in 0..g.w {
                let i = y * g.w + x;
                if a.land[i] < 0.5 {
                    s += annual[i].0 * annual[i].0.abs();
                    n += 1.0;
                }
            }
            if n > 0.0 { s / n } else { 0.0 }
        })
        .collect();
    // Basin-scale stress varies smoothly with latitude (eddies and storms average out the
    // idealised cell pattern's kinks).
    let sigma = (10.0 / g.res).max(1.0);
    let rad = (3.0 * sigma).ceil() as i64;
    let tau: Vec<f64> = (0..g.h as i64)
        .map(|y| {
            let (mut s, mut wsum) = (0.0, 0.0);
            for k in -rad..=rad {
                let yy = (y + k).clamp(0, g.h as i64 - 1) as usize;
                let w = (-0.5 * (k as f64 / sigma).powi(2)).exp();
                s += tau[yy] * w;
                wsum += w;
            }
            s / wsum
        })
        .collect();
    // Curl of the stress over the Coriolis gradient β ∝ cos(lat).
    let c: Vec<f64> = (0..g.h)
        .map(|y| {
            let (yn, ys) = (y.saturating_sub(1), (y + 1).min(g.h - 1));
            let curl = -(tau[yn] - tau[ys]) / ((ys - yn).max(1) as f64 * g.res);
            curl / g.lat(y).to_radians().cos().max(0.05)
        })
        .collect();
    let sig = (5.0 / g.res).max(1.0);
    let open_band: Vec<f64> = (0..g.h as i64)
        .map(|y| {
            let r = (2.0 * sig).ceil() as i64;
            let (mut s, mut w) = (0.0, 0.0);
            for k in -r..=r {
                let yy = (y + k).clamp(0, g.h as i64 - 1) as usize;
                let wt = (-0.5 * (k as f64 / sig).powi(2)).exp();
                s += wt * if b.open_row[yy] { 1.0 } else { 0.0 };
                w += wt;
            }
            s / w
        })
        .collect();
    // Sverdrup interior closed by a Stommel western boundary layer.
    let psi: Vec<f64> = (0..g.len())
        .map(|i| {
            let y = i / g.w;
            if a.land[i] >= 0.5 || b.open_row[y] || g.lat(y).abs() < 8.0 {
                return 0.0;
            }
            // Gyres span at most an ocean's width, and fade out next to a circumpolar band.
            let east = b.east_km[i].min(MAX_GYRE_KM);
            -c[y] * east * (1.0 - (-b.west_km[i] / p.boundary_current_km).exp()) * (1.0 - open_band[y])
        })
        .collect();
    // Meridional gyre flow v = ∂ψ/∂x: fast and narrow along western shores, slow and broad
    // across the interior. (Zonal gyre flow is left to the wind drift below; differentiating
    // ψ across rows with different coastlines gives spurious jets.)
    let mut cur: Vec<(f64, f64)> = (0..g.len())
        .map(|i| {
            let (x, y) = (i % g.w, i / g.w);
            if a.land[i] >= 0.5 {
                return (0.0, 0.0);
            }
            let dx = g.res.to_radians() * radius_km * g.lat(y).to_radians().cos().max(0.02);
            let (xe, xw) = ((x + 1) % g.w, (x + g.w - 1) % g.w);
            (0.0, (psi[y * g.w + xe] - psi[y * g.w + xw]) / (2.0 * dx))
        })
        .collect();
    let fronts = gyre_fronts(g, &tau);
    // Scale so the strongest western boundary currents (meridional flow) reach the set speed.
    let mut speeds: Vec<f64> = cur.iter().map(|v| v.1.abs()).filter(|s| *s > 0.0).collect();
    speeds.sort_by(f64::total_cmp);
    let top = speeds.get(speeds.len().saturating_sub(1 + speeds.len() / 500)).copied().unwrap_or(0.0);
    let k = if top > 0.0 { p.boundary_current_m_s / top } else { 0.0 };
    for (i, v) in cur.iter_mut().enumerate() {
        v.0 *= k;
        v.1 *= k;
        let y = i / g.w;
        if a.land[i] >= 0.5 {
            continue;
        }
        let lat = g.lat(y);
        let drift = wind_drift(annual[i], lat);
        v.0 += drift.0;
        v.1 += drift.1;
        // Equatorial currents are faster than drift alone; a countercurrent just north of the
        // equator.
        let eq = (-(lat / 7.0).powi(2)).exp();
        let counter = (-((lat - 6.0) / 2.5).powi(2)).exp();
        v.0 += 0.3 * annual[i].0.min(0.0) * eq + 0.35 * counter;
        // A circumpolar current where the ocean band has no shores (spread over the band).
        if open_band[y] > 0.0 {
            v.0 += 0.6 * open_band[y] * annual[i].0.max(0.0);
        }
        // Where the subtropical and subpolar gyres meet (under the strongest westerlies) the
        // western boundary current leaves the coast and crosses the basin, turning poleward
        // as it nears the far side and running on up the eastern shore (the Gulf Stream, North
        // Atlantic Drift and Norwegian Current; the Kuroshio, North Pacific and Alaska
        // Currents). Only in a basin of a gyre's width: across a wider stretch of open water
        // (the Atlantic and Indian Oceans together at 45°S) the current joins the circumpolar
        // flow instead.
        let hemi = lat.signum();
        let front = if hemi > 0.0 { fronts.0 } else { fronts.1 };
        if let Some(front) = front {
            let basin = 1.0 - open_band[y];
            let width = b.west_km[i] + b.east_km[i];
            let crossing = 1.0 - smoothstep(MAX_GYRE_KM, 1.5 * MAX_GYRE_KM, width);
            let x_rel = if b.open_row[y] { 0.0 } else { b.west_km[i] / width.max(1.0) };
            let jet = (-((lat.abs() - front) / 5.0).powi(2)).exp();
            let beyond = ((lat.abs() - front) / 5.0).clamp(0.0, 1.0) * ((front + 28.0 - lat.abs()) / 6.0).clamp(0.0, 1.0);
            // An eastern shore the current can run up: a coast that also reaches ~5° toward
            // the equator (Norway, Alaska, southern Chile), not Antarctica's, which runs
            // east–west.
            let (x, dx) = (i % g.w, g.res.to_radians() * radius_km * lat.to_radians().cos().max(0.02));
            let y_eq = (y as i64 + (hemi * 5.0 / g.res).round() as i64).clamp(0, g.h as i64 - 1) as usize;
            let runs_equatorward = !b.open_row[y] && a.land[y_eq * g.w + (x + (b.east_km[i] / dx).round() as usize + 1) % g.w] >= 0.5;
            let shore = if runs_equatorward { (-b.east_km[i] / 700.0).exp() } else { 0.0 };
            v.0 += basin * crossing * p.drift_m_s * jet * (1.0 - 0.5 * x_rel);
            v.1 += basin * p.drift_m_s * hemi * (0.6 * crossing * jet * x_rel + 0.8 * beyond * shore);
        }
    }
    cur
}

fn smoothstep(e0: f64, e1: f64, x: f64) -> f64 {
    let t = ((x - e0) / (e1 - e0)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

/// Wind drift (m/s): ~3% of the wind speed, turned 30° to the right of the wind in the north
/// and to the left in the south (Ekman).
fn wind_drift(wind: (f64, f64), lat: f64) -> (f64, f64) {
    let (wu, wv) = (wind.0 * WIND_REF_M_S * 0.03, wind.1 * WIND_REF_M_S * 0.03);
    let turn = -lat.signum() * 30f64.to_radians();
    (wu * turn.cos() - wv * turn.sin(), wu * turn.sin() + wv * turn.cos())
}

/// Latitude of the strongest westerly stress in each hemisphere (north, south), between 25°
/// and 65°: where the subtropical and subpolar gyres meet.
fn gyre_fronts(g: &AtmoGrid, tau: &[f64]) -> (Option<f64>, Option<f64>) {
    let best = |north: bool| {
        (0..g.h)
            .filter(|&y| {
                let lat = g.lat(y);
                (lat > 0.0) == north && (25.0..65.0).contains(&lat.abs())
            })
            .filter(|&y| tau[y] > 0.0)
            .max_by(|&a, &b| tau[a].total_cmp(&tau[b]))
            .map(|y| g.lat(y).abs())
    };
    (best(true), best(false))
}

/// How far offshore coastal upwelling cools the sea, km (e-folding).
const UPWELLING_KM: f64 = 400.0;
/// Freezing point of seawater, °C.
const SEA_FREEZE_C: f64 = -1.8;
/// Latitude span over which the sea-surface anomaly's zonal mean is taken off, degrees.
const ZONAL_MEAN_SMOOTH_DEG: f64 = 10.0;

/// Sea-surface temperature anomaly (°C) from currents crossing the temperature gradient,
/// carried downstream, plus coastal upwelling.
fn sst_anomaly(a: &Atmo, cur: &[(f64, f64)], seasons: [&[(f64, f64)]; 2], tz: &[f64], radius_km: f64, p: &CirculationParams) -> Vec<f64> {
    let g = &a.g;
    let b = basins(a, radius_km);
    let dy_km = g.res.to_radians() * radius_km;
    // The wind drift's part that is the same all round a latitude circle (the Ekman drift
    // under the westerlies) moves heat the zonal temperature already carries: counting it
    // painted warm and cold bands right round the globe. Only its departures from the
    // latitude's mean, and the gyres, make anomalies.
    let v_mean: Vec<f64> = (0..g.h)
        .map(|y| {
            let annual = |i: usize| (0.5 * (seasons[0][i].0 + seasons[1][i].0), 0.5 * (seasons[0][i].1 + seasons[1][i].1));
            let sea: Vec<usize> = (y * g.w..(y + 1) * g.w).filter(|&i| a.land[i] < 0.5).collect();
            if sea.is_empty() { 0.0 } else { sea.iter().map(|&i| wind_drift(annual(i), g.lat(y)).1).sum::<f64>() / sea.len() as f64 }
        })
        .collect();
    let source: Vec<f64> = (0..g.len())
        .map(|i| {
            if a.land[i] >= 0.5 {
                return 0.0;
            }
            let y = i / g.w;
            let (yn, ys) = (y.saturating_sub(1), (y + 1).min(g.h - 1));
            // The sea surface follows the zonal temperature down only to freezing: over the ice
            // at high latitudes its gradient is the air's, not the water's.
            let sst = |t: f64| t.max(SEA_FREEZE_C);
            let dtdy = (sst(tz[yn]) - sst(tz[ys])) / ((ys - yn).max(1) as f64 * dy_km * 1000.0); // °C per m northward
            let mut s = -(cur[i].1 - v_mean[y]) * dtdy * p.heat_memory_days * 86_400.0;
            let lat = g.lat(y);
            // Upwelling: equatorward wind along an eastern shore, in the subtropics, all year (a
            // monsoon's reversing winds upwell only half the year, and not off India).
            let equatorward = seasons.iter().all(|w| w[i].1 * lat.signum() < 0.0);
            if equatorward && (8.0..42.0).contains(&lat.abs()) && !b.open_row[y] {
                s -= p.upwelling_c * (-b.east_km[i] / UPWELLING_KM).exp();
            }
            s
        })
        .collect();
    // Carry the anomaly along the currents (one atmosphere cell per step at 0.5 m/s): water
    // keeps its warmth for months, so the Gulf Stream's reaches Norway.
    let flow: Vec<(f64, f64)> = cur.iter().map(|v| (v.0 / 0.5, v.1 / 0.5)).collect();
    let mut anom = carry_along(g, &flow, &source, &a.land, 200, 0.97);
    diffuse_sea(g, &mut anom, &a.land, 0.3);
    // A departure from the zonal temperature averages to zero round each latitude circle.
    // The mean taken off is smoothed over ~10° of latitude so that it cannot itself draw
    // bands where one row's mean differs from the next.
    let rows: Vec<(f64, f64)> = (0..g.h)
        .map(|y| (y * g.w..(y + 1) * g.w).filter(|&i| a.land[i] < 0.5).fold((0.0, 0.0), |(s, n), i| (s + anom[i], n + 1.0)))
        .collect();
    let half = (ZONAL_MEAN_SMOOTH_DEG / 2.0 / g.res).round() as i64;
    let mean: Vec<f64> = (0..g.h as i64)
        .map(|y| {
            let (s, n) = (y - half..=y + half).filter(|k| (0..g.h as i64).contains(k)).fold((0.0, 0.0), |(s, n), k| {
                let w = 1.0 - (k - y).abs() as f64 / (half + 1) as f64;
                (s + w * rows[k as usize].0, n + w * rows[k as usize].1)
            });
            if n > 0.0 { s / n } else { 0.0 }
        })
        .collect();
    for (i, v) in anom.iter_mut().enumerate() {
        if a.land[i] < 0.5 {
            *v -= mean[i / g.w];
        }
    }
    anom.iter().map(|v| v.clamp(-8.0, 8.0)).collect()
}

/// A quantity carried by a flow: each step it is replaced by `(1 − keep)·source +
/// keep·(value upstream)`, so it is a moving average of the sources along the path, with an
/// e-folding length of `1 / (1 − keep)` steps. Zero on land; water whose path starts at a
/// shore (upwelled from below) takes the local source rather than the land's zero, which
/// had washed out the cold along every coast.
fn carry_along(g: &AtmoGrid, flow: &[(f64, f64)], source: &[f64], land: &[f64], steps: usize, keep: f64) -> Vec<f64> {
    let sea: Vec<f64> = land.iter().map(|&l| if l >= 0.5 { 0.0 } else { 1.0 }).collect();
    let sea_up = semi_lagrangian(g, flow, &sea, 0.0);
    let mut a: Vec<f64> = source.iter().zip(&sea).map(|(s, w)| s * w).collect();
    for _ in 0..steps {
        let up = semi_lagrangian(g, flow, &a, 0.0);
        a = (0..g.len())
            .map(|i| {
                if sea[i] == 0.0 {
                    return 0.0;
                }
                let upstream = if sea_up[i] > 0.05 { up[i] / sea_up[i] } else { source[i] };
                (1.0 - keep) * source[i] + keep * upstream
            })
            .collect();
    }
    a
}

/// Diffusion over the sea alone: land neither gives nor takes.
fn diffuse_sea(g: &AtmoGrid, f: &mut [f64], land: &[f64], k: f64) {
    let src = f.to_vec();
    let sea = |i: usize| land[i] < 0.5;
    for y in 0..g.h {
        for x in 0..g.w {
            let i = y * g.w + x;
            if !sea(i) {
                continue;
            }
            let nbrs = [y * g.w + (x + g.w - 1) % g.w, y * g.w + (x + 1) % g.w, y.saturating_sub(1) * g.w + x, (y + 1).min(g.h - 1) * g.w + x];
            let (s, n) = nbrs.iter().filter(|&&j| sea(j)).fold((0.0, 0.0), |(s, n), &j| (s + src[j], n + 1.0));
            if n > 0.0 {
                f[i] = (1.0 - k) * src[i] + k * s / n;
            }
        }
    }
}

pub struct Circulation;

impl Generator for Circulation {
    type Params = CirculationParams;
    const NAME: &'static str = "circulation";
    const VERSION: u32 = 6;
    const LAYER: LayerId = LayerId::Circulation;

    fn needs() -> &'static [ProductId] {
        &[ProductId::Elevation]
    }
    fn applies(body: &BodyContext) -> Result<(), Precondition> {
        crate::preconditions::surface_air(body)
    }

    fn run(ctx: &RunContext<'_>, inputs: &Products, p: &CirculationParams) -> Result<Products, Error> {
        let mut out = Products::new();
        if [ProductId::WindJan, ProductId::WindJul, ProductId::OceanCurrents, ProductId::SstAnomaly].iter().all(|x| inputs.contains(*x)) {
            return Ok(out);
        }
        let elev = &inputs.elevation().ok_or(Error::MissingInput(ProductId::Elevation))?.values;
        let grid = ctx.grid;
        let body = ctx.body;
        let r = ctx.radius_km();
        let wet = body.interior.has_surface_liquid_water();
        let atmo = Atmo::from_world(grid, elev, p.resolution_deg, wet);
        let g: AtmoGrid = atmo.g;
        let obl = body.irradiance.obliquity_deg;
        let rot = body.exterior.rotation_h;
        let jul = season_winds(&atmo, &Season::new(&atmo, rot, obl, 1.0, p.monsoon), p);
        let jan = season_winds(&atmo, &Season::new(&atmo, rot, obl, -1.0, p.monsoon), p);
        let annual: Vec<(f64, f64)> = jul.iter().zip(&jan).map(|(a, b)| (0.5 * (a.0 + b.0), 0.5 * (a.1 + b.1))).collect();
        let has_sea = atmo.land.iter().any(|l| *l < 0.5);
        let cur = if has_sea { currents(&atmo, &annual, r, p) } else { vec![(0.0, 0.0); g.len()] };
        let lats: Vec<f64> = (0..g.h).map(|y| g.lat(y)).collect();
        let cp = ClimateParams::default();
        let (tz, _) = Radiation::for_body(body, cp.heat_transport, cp.greenhouse_offset_c).zonal_temperature(&lats);
        let sst = if has_sea { sst_anomaly(&atmo, &cur, [&jan, &jul], &tz, r, p) } else { vec![0.0; g.len()] };

        // To the world grid.
        let n = grid.len();
        struct Flow {
            jan: [f32; 2],
            jul: [f32; 2],
            current: [f32; 2],
            sst: f32,
        }
        let fine: Vec<Flow> = (0..n)
            .into_par_iter()
            .map(|i| {
                let ll = grid.center(i as Cell).to_lat_lon();
                let s2 = |f: &[(f64, f64)], k: f64| {
                    let (x, y) = sample2(&g, f, ll.lat_deg, ll.lon_deg);
                    [(x * k) as f32, (y * k) as f32]
                };
                let sea = is_sea(wet, f64::from(elev[i]));
                let cu = if sea { s2(&cur, 1.0) } else { [0.0, 0.0] };
                let sa = if sea { g.sample(&sst, ll.lat_deg, ll.lon_deg) as f32 } else { 0.0 };
                Flow { jan: s2(&jan, WIND_REF_M_S), jul: s2(&jul, WIND_REF_M_S), current: cu, sst: sa }
            })
            .collect();
        let spec = grid.spec();
        if !inputs.contains(ProductId::WindJan) {
            out.insert(ProductData::WindJan(Field::from_values(spec, fine.iter().map(|f| f.jan).collect())));
        }
        if !inputs.contains(ProductId::WindJul) {
            out.insert(ProductData::WindJul(Field::from_values(spec, fine.iter().map(|f| f.jul).collect())));
        }
        if !inputs.contains(ProductId::OceanCurrents) {
            out.insert(ProductData::OceanCurrents(Field::from_values(spec, fine.iter().map(|f| f.current).collect())));
        }
        if !inputs.contains(ProductId::SstAnomaly) {
            out.insert(ProductData::SstAnomaly(Field::from_values(spec, fine.iter().map(|f| f.sst).collect())));
        }
        Ok(out)
    }
}

/// Bilinear sample of a vector field on the atmosphere grid.
fn sample2(g: &AtmoGrid, f: &[(f64, f64)], lat: f64, lon: f64) -> (f64, f64) {
    let fy = ((90.0 - lat) / g.res - 0.5).clamp(0.0, (g.h - 1) as f64);
    let fx = (lon + 180.0) / g.res - 0.5;
    let (y0, ty) = (fy.floor() as usize, fy.fract());
    let y1 = (y0 + 1).min(g.h - 1);
    let x0f = fx.floor();
    let tx = fx - x0f;
    let x0 = (x0f as i64).rem_euclid(g.w as i64) as usize;
    let x1 = (x0 + 1) % g.w;
    let at = |y: usize, x: usize| f[y * g.w + x];
    let lerp = |a: (f64, f64), b: (f64, f64), t: f64| (a.0 + (b.0 - a.0) * t, a.1 + (b.1 - a.1) * t);
    let a = lerp(at(y0, x0), at(y0, x1), tx);
    let b = lerp(at(y1, x0), at(y1, x1), tx);
    lerp(a, b, ty)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atmosphere::Cells;

    /// One ocean basin between two meridional continents, Earth-like body: a clockwise
    /// subtropical gyre in the north with a fast poleward western boundary current, an
    /// anticlockwise one in the south, warm water on the western side and cold upwelling on
    /// the eastern side.
    #[test]
    fn gyre_has_western_intensification() {
        let g = AtmoGrid::new(1.0);
        let land: Vec<f64> = (0..g.len()).map(|i| if (i % g.w) < 40 || (i % g.w) > 120 { 1.0 } else { 0.0 }).collect();
        let a = Atmo { g, wet: true, land, elev: vec![0.0; g.len()] };
        let p = CirculationParams::default();
        let annual = season_winds(&a, &Season::zonal(Cells::for_body(24.0, 23.0, 0.0), 0.0, g.w), &p);
        let cur = currents(&a, &annual, 6371.0, &p);
        let at = |lat: f64, x: usize| cur[g.index(lat, g.lon(x))];
        // Western boundary (x just east of 40) at 30°N: strong northward flow.
        let wbc = at(31.0, 41);
        assert!(wbc.1 > 0.8, "western boundary current {wbc:?}");
        // Interior at 30°N: weak southward flow.
        let interior = at(31.0, 80);
        assert!(interior.1 < 0.0 && interior.1.abs() < wbc.1 / 10.0, "interior {interior:?}");
        // Southern hemisphere western boundary current runs south.
        assert!(at(-31.0, 41).1 < -0.8);
        let lats: Vec<f64> = (0..g.h).map(|y| g.lat(y)).collect();
        let earth = geoforge_core::Archetype::Earth.context(geoforge_core::Key(1), "Earth");
        let (tz, _) = Radiation::for_body(&earth, 3.2, 0.0).zonal_temperature(&lats);
        let sst = sst_anomaly(&a, &cur, [&annual, &annual], &tz, 6371.0, &p);
        let s = |lat: f64, x: usize| sst[g.index(lat, g.lon(x))];
        assert!(s(38.0, 43) > 1.0, "warm western boundary {}", s(38.0, 43));
        assert!(s(25.0, 119) < -0.5, "cool eastern boundary {}", s(25.0, 119));
    }
}
