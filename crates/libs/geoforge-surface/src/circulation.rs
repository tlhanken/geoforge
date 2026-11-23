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
//!   narrow western boundary current (Stommel), which leaves the coast short of where the
//!   gyres meet, climbs across the basin and runs up the far shore (the Gulf Stream and North
//!   Atlantic Drift), met at the coast by a cold current from the subpolar gyre (the
//!   Labrador Current); an equatorward jet along shores where the winds upwell (the
//!   California Current); an overturning branch drawing warm water poleward along the
//!   eastern side of a basin that opens wide onto a polar sea (the North Atlantic Current);
//!   plus equatorial currents and a circumpolar current where an ocean band has no shores.
//! - **SST anomaly**: currents carry water across the temperature gradient (poleward
//!   western boundary currents run warm, equatorward eastern ones cold), as far as a
//!   boundary current runs before it mixes, and that heat is carried downstream along the
//!   current; upwelling cools eastern margins under equatorward winds, in any season that
//!   has them, and the equator under the trades, most in a basin's east.
//!
//! Climate reads these products, so editing a wind or current map changes rainfall and
//! temperature downstream.

use crate::atmosphere::flow::{self, Keep};
use crate::atmosphere::grid::is_sea;
use crate::atmosphere::transport::semi_lagrangian;
use crate::atmosphere::{winds, Atmo, AtmoGrid, Radiation, Season};
use crate::climate::{Baseline, ClimateParams};
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
        Self { resolution_deg: 1.0, monsoon: 0.7, subtropical_highs: 0.5, boundary_current_m_s: 1.5, boundary_current_km: 150.0, drift_m_s: 0.3, heat_memory_days: 60.0, upwelling_c: 12.0 }
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
    feel_terrain(a, &season_winds_flat(a, season, p))
}

/// The season's winds as they would blow if the ground were flat: the cells, highs and
/// monsoons alone. Their rising and sinking is what the large-scale pattern sets.
#[must_use]
pub fn season_winds_flat(a: &Atmo, season: &Season, p: &CirculationParams) -> Vec<(f64, f64)> {
    let g = &a.g;
    let cells = &season.cells;
    let mut w = winds(a, season, p.monsoon);
    if p.subtropical_highs > 0.0 && season.locked.is_none() {
        // Anticyclones over each basin: poleward on the western side, equatorward on the east.
        let b = basins(a, 6371.0);
        let mut v_high = vec![0.0; g.len()];
        for (i, vh) in v_high.iter_mut().enumerate() {
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
            *vh = hemi * p.subtropical_highs * band * (std::f64::consts::PI * x_rel).cos();
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

/// Winds near the ground going round high ranges rather than over them. The stable air near
/// the ground is a layer up to ~3.5 km above sea level, thinning over ground higher than the
/// lowlands; the rising and sinking the cells and monsoons set up over flat ground stay.
#[must_use]
pub fn feel_terrain(a: &Atmo, w: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let g = &a.g;
    let depth: Vec<f64> = a.elev.iter().map(|&z| (AIR_LAYER_TOP_M - z).clamp(AIR_LAYER_MIN_SHARE * AIR_LAYER_M, AIR_LAYER_M)).collect();
    flow::conserve(g, 6371.0, &vec![true; g.len()], &depth, w, Keep::FlatLayer(AIR_LAYER_M))
}

/// Ocean surface currents (m/s, east and north) on the atmosphere grid.
/// Surface currents (m/s, east and north) and the share of each cell's flow that is the
/// overturning's warm poleward branch (0–1).
/// The ocean's own superrotation on a tidally locked world: an eastward current round the
/// equator, O(1 m/s), kept up by the stationary Rossby–Kelvin waves the day side's heating
/// raises even though the winds converge on the substellar point (Zeng & Yang 2021; the strong
/// eastward equatorial currents of Hu & Yang 2014). It is as wide as the equatorial
/// deformation radius, which grows as √(period): ±40° at a 30-day turn (an e-folding width
/// of 25° here). Slow rotators (60 days or more) keep it on the day side but lose most of it
/// on the night side.
#[derive(Debug, Clone, Copy)]
pub(crate) struct OceanJet {
    star: crate::atmosphere::Substellar,
    width_deg: f64,
    night_share: f64,
}

/// Speed of the ocean's superrotation, m/s.
const OCEAN_JET_M_S: f64 = 1.0;

impl OceanJet {
    pub(crate) fn for_body(body: &BodyContext) -> Option<Self> {
        let star = crate::atmosphere::Substellar::for_body(body)?;
        let period_d = body.exterior.rotation_h / 24.0;
        Some(Self { star, width_deg: (25.0 * (period_d / 30.0).sqrt()).clamp(8.0, 30.0), night_share: ((80.0 - period_d) / 30.0).clamp(0.2, 1.0) })
    }

    /// Eastward speed at a place, m/s.
    fn east_m_s(&self, lat: f64, lon: f64) -> f64 {
        let night = if self.star.zenith(lat, lon) > 90.0 { self.night_share } else { 1.0 };
        OCEAN_JET_M_S * night * (-(lat / self.width_deg).powi(2)).exp()
    }
}

fn currents(a: &Atmo, annual: &[(f64, f64)], radius_km: f64, p: &CirculationParams, jet: Option<&OceanJet>) -> (Vec<(f64, f64)>, Vec<f64>) {
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
            if n > 0.0 {
                s / n
            } else {
                0.0
            }
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
    let fronts = gyre_fronts(g, &tau);
    // Sverdrup interior closed by a Stommel western boundary layer.
    let psi: Vec<f64> = (0..g.len())
        .map(|i| {
            let y = i / g.w;
            let lat = g.lat(y);
            if a.land[i] >= 0.5 || b.open_row[y] || lat.abs() < 8.0 {
                return 0.0;
            }
            // Gyres span at most an ocean's width, and fade out next to a circumpolar band.
            let east = b.east_km[i].min(MAX_GYRE_KM);
            -c[y] * east * (1.0 - (-b.west_km[i] / p.boundary_current_km).exp()) * (1.0 - open_band[y])
        })
        .collect();
    let sea: Vec<bool> = a.land.iter().map(|&l| l < 0.5).collect();
    // The gyres' flow follows ψ's streamlines: fast and narrow along western shores, slow and
    // broad across the interior, and along a shore wherever one stands across its path, so
    // whatever flows into a stretch of sea flows out and none crosses a coast.
    let cur = flow::from_streamfunction(g, radius_km, &sea, &psi);
    let reach = overturning_reach(a, &b, fronts, radius_km);
    let mut overturning = vec![0.0; g.len()];
    // Scale so the strongest western boundary currents (meridional flow) reach the set speed.
    let mut speeds: Vec<f64> = cur.iter().map(|v| v.1.abs()).filter(|s| *s > 0.0).collect();
    speeds.sort_by(f64::total_cmp);
    let top = speeds.get(speeds.len().saturating_sub(1 + speeds.len() / 500)).copied().unwrap_or(0.0);
    let k = if top > 0.0 { p.boundary_current_m_s / top } else { 0.0 };
    let mut cur: Vec<(f64, f64)> = cur.iter().map(|v| (v.0 * k, v.1 * k)).collect();
    let mut drifts = vec![(0.0, 0.0); g.len()];
    let mut jets = vec![(0.0, 0.0); g.len()];
    for (i, v) in cur.iter_mut().enumerate() {
        let y = i / g.w;
        if a.land[i] >= 0.5 {
            continue;
        }
        let lat = g.lat(y);
        let hemi = lat.signum();
        let gyre = if hemi > 0.0 { fronts.0 } else { fronts.1 };
        // Poleward of where the western boundary current separates, the cold current coming
        // down from the subpolar gyre takes its place along the coast (the Labrador Current,
        // Oyashio and Malvinas Current), slower than the warm current it meets. (Where the
        // warm current stops, its water has turned offshore into the jet below.)
        if let Some(f) = gyre {
            if v.1 * hemi > 0.0 && lat.abs() < f.front + 1.0 {
                v.1 *= 1.0 - (1.0 + SUBPOLAR_COAST_SHARE) * smoothstep(f.separation, f.separation + 3.0, lat.abs());
            }
        }
        let drift = wind_drift(annual[i], lat);
        let d = &mut drifts[i];
        d.0 += drift.0;
        d.1 += drift.1;
        // Equatorial currents are faster than drift alone; a countercurrent just north of the
        // equator.
        let eq = (-(lat / 7.0).powi(2)).exp();
        let counter = (-((lat - 6.0) / 2.5).powi(2)).exp();
        d.0 += 0.3 * annual[i].0.min(0.0) * eq + 0.35 * counter;
        let j = &mut jets[i];
        // An equatorward jet along shores where the winds upwell: the eastern legs that close
        // the subtropical gyres.
        if !b.open_row[y] {
            j.1 -= lat.signum() * EASTERN_CURRENT_M_S * upwelling_drive([annual[i]; 2], lat, b.west_km[i] + b.east_km[i]) * (-b.east_km[i] / EASTERN_CURRENT_KM).exp();
        }
        // A circumpolar current where the ocean band has no shores (spread over the band).
        if open_band[y] > 0.0 {
            d.0 += 0.6 * open_band[y] * annual[i].0.max(0.0);
        }
        // Where the subtropical and subpolar gyres meet (under the strongest westerlies) the
        // western boundary current leaves the coast and crosses the basin, turning poleward
        // as it nears the far side and running on up the eastern shore (the Gulf Stream, North
        // Atlantic Drift and Norwegian Current; the Kuroshio, North Pacific and Alaska
        // Currents). Only in a basin of a gyre's width: across a wider stretch of open water
        // (the Atlantic and Indian Oceans together at 45°S) the current joins the circumpolar
        // flow instead.
        if let Some(GyreFront { front, separation }) = gyre {
            let basin = 1.0 - open_band[y];
            let width = b.west_km[i] + b.east_km[i];
            let crossing = 1.0 - smoothstep(MAX_GYRE_KM, 1.5 * MAX_GYRE_KM, width);
            let x_rel = if b.open_row[y] { 0.0 } else { b.west_km[i] / width.max(1.0) };
            // The jet leaves the coast at the separation latitude and climbs to the front over
            // the first part of the basin (Cape Hatteras to the Grand Banks), flowing along
            // its path.
            let path = separation + (front - separation) * smoothstep(0.0, JET_RISE_SHARE, x_rel);
            let rise_km = (front - separation) * smoothstep_slope(0.0, JET_RISE_SHARE, x_rel) / width.max(1.0) * radius_km.to_radians();
            let jet = (-((lat.abs() - path) / 5.0).powi(2)).exp();
            j.1 += hemi * basin * crossing * p.drift_m_s * jet * (1.0 - 0.5 * x_rel) * rise_km.min(1.0);
            let beyond = ((lat.abs() - front) / 5.0).clamp(0.0, 1.0) * ((front + 28.0 - lat.abs()) / 6.0).clamp(0.0, 1.0);
            // An eastern shore the current can run up: a coast that also reaches ~5° toward
            // the equator (Norway, Alaska, southern Chile), not Antarctica's, which runs
            // east–west.
            let (x, dx) = (i % g.w, g.res.to_radians() * radius_km * lat.to_radians().cos().max(0.02));
            let y_eq = (y as i64 + (hemi * 5.0 / g.res).round() as i64).clamp(0, g.h as i64 - 1) as usize;
            let runs_equatorward = !b.open_row[y] && a.land[y_eq * g.w + (x + (b.east_km[i] / dx).round() as usize + 1) % g.w] >= 0.5;
            let shore = if runs_equatorward { (-b.east_km[i] / 700.0).exp() } else { 0.0 };
            j.0 += basin * crossing * p.drift_m_s * jet * (1.0 - 0.5 * x_rel);
            j.1 += basin * p.drift_m_s * hemi * (0.6 * crossing * jet * x_rel + 0.8 * beyond * shore);
            // Where the subpolar sea opens wide onto a polar sea, water sinking there draws
            // warm water poleward along the basin's eastern side (the North Atlantic Current
            // into the Norwegian Sea).
            if reach[i] {
                let band = smoothstep(front - 5.0, front + 2.0, lat.abs()) * (1.0 - smoothstep(front + 18.0, front + 25.0, lat.abs()));
                overturning[i] = band * (-b.east_km[i] / OVERTURNING_EAST_KM).exp();
                j.1 += hemi * OVERTURNING_M_S * overturning[i];
            }
        }
    }
    // The jets are geostrophic currents like the gyres: none runs into a coast; what one
    // pushes against a shore is turned along it, and only the overturning's sinking stays. The
    // wind drift is the Ekman layer: water it piles against a coast sinks there and water it
    // draws off one wells up, so its convergence and divergence stay.
    if let Some(j) = jet {
        for (i, v) in jets.iter_mut().enumerate() {
            v.0 += j.east_m_s(g.lat(i / g.w), g.lon(i % g.w));
        }
    }
    let jets = flow::conserve(g, radius_km, &sea, &vec![1.0; g.len()], &jets, Keep::Interior);
    let cur: Vec<(f64, f64)> = (0..g.len()).map(|i| (cur[i].0 + drifts[i].0 + jets[i].0, cur[i].1 + drifts[i].1 + jets[i].1)).collect();
    (cur, overturning)
}

fn smoothstep(e0: f64, e1: f64, x: f64) -> f64 {
    let t = ((x - e0) / (e1 - e0)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

/// d[`smoothstep`]/dx.
fn smoothstep_slope(e0: f64, e1: f64, x: f64) -> f64 {
    let t = (x - e0) / (e1 - e0);
    if (0.0..=1.0).contains(&t) {
        6.0 * t * (1.0 - t) / (e1 - e0)
    } else {
        0.0
    }
}

/// Wind drift (m/s): ~3% of the wind speed, turned 30° to the right of the wind in the north
/// and to the left in the south (Ekman).
fn wind_drift(wind: (f64, f64), lat: f64) -> (f64, f64) {
    let (wu, wv) = (wind.0 * WIND_REF_M_S * 0.03, wind.1 * WIND_REF_M_S * 0.03);
    let turn = -lat.signum() * 30f64.to_radians();
    (wu * turn.cos() - wv * turn.sin(), wu * turn.sin() + wv * turn.cos())
}

/// Narrowest strait the overturning flows through, km: a sea cell nearer land on both sides
/// than this is shut (the Bering Strait is 85 km wide; the Denmark Strait 290).
const STRAIT_KM: f64 = 100.0;
/// Latitude beyond which open water is the polar sea, degrees.
const POLAR_SEA_LAT: f64 = 75.0;
/// Speed of the warm poleward flow the overturning draws, m/s.
const OVERTURNING_M_S: f64 = 0.25;
/// Warmth the overturning's deep warm water gives the surface it flows under, °C.
const OVERTURNING_HEAT_C: f64 = 4.0;
/// How far from the basin's eastern shore that flow reaches, km (e-folding).
const OVERTURNING_EAST_KM: f64 = 1500.0;

/// Subpolar sea that the overturning reaches: water from the gyre front poleward that joins a
/// polar sea through straits wider than [`STRAIT_KM`] without crossing an open circumpolar
/// band. Deep water forms where a salty subpolar sea opens onto a polar one (the Nordic Seas);
/// a sea closed off by a narrow strait (the North Pacific behind the Bering Strait) or a polar
/// ocean ringing a continent (the Southern Ocean) draws no warm water poleward.
fn overturning_reach(a: &Atmo, b: &Basins, fronts: (Option<GyreFront>, Option<GyreFront>), radius_km: f64) -> Vec<bool> {
    let g = &a.g;
    let sea = |x: usize, y: usize| a.land[y * g.w + x] < 0.5;
    let dy_km = g.res.to_radians() * radius_km;
    // Passable: sea with sea within STRAIT_KM on every side.
    let passable: Vec<bool> = (0..g.len())
        .map(|i| {
            let (x, y) = (i % g.w, i / g.w);
            let ky = (STRAIT_KM / dy_km).round().max(1.0) as i64;
            let kx = (STRAIT_KM / (dy_km * g.lat(y).to_radians().cos().max(0.05))).round().max(1.0) as i64;
            sea(x, y)
                && (-ky..=ky).all(|oy| {
                    let yy = (y as i64 + oy).clamp(0, g.h as i64 - 1) as usize;
                    (-kx..=kx).all(|ox| sea((x as i64 + ox).rem_euclid(g.w as i64) as usize, yy))
                })
        })
        .collect();
    let mut reach = vec![false; g.len()];
    for (north, front) in [(true, fronts.0), (false, fronts.1)] {
        let Some(front) = front else { continue };
        let inside = |y: usize| (g.lat(y) > 0.0) == north && g.lat(y).abs() >= front.separation && !b.open_row[y];
        let mut seen = vec![false; g.len()];
        let mut stack: Vec<usize> = (0..g.len()).filter(|&i| passable[i] && inside(i / g.w) && g.lat(i / g.w).abs() >= POLAR_SEA_LAT).collect();
        for &i in &stack {
            seen[i] = true;
        }
        let mut component = Vec::new();
        while let Some(i) = stack.pop() {
            component.push(i);
            let (x, y) = ((i % g.w) as i64, (i / g.w) as i64);
            for (ox, oy) in [(-1, 0), (1, 0), (0, -1), (0, 1)] {
                let ny = y + oy;
                if !(0..g.h as i64).contains(&ny) {
                    continue;
                }
                let j = ny as usize * g.w + (x + ox).rem_euclid(g.w as i64) as usize;
                if !seen[j] && passable[j] && inside(ny as usize) {
                    seen[j] = true;
                    stack.push(j);
                }
            }
        }
        // Only a polar sea that the gyre front's own water reaches.
        if component.iter().any(|&i| g.lat(i / g.w).abs() <= front.front) {
            for i in component {
                reach[i] = true;
            }
        }
    }
    // The coastal cells the strait test shut out still carry the flow along the shore.
    let mut out = reach.clone();
    for i in (0..g.len()).filter(|&i| reach[i]) {
        let (x, y) = ((i % g.w) as i64, (i / g.w) as i64);
        for oy in -1..=1 {
            for ox in -1..=1 {
                let ny = y + oy;
                if (0..g.h as i64).contains(&ny) {
                    let j = ny as usize * g.w + (x + ox).rem_euclid(g.w as i64) as usize;
                    if a.land[j] < 0.5 {
                        out[j] = true;
                    }
                }
            }
        }
    }
    out
}

/// Where the subtropical and subpolar gyres meet in a hemisphere, and where the western
/// boundary current leaves the coast, degrees of latitude.
#[derive(Debug, Clone, Copy)]
struct GyreFront {
    /// Under the strongest westerly stress, between 25° and 65°.
    front: f64,
    /// Where the western boundary current separates, short of the front: its momentum carries
    /// it off the coast before the gyres' boundary (the Gulf Stream at Cape Hatteras, 35°N,
    /// the Brazil Current at 38°S).
    separation: f64,
}

/// Each hemisphere's (north, south) gyre front.
fn gyre_fronts(g: &AtmoGrid, tau: &[f64]) -> (Option<GyreFront>, Option<GyreFront>) {
    let best = |north: bool| {
        let rows = || (0..g.h).filter(move |&y| (g.lat(y) > 0.0) == north);
        let front = rows().filter(|&y| (25.0..65.0).contains(&g.lat(y).abs()) && tau[y] > 0.0).max_by(|&a, &b| tau[a].total_cmp(&tau[b])).map(|y| g.lat(y).abs())?;
        // The westerlies' equatorward edge: the last row with easterly stress short of the front.
        let start = rows().filter(|&y| g.lat(y).abs() < front && tau[y] <= 0.0).map(|y| g.lat(y).abs()).fold(0.0, f64::max);
        Some(GyreFront { front, separation: start + SEPARATION_SHARE * (front - start) })
    };
    (best(true), best(false))
}

/// Depth of the air near the ground that high ranges block, m, and the least share of it
/// left over the highest ground.
const AIR_LAYER_M: f64 = 3000.0;
/// Height of the top of that layer above sea level, m: ground up to ~500 m does not thin it.
const AIR_LAYER_TOP_M: f64 = 3500.0;
const AIR_LAYER_MIN_SHARE: f64 = 0.1;
/// How far offshore coastal upwelling cools the sea, km (e-folding).
const UPWELLING_KM: f64 = 400.0;
/// Farthest a current carries water across the temperature gradient before it mixes with
/// the water beside it, km (soft limit).
const MAX_REACH_KM: f64 = 1500.0;
/// Equatorward wind along the shore, averaged over the year's upwelling seasons, at which
/// coastal upwelling and its jet are at full strength, in units of [`WIND_REF_M_S`] (3 m/s).
const UPWELLING_WIND: f64 = 0.375;
/// Speed of the equatorward jet along an upwelling shore (the California, Humboldt, Canary
/// and Benguela Currents), m/s.
const EASTERN_CURRENT_M_S: f64 = 0.2;
/// Width of that jet, km (e-folding from the shore).
const EASTERN_CURRENT_KM: f64 = 500.0;
/// Speed of the cold coastal current poleward of a western boundary current's separation,
/// as a share of the warm current's.
const SUBPOLAR_COAST_SHARE: f64 = 0.25;
/// How far from where the westerlies begin toward the gyre front a western boundary current
/// separates from the coast.
const SEPARATION_SHARE: f64 = 0.5;
/// Share of a basin's width over which the separated jet climbs from its separation latitude
/// to the gyre front.
const JET_RISE_SHARE: f64 = 0.5;
/// Half-width of the band where the trades upwell cold water along the equator, degrees.
const EQUATORIAL_UPWELLING_DEG: f64 = 4.0;
/// Equatorial upwelling's cooling as a share of a coast's.
const EQUATORIAL_UPWELLING_SHARE: f64 = 0.5;

/// How hard the seasons' winds drive upwelling along the eastern shore of a sea `width_km`
/// across (0–1): the equatorward wind along the shore, which drags the water offshore, from
/// ~6° to ~45° of latitude, off an ocean (a narrow sea has no deep cold water to bring up).
/// Each season upwells on its own and the opposite wind does not undo it, so a summer of
/// upwelling chills a coast all year (California, Iberia).
fn upwelling_drive(seasons: [(f64, f64); 2], lat: f64, width_km: f64) -> f64 {
    let equatorward = |w: (f64, f64)| (-w.1 * lat.signum()).max(0.0);
    let band = smoothstep(4.0, 8.0, lat.abs()) * (1.0 - smoothstep(40.0, 48.0, lat.abs()));
    let year = 0.5 * (equatorward(seasons[0]) + equatorward(seasons[1]));
    (year / UPWELLING_WIND).min(1.0) * band * smoothstep(1000.0, 3000.0, width_km)
}

/// How hard the trades drive upwelling along the equator (0–1): annual mean easterly wind,
/// near the equator, in an ocean wide enough for the trades to tilt its warm layer, growing
/// toward the basin's east where that layer is thin (the cold tongues of the eastern Pacific
/// and Atlantic; none among the Indonesian islands).
fn equatorial_drive(annual: (f64, f64), lat: f64, west_km: f64, east_km: f64) -> f64 {
    let easterly = (-annual.0 / UPWELLING_WIND).clamp(0.0, 1.0);
    let width = west_km + east_km;
    let x_rel = west_km / width.max(1.0);
    easterly * (-(lat / EQUATORIAL_UPWELLING_DEG).powi(2)).exp() * smoothstep(0.2, 0.9, x_rel) * smoothstep(2000.0, 5000.0, width)
}

/// Freezing point of seawater, °C.
const SEA_FREEZE_C: f64 = -1.8;
/// Latitude span over which the sea-surface anomaly's zonal mean is taken off, degrees.
const ZONAL_MEAN_SMOOTH_DEG: f64 = 10.0;

/// Sea-surface temperature anomaly (°C) from currents crossing the temperature gradient,
/// carried downstream, plus coastal upwelling. The gradient is the baseline's: north–south on a
/// world that turns, and on a locked world from the day side to the night side in every
/// direction, so its currents carry the day side's heat onto the night side where they cross
/// the terminator (the warm tongue east of the substellar point, Hu & Yang 2014).
fn sst_anomaly(a: &Atmo, cur: &[(f64, f64)], overturning: &[f64], seasons: [&[(f64, f64)]; 2], base: &Baseline, radius_km: f64, p: &CirculationParams) -> Vec<f64> {
    let g = &a.g;
    let b = basins(a, radius_km);
    let dy_km = g.res.to_radians() * radius_km;
    // The wind drift's part that is the same all round a band (the Ekman drift under the
    // westerlies, round a latitude circle) moves heat the band's temperature already carries:
    // counting it painted warm and cold bands right round the globe. Only its departures from
    // the band's mean, and the gyres, make anomalies.
    let annual = |i: usize| (0.5 * (seasons[0][i].0 + seasons[1][i].0), 0.5 * (seasons[0][i].1 + seasons[1][i].1));
    let mut drift_sum = vec![(0.0, 0.0, 0.0); g.h];
    for i in (0..g.len()).filter(|&i| a.land[i] < 0.5) {
        let (u, v) = wind_drift(annual(i), g.lat(i / g.w));
        let s = &mut drift_sum[base.band[i]];
        *s = (s.0 + u, s.1 + v, s.2 + 1.0);
    }
    let drift_mean: Vec<(f64, f64)> = drift_sum.iter().map(|&(u, v, n)| if n > 0.0 { (u / n, v / n) } else { (0.0, 0.0) }).collect();
    // What the currents carry across the temperature gradient, and what the winds bring up.
    let (carried, upwelled): (Vec<f64>, Vec<f64>) = (0..g.len())
        .map(|i| {
            if a.land[i] >= 0.5 {
                return (0.0, 0.0);
            }
            let (y, x) = (i / g.w, i % g.w);
            let (yn, ys) = (y.saturating_sub(1), (y + 1).min(g.h - 1));
            let (xw, xe) = ((x + g.w - 1) % g.w, (x + 1) % g.w);
            // The sea surface follows the baseline down only to freezing: over the ice at high
            // latitudes (or on the night side) its gradient is the air's, not the water's.
            let sst = |j: usize| base.t[j].max(SEA_FREEZE_C);
            let lat = g.lat(y);
            // °C per m northward and eastward.
            let dtdy = (sst(yn * g.w + x) - sst(ys * g.w + x)) / ((ys - yn).max(1) as f64 * dy_km * 1000.0);
            let dtdx = (sst(y * g.w + xe) - sst(y * g.w + xw)) / (2.0 * dy_km * lat.to_radians().cos().max(0.05) * 1000.0);
            // How far the water has come across the gradient: a current carries it for
            // `heat_memory_days`, but no further than a boundary current runs before eddies
            // mix it into the water beside it.
            let reach = |v: f64| {
                let m = v * p.heat_memory_days * 86_400.0;
                MAX_REACH_KM * 1000.0 * (m / (MAX_REACH_KM * 1000.0)).tanh()
            };
            let mean = drift_mean[base.band[i]];
            let across = reach(cur[i].0 - mean.0) * dtdx + reach(cur[i].1 - mean.1) * dtdy;
            // Upwelling: along an eastern shore where the winds blow equatorward, and along the
            // equator under the trades.
            let mut up = 0.0;
            if !b.open_row[y] {
                let (winds, width) = ([seasons[0][i], seasons[1][i]], b.west_km[i] + b.east_km[i]);
                let annual = (0.5 * (winds[0].0 + winds[1].0), 0.5 * (winds[0].1 + winds[1].1));
                up += upwelling_drive(winds, lat, width) * (-b.east_km[i] / UPWELLING_KM).exp();
                up += EQUATORIAL_UPWELLING_SHARE * equatorial_drive(annual, lat, b.west_km[i], b.east_km[i]);
            }
            // The overturning's branch is deep warm water that keeps giving up heat as it goes
            // poleward, past where the zonal gradient flattens at the freezing point.
            (-across + OVERTURNING_HEAT_C * overturning[i], -p.upwelling_c * up)
        })
        .unzip();
    // Carry the anomaly along the currents (one atmosphere cell per step at 0.5 m/s): water
    // keeps its warmth for months, so the Gulf Stream's reaches Norway.
    let flow: Vec<(f64, f64)> = cur.iter().map(|v| (v.0 / 0.5, v.1 / 0.5)).collect();
    let source: Vec<f64> = carried.iter().zip(&upwelled).map(|(c, u)| c + u).collect();
    let mut anom = carry_along(g, &flow, &source, &a.land, 200, 0.97);
    diffuse_sea(g, &mut anom, &a.land, 0.3);
    // A departure from the baseline averages to zero round each band. The mean taken off is
    // smoothed over ~10° of latitude so that it cannot itself draw bands where one row's mean
    // differs from the next.
    let mut rows = vec![(0.0, 0.0); g.h];
    for i in (0..g.len()).filter(|&i| a.land[i] < 0.5) {
        let r = &mut rows[base.band[i]];
        *r = (r.0 + anom[i], r.1 + 1.0);
    }
    let half = (ZONAL_MEAN_SMOOTH_DEG / 2.0 / g.res).round() as i64;
    let mean: Vec<f64> = (0..g.h as i64)
        .map(|y| {
            let (s, n) = (y - half..=y + half).filter(|k| (0..g.h as i64).contains(k)).fold((0.0, 0.0), |(s, n), k| {
                let w = 1.0 - (k - y).abs() as f64 / (half + 1) as f64;
                (s + w * rows[k as usize].0, n + w * rows[k as usize].1)
            });
            if n > 0.0 {
                s / n
            } else {
                0.0
            }
        })
        .collect();
    for (i, v) in anom.iter_mut().enumerate() {
        if a.land[i] < 0.5 {
            *v -= mean[base.band[i]];
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
    const VERSION: u32 = 12;
    const LAYER: LayerId = LayerId::Circulation;

    fn needs() -> &'static [ProductId] {
        &[ProductId::RoughElevation]
    }
    fn applies(body: &BodyContext) -> Result<(), Precondition> {
        crate::preconditions::surface_air(body)
    }

    fn run(ctx: &RunContext<'_>, inputs: &Products, p: &CirculationParams) -> Result<Products, Error> {
        let mut out = Products::new();
        if [ProductId::WindJan, ProductId::WindJul, ProductId::OceanCurrents, ProductId::SstAnomaly].iter().all(|x| inputs.contains(*x)) {
            return Ok(out);
        }
        let elev = &inputs.rough_elevation().ok_or(Error::MissingInput(ProductId::RoughElevation))?.values;
        let grid = ctx.grid;
        let body = ctx.body;
        let r = ctx.radius_km();
        let wet = body.interior.has_surface_liquid_water();
        let atmo = Atmo::from_world(grid, elev, p.resolution_deg, wet);
        let g: AtmoGrid = atmo.g;
        let jul = season_winds(&atmo, &Season::for_body(&atmo, body, 1.0, p.monsoon), p);
        let jan = season_winds(&atmo, &Season::for_body(&atmo, body, -1.0, p.monsoon), p);
        let annual: Vec<(f64, f64)> = jul.iter().zip(&jan).map(|(a, b)| (0.5 * (a.0 + b.0), 0.5 * (a.1 + b.1))).collect();
        let has_sea = atmo.land.iter().any(|l| *l < 0.5);
        let jet = OceanJet::for_body(body);
        let (cur, overturning) = if has_sea { currents(&atmo, &annual, r, p, jet.as_ref()) } else { (vec![(0.0, 0.0); g.len()], vec![0.0; g.len()]) };
        let lats: Vec<f64> = (0..g.h).map(|y| g.lat(y)).collect();
        let cp = ClimateParams::default();
        let radiation = Radiation::for_body(body, cp.heat_transport, cp.greenhouse_offset_c);
        let locked = crate::atmosphere::Substellar::for_body(body);
        let (tz, _) = if locked.is_some() { radiation.substellar_temperature(&lats) } else { radiation.zonal_temperature(&lats) };
        let base = Baseline::new(&g, &lats, &tz, locked.as_ref());
        let sst = if has_sea { sst_anomaly(&atmo, &cur, &overturning, [&jan, &jul], &base, r, p) } else { vec![0.0; g.len()] };

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

    /// East and north components on the atmosphere grid.
    type Vectors = Vec<(f64, f64)>;

    /// One ocean basin between two meridional continents (x 40–120) on an Earth-like body,
    /// with its annual winds and currents.
    fn one_basin() -> (Atmo, CirculationParams, Vectors, (Vectors, Vec<f64>)) {
        basin_with(|_, _| false)
    }

    /// The same basin with extra land where `land(x, y)` says.
    fn basin_with(extra: impl Fn(usize, usize) -> bool) -> (Atmo, CirculationParams, Vectors, (Vectors, Vec<f64>)) {
        let g = AtmoGrid::new(1.0);
        let land: Vec<f64> = (0..g.len()).map(|i| if (i % g.w) < 40 || (i % g.w) > 120 || extra(i % g.w, i / g.w) { 1.0 } else { 0.0 }).collect();
        let a = Atmo { g, wet: true, land, elev: vec![0.0; g.len()] };
        let p = CirculationParams::default();
        // The year's mean of July's and January's winds, as the generator takes it.
        let [jul, jan] = [1.0, -1.0].map(|sign| season_winds(&a, &Season::new(&a, 24.0, 6371.0, 23.0, sign, p.monsoon), &p));
        let annual: Vectors = jul.iter().zip(&jan).map(|(n, s)| (0.5 * (n.0 + s.0), 0.5 * (n.1 + s.1))).collect();
        let ocean = currents(&a, &annual, 6371.0, &p, None);
        (a, p, annual, ocean)
    }

    fn sst_of(a: &Atmo, p: &CirculationParams, annual: &[(f64, f64)], (cur, overturning): &(Vectors, Vec<f64>)) -> Vec<f64> {
        let lats: Vec<f64> = (0..a.g.h).map(|y| a.g.lat(y)).collect();
        let earth = geoforge_core::Archetype::Earth.context(geoforge_core::Key(1), "Earth");
        let (tz, _) = Radiation::for_body(&earth, 3.2, 0.0).zonal_temperature(&lats);
        sst_anomaly(a, cur, overturning, [annual, annual], &Baseline::new(&a.g, &lats, &tz, None), 6371.0, p)
    }

    /// A locked ocean world turning once in ten days: its currents carry the day side's warmth
    /// across the terminators onto the night side (Hu & Yang 2014), and what they carry
    /// averages to nothing round each ring about the point under the star.
    #[test]
    fn a_locked_ocean_carries_heat_to_the_night_side() {
        let mut body = geoforge_core::Archetype::Earth.context(geoforge_core::Key(1), "Earth");
        (body.irradiance.year_days, body.exterior.rotation_h) = (10.0, 240.0);
        let g = AtmoGrid::new(2.0);
        let a = Atmo { g, wet: true, land: vec![0.0; g.len()], elev: vec![0.0; g.len()] };
        let p = CirculationParams::default();
        let [jul, jan] = [1.0, -1.0].map(|sign| season_winds(&a, &Season::for_body(&a, &body, sign, p.monsoon), &p));
        let annual: Vectors = jul.iter().zip(&jan).map(|(n, s)| (0.5 * (n.0 + s.0), 0.5 * (n.1 + s.1))).collect();
        let (cur, overturning) = currents(&a, &annual, 6371.0, &p, OceanJet::for_body(&body).as_ref());
        let star = crate::atmosphere::Substellar::for_body(&body).expect("locked");
        let lats: Vec<f64> = (0..g.h).map(|y| g.lat(y)).collect();
        let (tz, _) = Radiation::for_body(&body, 3.2, 0.0).substellar_temperature(&lats);
        let base = Baseline::new(&g, &lats, &tz, Some(&star));
        let sst = sst_anomaly(&a, &cur, &overturning, [&annual, &annual], &base, 6371.0, &p);
        let night_warmth = (0..g.len()).filter(|&i| star.zenith(g.lat(i / g.w), g.lon(i % g.w)) > 95.0).map(|i| sst[i]).fold(f64::MIN, f64::max);
        let spread = sst.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        let (east, west) = (sst[g.index(0.0, 100.0)], sst[g.index(0.0, -100.0)]);
        assert!(night_warmth > 1.0, "warm water on the night side: {night_warmth} °C");
        // The jet runs east: past the evening terminator the water is warm, past the morning
        // one it is not.
        assert!(east > 1.0 && east > west + 1.0, "east of the terminator {east}, west {west}");
        let mut rings = vec![(0.0, 0.0); g.h];
        for i in 0..g.len() {
            let r = &mut rings[base.band[i]];
            *r = (r.0 + sst[i], r.1 + 1.0);
        }
        assert!(rings.iter().filter(|r| r.1 > 50.0).all(|r| (r.0 / r.1).abs() < 0.5 * spread), "rings {rings:?}");
    }

    /// A clockwise subtropical gyre in the north with a fast poleward western boundary
    /// current, an anticlockwise one in the south, warm water on the western side and cold
    /// upwelling on the eastern side.
    #[test]
    fn gyre_has_western_intensification() {
        let (a, p, annual, ocean) = one_basin();
        let g = a.g;
        let at = |lat: f64, x: usize| ocean.0[g.index(lat, g.lon(x))];
        // Western boundary (x just east of 40) at 30°N: strong northward flow.
        let wbc = at(31.0, 41);
        assert!(wbc.1 > 0.8, "western boundary current {wbc:?}");
        // Interior at 30°N: weak southward flow.
        let interior = at(31.0, 80);
        assert!(interior.1 < 0.0 && interior.1.abs() < wbc.1 / 10.0, "interior {interior:?}");
        // Southern hemisphere western boundary current runs south.
        assert!(at(-31.0, 41).1 < -0.8);
        let sst = sst_of(&a, &p, &annual, &ocean);
        let s = |lat: f64, x: usize| sst[g.index(lat, g.lon(x))];
        // Warm along the western shore below where the current leaves it.
        assert!(s(33.0, 43) > 1.0 && s(33.0, 43) > s(33.0, 80) + 1.0, "warm western boundary {}", s(33.0, 43));
        assert!(s(25.0, 119) < -0.5, "cool eastern boundary {}", s(25.0, 119));
    }

    /// The gyres close: an equatorward current along the eastern shore in both hemispheres,
    /// and, poleward of where the western boundary current leaves the coast, a cold current
    /// coming down the western shore to meet it. The trades upwell cold water along the
    /// equator, most in the east of the basin.
    #[test]
    fn gyres_close_and_the_equator_upwells() {
        let (a, p, annual, ocean) = one_basin();
        let g = a.g;
        let at = |lat: f64, x: usize| ocean.0[g.index(lat, g.lon(x))];
        for lat in [25.0, 32.0] {
            assert!(at(lat, 118).1 < -0.1, "eastern leg at {lat}°N {:?}", at(lat, 118));
            assert!(at(-lat, 118).1 > 0.1, "eastern leg at {lat}°S {:?}", at(-lat, 118));
        }
        let (north, _) = gyre_fronts(&g, &(0..g.h).map(|y| annual[y * g.w + 80].0 * annual[y * g.w + 80].0.abs()).collect::<Vec<_>>());
        let f = north.expect("a gyre front");
        assert!(f.separation > 30.0 && f.separation < f.front - 5.0, "{f:?}");
        let coast = at(f.separation + 4.0, 41);
        assert!(coast.1 < 0.0, "cold coastal current poleward of separation {coast:?}");
        let sst = sst_of(&a, &p, &annual, &ocean);
        let s = |lat: f64, x: usize| sst[g.index(lat, g.lon(x))];
        assert!(s(f.separation + 4.0, 41) < s(f.separation - 4.0, 41) - 2.0, "warm current meets cold at the separation");
        assert!(s(0.5, 112) < s(0.5, 50) - 1.0, "equatorial cold tongue in the east: {} vs {}", s(0.5, 112), s(0.5, 50));
    }

    /// A basin that opens wide onto its polar sea draws warm water poleward along its eastern
    /// side; the same basin closed off by a one-cell strait at 66°N does not.
    #[test]
    fn overturning_needs_an_open_polar_sea() {
        let (a, p, annual, open) = one_basin();
        let strait = |x: usize, y: usize| (66..=68).contains(&(90 - y as i64)) && x != 80;
        let (b, _, annual_b, closed) = basin_with(strait);
        let g = a.g;
        let east = g.index(58.0, g.lon(115));
        assert!(open.1[east] > 0.3, "overturning off the open basin's eastern shore {}", open.1[east]);
        assert!(closed.1[east] == 0.0, "none behind a strait {}", closed.1[east]);
        let (so, sc) = (sst_of(&a, &p, &annual, &open), sst_of(&b, &p, &annual_b, &closed));
        // One basin is the whole latitude circle here, so most of the warmth goes into the
        // zonal mean the anomaly is taken from; the eastern side keeps its lead.
        assert!(so[east] > sc[east] + 0.3, "warmer for it: {} vs {}", so[east], sc[east]);
    }
}
