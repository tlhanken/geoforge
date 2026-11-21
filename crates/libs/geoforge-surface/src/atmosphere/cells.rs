//! The general circulation by rotation: Hadley, Ferrel and polar cells, and monsoons.
//!
//! [`Cells`] is the zonal pattern of one season. [`Season`] adds the longitudes: the ITCZ
//! follows the sun over the oceans but is pulled far into the summer hemisphere over large
//! heated land masses (India in July, southern Africa and Australia in January). Between
//! the equator and such an ITCZ the trades cross the equator and turn into monsoon
//! westerlies.

use super::grid::Atmo;
use std::f64::consts::PI;

/// How far the ITCZ moves over land, as a share of the obliquity per unit of summer-side
/// land fraction (on top of the zonal 0.3 × obliquity). Fitted to Earth's July monsoon
/// trough over India and Africa (~18–20°N) and January's over Australia and Africa.
const LAND_ITCZ_SHIFT: f64 = 0.5;
/// Rain falls on the equatorward side of an ITCZ pulled over land: the West African rain
/// belt peaks ~8° south of the surface convergence, the Amazon's stays near the equator
/// while January's convergence reaches the Chaco.
const RAIN_LAG: f64 = 0.7;

/// Extra pull of a summer heat low per km of ground height.
const ELEVATED_HEATING: f64 = 0.4;

/// Hadley, Ferrel and polar cells of one season.
#[derive(Debug, Clone, Copy)]
pub struct Cells {
    /// Latitude of the ITCZ over open ocean, degrees.
    pub itcz: f64,
    /// Width of each Hadley cell, degrees from the ITCZ.
    pub hadley: f64,
    /// Poleward edge of the Ferrel cells, degrees.
    pub ferrel: f64,
}

/// 0 below `lo`, 1 inside, 0 above `hi`, with smooth shoulders `soft` wide.
fn plateau(x: f64, lo: f64, hi: f64, soft: f64) -> f64 {
    let s = |t: f64| {
        let t = t.clamp(0.0, 1.0);
        t * t * (3.0 - 2.0 * t)
    };
    s((x - lo) / soft + 0.5) * (1.0 - s((x - hi) / soft + 0.5))
}

impl Cells {
    #[must_use]
    pub fn for_body(rotation_h: f64, obliquity: f64, season: f64) -> Self {
        // Hadley width ~ Ω^(-2/3) (GCM scaling); Earth: 30°.
        let hadley = (30.0 * (rotation_h.abs().max(1.0) / 24.0).powf(2.0 / 3.0)).clamp(12.0, 80.0);
        let ferrel = (hadley * 2.0).min(85.0);
        Self { itcz: season * 0.3 * obliquity.min(45.0), hadley, ferrel }
    }

    /// Subtropical edge of the Hadley cell on the side of `lat`.
    fn edge(&self, lat: f64) -> f64 {
        if lat >= self.itcz { self.itcz + self.hadley } else { self.itcz - self.hadley }
    }

    /// Poleward edge of the Ferrel cell (the polar front) on the side of `lat`: it follows
    /// the sun like the rest of the pattern.
    fn polar_front(&self, lat: f64) -> f64 {
        (self.ferrel + 4.0 + 0.5 * self.itcz * lat.signum()).clamp((self.hadley + 10.0).min(88.0), 88.0)
    }

    /// Surface wind (east, north), in units of a typical trade wind, of the zonal pattern.
    #[must_use]
    pub fn wind(&self, lat: f64) -> (f64, f64) {
        self.wind_at(lat, self.itcz)
    }

    /// Surface wind where the column's ITCZ lies at `itcz`.
    #[must_use]
    pub fn wind_at(&self, lat: f64, itcz: f64) -> (f64, f64) {
        let (north, south) = (self.itcz + self.hadley, self.itcz - self.hadley);
        let itcz = itcz.clamp(south + 5.0, north - 5.0);
        let hemi = if lat >= 0.0 { 1.0 } else { -1.0 };
        if lat > south && lat < north {
            let side = if lat >= itcz { 1.0 } else { -1.0 };
            let edge = if side > 0.0 { north } else { south };
            let x = ((lat - itcz) / (edge - itcz)).clamp(0.0, 1.0);
            let s = (PI * x).sin().max(0.0);
            // Trades blow toward the ITCZ, turned west where they flow toward the equator; once
            // across it, toward an ITCZ in the summer hemisphere, they turn east over ~10°.
            let poleward = lat * side < 0.0;
            let u = if poleward {
                let q = (lat.abs() / 10.0).min(1.0);
                s * (0.7 * q - 0.9 * (1.0 - q))
            } else {
                -0.9 * s * (lat.abs() / 5.0).min(1.0)
            };
            return (u, -0.45 * s * side);
        }
        let edge = self.edge(lat);
        let front = self.polar_front(lat);
        let ferrel_w = (front - edge.abs()).max(5.0);
        let x2 = (lat.abs() - edge.abs()) / ferrel_w;
        if x2 < 1.0 {
            let s = (PI * x2.max(0.0)).sin();
            return (1.0 * s, 0.3 * s * hemi);
        }
        let x3 = ((lat.abs() - front) / (90.0 - front).max(1.0)).clamp(0.0, 1.0);
        let s = (PI * x3).sin();
        (-0.6 * s, -0.35 * s * hemi)
    }

    /// Rain-out multiplier from large-scale vertical motion over open sea, zonal pattern.
    #[must_use]
    pub fn convergence(&self, lat: f64) -> f64 {
        self.convergence_at(lat, self.itcz, 0.0)
    }

    /// Rain-out multiplier from large-scale vertical motion where the column's ITCZ lies at
    /// `itcz`, over ground with land fraction `land` (heated land spreads the tropical rain
    /// belt wider than over the sea).
    #[must_use]
    pub fn convergence_at(&self, lat: f64, itcz: f64, land: f64) -> f64 {
        let rain_belt = itcz - RAIN_LAG * (itcz - self.itcz);
        let rel = (lat - rain_belt).abs();
        let ascent = 0.9 * (-(rel / (7.0 + 3.0 * land)).powi(2)).exp();
        // The sinking branch of each Hadley cell: a broad belt inside its subtropical edge,
        // strongest in the winter hemisphere (whose cell reaches across the equator). Each
        // cell's belt is judged on its own side, so the rain-out is continuous across the ITCZ.
        let belt = |poleward: f64, edge: f64, winter: f64| {
            (0.525 + 0.225 * winter.clamp(-1.0, 1.0)) * plateau(poleward, edge - 0.55 * self.hadley, edge + 0.1 * self.hadley, 8.0)
        };
        let north = belt(lat, self.itcz + self.hadley, -self.itcz / 7.0);
        let south = belt(-lat, self.hadley - self.itcz, self.itcz / 7.0);
        let sink = north.max(south);
        (1.0 + ascent - sink + self.storminess(lat)).max(0.2)
    }

    /// Mid-latitude storm tracks, and frequent frontal weather poleward of them (0–1.3). The
    /// storm tracks follow the sun: poleward in summer, equatorward in winter, bringing winter
    /// rain to the subtropical west coasts.
    #[must_use]
    pub fn storminess(&self, lat: f64) -> f64 {
        let storm_lat = ((self.hadley + 18.0).min(self.ferrel) + 0.5 * self.itcz * lat.signum()).max(self.hadley);
        let storm = 0.8 * (-((lat.abs() - storm_lat) / 12.0).powi(2)).exp();
        // Fronts come often poleward of the track, fading in over ~12° (a step here drew
        // rain lines along the latitude circles).
        let subpolar = 0.5 * plateau(lat.abs(), storm_lat, 200.0, 12.0);
        storm + subpolar
    }
}

/// One solstice season: the zonal cells and the ITCZ in each column of the atmosphere grid.
#[derive(Debug, Clone)]
pub struct Season {
    pub cells: Cells,
    /// +1 northern summer (July), −1 northern winter (January), 0 equinox.
    pub sign: f64,
    /// ITCZ latitude in each grid column, degrees.
    pub itcz: Vec<f64>,
}

impl Season {
    /// The ITCZ sits at the zonal latitude over open ocean and moves into the summer
    /// hemisphere over land, in proportion to the land (high land counts more: it heats the
    /// air higher up) on the summer side of the tropics. `monsoon` scales the shift
    /// (0 = zonal everywhere; 0.7 = Earth).
    #[must_use]
    pub fn new(a: &Atmo, rotation_h: f64, obliquity: f64, sign: f64, monsoon: f64) -> Self {
        let cells = Cells::for_body(rotation_h, obliquity, sign);
        let g = &a.g;
        let (lo, hi) = (8.0 * cells.hadley / 30.0, 35.0 * cells.hadley / 30.0);
        let heat: Vec<f64> = (0..g.w)
            .map(|x| {
                let (mut s, mut w) = (0.0, 0.0);
                for y in 0..g.h {
                    let lat = g.lat(y) * sign;
                    if sign != 0.0 && lat >= lo && lat <= hi {
                        let i = y * g.w + x;
                        let c = g.lat(y).to_radians().cos();
                        s += c * a.land[i] * (1.0 + 0.25 * a.elev[i].max(0.0) / 1000.0).min(2.0);
                        w += c;
                    }
                }
                if w > 0.0 { (s / w).min(1.2) } else { 0.0 }
            })
            .collect();
        // Smooth a little along the latitude circle; wider, and the trough would cross onto
        // the oceans beside a continent (off Peru, the trades blow all year).
        let r = (5.0 / g.res).round().max(1.0) as i64;
        let shift = LAND_ITCZ_SHIFT * obliquity.min(45.0) * (monsoon / 0.7).clamp(0.0, 2.0);
        let itcz = (0..g.w as i64)
            .map(|x| {
                let l: f64 = (-r..=r).map(|k| heat[(x + k).rem_euclid(g.w as i64) as usize]).sum::<f64>() / (2 * r + 1) as f64;
                cells.itcz + sign * shift * l
            })
            .collect();
        Self { cells, sign, itcz }
    }

    /// The zonal pattern alone (no land): for winds without a surface, and tests.
    #[must_use]
    pub fn zonal(cells: Cells, sign: f64, columns: usize) -> Self {
        Self { cells, sign, itcz: vec![cells.itcz; columns] }
    }

    #[must_use]
    pub fn wind(&self, x: usize, lat: f64) -> (f64, f64) {
        self.cells.wind_at(lat, self.itcz[x])
    }

    #[must_use]
    pub fn convergence(&self, x: usize, lat: f64, land: f64) -> f64 {
        self.cells.convergence_at(lat, self.itcz[x], land)
    }
}

/// How much of a low-level flow starting at (`lat`, `lon`) in direction (`du`, `dv`) gets
/// past the terrain ahead (~600 km): 1 over open ground, 0 against a range ≥ 2.3 km higher.
fn blocked_by_terrain(a: &Atmo, lat: f64, lon: f64, du: f64, dv: f64) -> f64 {
    let g = &a.g;
    let here = g.sample(&a.elev, lat, lon);
    let c = lat.to_radians().cos().max(0.15);
    let rise = (1..=6)
        .map(|k| {
            let d = k as f64 * g.res.max(1.0);
            g.sample(&a.elev, (lat + dv * d).clamp(-89.9, 89.9), lon + du * d / c) - here
        })
        .fold(0.0, f64::max);
    (1.0 - (rise - 800.0) / 1500.0).clamp(0.0, 1.0)
}

/// Winds for one season: the cells, plus monsoon flow into the summer heat lows over land
/// (strongest over the subtropics, which heat most) and out of the winter highs, turned by
/// the Coriolis force (cyclonically into a low, anticyclonically out of a high).
pub fn winds(a: &Atmo, season: &Season, monsoon: f64) -> Vec<(f64, f64)> {
    let g = &a.g;
    let cells: Vec<(f64, f64)> = (0..g.len()).map(|i| season.wind(i % g.w, g.lat(i / g.w))).collect();
    let m = monsoon_flow(a, season, monsoon);
    cells.iter().zip(&m).map(|(c, m)| (c.0 + m.0, c.1 + m.1)).collect()
}

/// The monsoon part of [`winds`].
pub fn monsoon_flow(a: &Atmo, season: &Season, monsoon: f64) -> Vec<(f64, f64)> {
    let g = &a.g;
    if monsoon <= 0.0 {
        return vec![(0.0, 0.0); g.len()];
    }
    let rad = (12.0 / g.res).round().max(1.0) as usize;
    // High ground heats the air higher up, so summer lows sit over plateaus and high deserts
    // (Tibet, the Mexican and Colorado plateaus, Iran), not over the middle of a continent.
    let heated: Vec<f64> = (0..g.len())
        .map(|i| {
            let lat = g.lat(i / g.w) * season.sign;
            let lift = 1.0 + ELEVATED_HEATING * (a.elev[i].max(0.0) / 1000.0).min(4.0);
            a.land[i] * (0.35 + 0.65 * (-((lat - 25.0) / 15.0).powi(2)).exp()) * lift
        })
        .collect();
    let summer_low = g.blur(&heated, rad);
    let winter_high = g.blur(&a.land, rad);
    (0..g.len())
        .map(|i| {
            let (x, y) = (i % g.w, i / g.w);
            let lat = g.lat(y);
            let summer = lat * season.sign > 0.0;
            let ls = if summer { &summer_low } else { &winter_high };
            let c = lat.to_radians().cos().max(0.1);
            let dx = (ls[y * g.w + (x + 1) % g.w] - ls[y * g.w + (x + g.w - 1) % g.w]) / (2.0 * g.res * c);
            let (yn, ys) = (y.saturating_sub(1), (y + 1).min(g.h - 1));
            let dy = (ls[yn * g.w + x] - ls[ys * g.w + x]) / ((ys - yn).max(1) as f64 * g.res);
            let mag = (dx * dx + dy * dy).sqrt();
            if mag <= 1e-9 {
                return (0.0, 0.0);
            }
            let taper = 1.0 - ((lat.abs() - 25.0) / 20.0).clamp(0.0, 1.0);
            let strength = monsoon * taper * (mag / 0.02).min(1.0) * if summer { 1.0 } else { -0.4 };
            let (u, v) = (strength * dx / mag, strength * dy / mag);
            // High ranges stand in the way: the Andes keep Amazonian heat from pulling the
            // Pacific ashore at Peru.
            let ahead = blocked_by_terrain(a, lat, g.lon(x), u / strength.abs(), v / strength.abs());
            let (u, v) = (u * ahead, v * ahead);
            // Turned right of the pressure gradient in the north, left in the south, little
            // near the equator.
            let turn = -50f64.to_radians() * (lat / 15.0).clamp(-1.0, 1.0);
            (u * turn.cos() - v * turn.sin(), u * turn.sin() + v * turn.cos())
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atmosphere::AtmoGrid;

    #[test]
    fn hadley_cell_scales_with_rotation() {
        let earth = Cells::for_body(24.0, 23.0, 0.0);
        let slow = Cells::for_body(240.0, 23.0, 0.0);
        assert!((earth.hadley - 30.0).abs() < 1e-9);
        assert!(slow.hadley > 60.0);
        // Trades in the tropics, westerlies in mid-latitudes.
        assert!(earth.wind(15.0).0 < -0.5 && earth.wind(-15.0).0 < -0.5);
        assert!(earth.wind(45.0).0 > 0.5 && earth.wind(-45.0).0 > 0.5);
        // Wet at the ITCZ, dry under the subtropical highs.
        assert!(earth.convergence(0.0) > 1.5 && earth.convergence(22.0) < 0.6 && earth.convergence(-22.0) < 0.6);
    }

    /// A continent filling the northern subtropics of half the globe pulls July's ITCZ north
    /// over it, and the trades that cross the equator to reach it blow from the west.
    #[test]
    fn monsoon_itcz_follows_land() {
        let g = AtmoGrid::new(2.0);
        let land: Vec<f64> = (0..g.len()).map(|i| if i % g.w < g.w / 2 && g.lat(i / g.w) > 5.0 && g.lat(i / g.w) < 40.0 { 1.0 } else { 0.0 }).collect();
        let a = Atmo { g, wet: true, land, elev: vec![0.0; g.len()] };
        let jul = Season::new(&a, 24.0, 23.4, 1.0, 0.7);
        let (over_land, over_sea) = (jul.itcz[g.w / 4], jul.itcz[3 * g.w / 4]);
        assert!(over_land > 15.0 && (over_sea - 7.0).abs() < 0.5, "land {over_land} sea {over_sea}");
        let (u, v) = jul.wind(g.w / 4, 8.0);
        assert!(u > 0.2 && v > 0.0, "monsoon westerlies {u} {v}");
        let jan = Season::new(&a, 24.0, 23.4, -1.0, 0.7);
        assert!((jan.itcz[g.w / 4] + 7.0).abs() < 0.5, "no southern land, zonal ITCZ in January");
    }
}
