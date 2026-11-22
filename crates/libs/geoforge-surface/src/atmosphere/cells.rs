//! The general circulation by rotation and size: Hadley cells, eddy-driven jets with the
//! cells between them, and monsoons.
//!
//! The Hadley cell narrows as rotation quickens (width ∝ Ω^(−2/3); Earth 30°). Poleward of
//! it, eddy-driven jets are spaced by the Rhines scale (Chemke & Kaspi 2015), so a world
//! that spins faster or is larger has more of them: at four times Earth's rotation, or with
//! a larger radius, idealized GCMs show about three per hemisphere (GFDL). Each jet has
//! surface westerlies under it and easterlies on its poleward side, so a hemisphere has
//! 1 + 2·jets cells: Earth one jet and three cells (Hadley, Ferrel, polar), a fast or big
//! world five or seven, and a slow one like Venus or Titan no jet and a Hadley cell reaching
//! to a polar cap.
//!
//! A world that turns once per orbit keeps one face to its star ([`Substellar`]): the air
//! rises over the substellar point and flows back toward it along the ground from every
//! side, the night side under sinking air; faster-turning ones add an equatorial
//! superrotating jet (Haqq-Misra et al. 2018: slow rotators beyond ~20-day periods, Rhines
//! rotators below).
//!
//! [`Cells`] is the zonal pattern of one season. [`Season`] adds the longitudes: the ITCZ
//! follows the sun over the oceans but is pulled far into the summer hemisphere over large
//! heated land masses (India in July, southern Africa and Australia in January). Between
//! the equator and such an ITCZ the trades cross the equator and turn into monsoon
//! westerlies.

use super::grid::Atmo;
use geoforge_core::BodyContext;
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
    /// Eddy-driven jets per hemisphere ([`jets_for`]).
    pub jets: u32,
}

/// Eddy-driven jets per hemisphere: jet spacing follows the Rhines scale, ∝ (a·Ω)^(−1/2) at
/// fixed wind speed, and the winds themselves weaken with faster rotation; fitted as
/// round((a·Ω / a⊕·Ω⊕)^0.75) through Earth's one jet and the three of a world spinning four
/// times faster. Mars: 1; a 1.8 R⊕ world with Earth's day: 2; Venus and Titan: 0.
#[must_use]
pub fn jets_for(rotation_h: f64, radius_km: f64) -> u32 {
    let a_omega = (radius_km / 6371.0) * (24.0 / rotation_h.abs().max(1.0));
    a_omega.powf(0.75).round().min(12.0) as u32
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
    pub fn for_body(rotation_h: f64, radius_km: f64, obliquity: f64, season: f64) -> Self {
        // Hadley width ~ Ω^(-2/3) (GCM scaling); Earth: 30°.
        let hadley = (30.0 * (rotation_h.abs().max(1.0) / 24.0).powf(2.0 / 3.0)).clamp(12.0, 80.0);
        let ferrel = (hadley * 2.0).min(85.0);
        Self { itcz: season * 0.3 * obliquity.min(45.0), hadley, ferrel, jets: jets_for(rotation_h, radius_km) }
    }

    /// Cells per hemisphere: the Hadley cell, then a pair per jet; with no jet, a polar cap.
    #[must_use]
    pub fn count(&self) -> u32 {
        if self.jets == 0 {
            2
        } else {
            1 + 2 * self.jets
        }
    }

    /// Where `lat` (poleward of the Hadley cell) falls among the jets: westerly band or not,
    /// position across the band (0–1), and the band's width (degrees). One jet's westerlies run
    /// from the Hadley edge to the polar front, as on Earth; more jets repeat that pair of
    /// bands, westerly and easterly, in the same proportions out to the pole; no jet leaves
    /// one easterly polar band.
    fn band(&self, lat: f64) -> (bool, f64, f64) {
        let e0 = self.edge(lat).abs();
        let span = (90.0 - e0).max(1.0);
        if self.jets == 0 {
            return (false, ((lat.abs() - e0) / span).clamp(0.0, 1.0), span);
        }
        let f = ((self.polar_front(lat) - e0) / span).clamp(0.2, 0.8);
        let pair = span / f64::from(self.jets);
        let r = (lat.abs() - e0).clamp(0.0, span - 1e-9);
        let q = r - (r / pair).floor() * pair;
        if q < f * pair {
            (true, q / (f * pair), f * pair)
        } else {
            (false, (q - f * pair) / ((1.0 - f) * pair), (1.0 - f) * pair)
        }
    }

    /// Centres of the westerly bands on the side of `lat` (where the storm tracks run), and
    /// the band width, degrees.
    fn storm_tracks(&self, lat: f64) -> (Vec<f64>, f64) {
        let e0 = self.edge(lat).abs();
        let span = (90.0 - e0).max(1.0);
        let f = ((self.polar_front(lat) - e0) / span).clamp(0.2, 0.8);
        let pair = span / f64::from(self.jets.max(1));
        ((0..self.jets).map(|k| e0 + f64::from(k) * pair + 0.5 * f * pair).collect(), f * pair)
    }

    /// Subtropical edge of the Hadley cell on the side of `lat`.
    fn edge(&self, lat: f64) -> f64 {
        if lat >= self.itcz {
            self.itcz + self.hadley
        } else {
            self.itcz - self.hadley
        }
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
        if self.jets != 1 {
            let (westerly, x, _) = self.band(lat);
            let s = (PI * x).sin();
            return if westerly { (1.0 * s, 0.3 * s * hemi) } else { (-0.6 * s, -0.35 * s * hemi) };
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
        let belt = |poleward: f64, edge: f64, winter: f64| (0.525 + 0.225 * winter.clamp(-1.0, 1.0)) * plateau(poleward, edge - 0.55 * self.hadley, edge + 0.1 * self.hadley, 8.0);
        let north = belt(lat, self.itcz + self.hadley, -self.itcz / 7.0);
        let south = belt(-lat, self.hadley - self.itcz, self.itcz / 7.0);
        let mut sink = north.max(south);
        // With more jets, each easterly band hands over to the next westerlies under a high.
        if self.jets > 1 {
            let (tracks, w) = self.storm_tracks(lat);
            for c in tracks.iter().skip(1) {
                let high = c - 0.5 * w;
                sink = sink.max(0.4 * (-((lat.abs() - high) / (0.3 * w)).powi(2)).exp());
            }
        }
        (1.0 + ascent - sink + self.storminess(lat)).max(0.2)
    }

    /// Mid-latitude storm tracks, and frequent frontal weather poleward of them (0–1.3). The
    /// storm tracks follow the sun: poleward in summer, equatorward in winter, bringing winter
    /// rain to the subtropical west coasts.
    #[must_use]
    pub fn storminess(&self, lat: f64) -> f64 {
        if self.jets == 0 {
            // No baroclinic eddies to speak of on a slow rotator.
            return 0.0;
        }
        if self.jets > 1 {
            let (tracks, w) = self.storm_tracks(lat);
            let scale = w / 34.0;
            let storm = tracks.iter().map(|c| 0.8 * (-((lat.abs() - c) / (12.0 * scale)).powi(2)).exp()).fold(0.0, f64::max);
            return storm + 0.5 * plateau(lat.abs(), tracks[0], 200.0, 12.0 * scale);
        }
        let storm_lat = ((self.hadley + 18.0).min(self.ferrel) + 0.5 * self.itcz * lat.signum()).max(self.hadley);
        let storm = 0.8 * (-((lat.abs() - storm_lat) / 12.0).powi(2)).exp();
        // Fronts come often poleward of the track, fading in over ~12° (a step here drew
        // rain lines along the latitude circles).
        let subpolar = 0.5 * plateau(lat.abs(), storm_lat, 200.0, 12.0);
        storm + subpolar
    }
}

/// The circulation of a tidally locked world, organised around the point under its star.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Substellar {
    /// Longitude of the substellar point (on the equator), degrees.
    pub lon: f64,
    /// Surface equatorial superrotation, in trade-wind units: none for periods over 20 days.
    pub superrotation: f64,
}

/// Longitude under the star of a locked world: the map is centred on its day side.
pub const SUBSTELLAR_LON: f64 = 0.0;

impl Substellar {
    /// `Some` when the body turns once per orbit (to within 1%).
    #[must_use]
    pub fn for_body(body: &BodyContext) -> Option<Self> {
        let locked = body.is_tidally_locked();
        let period_days = body.exterior.rotation_h / 24.0;
        locked.then(|| Self { lon: SUBSTELLAR_LON, superrotation: 0.4 * ((20.0 - period_days) / 15.0).clamp(0.0, 1.0) })
    }

    /// Angle from the substellar point, degrees (0 under the star, 180 opposite).
    #[must_use]
    pub fn zenith(&self, lat: f64, lon: f64) -> f64 {
        (lat.to_radians().cos() * (lon - self.lon).to_radians().cos()).clamp(-1.0, 1.0).acos().to_degrees()
    }

    /// Surface wind (east, north): toward the substellar point, fastest at the terminator,
    /// plus the equatorial jet.
    #[must_use]
    pub fn wind(&self, lat: f64, lon: f64) -> (f64, f64) {
        let (phi, dl) = (lat.to_radians(), (self.lon - lon).to_radians());
        let (e, n) = (dl.sin(), -phi.sin() * dl.cos());
        let norm = e.hypot(n);
        let speed = 0.9 * self.zenith(lat, lon).to_radians().sin();
        let jet = self.superrotation * (-(lat / 25.0).powi(2)).exp();
        if norm < 1e-9 {
            return (jet, 0.0);
        }
        (speed * e / norm + jet, speed * n / norm)
    }

    /// Rain-out multiplier: a wide convective region under the star, sinking air elsewhere.
    #[must_use]
    pub fn convergence(&self, lat: f64, lon: f64) -> f64 {
        0.3 + 1.5 * (-(self.zenith(lat, lon) / 55.0).powi(2)).exp()
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
    /// A tidally locked world's pattern, which replaces the cells.
    pub locked: Option<Substellar>,
    /// Longitude of each grid column, degrees (used when locked).
    pub lon: Vec<f64>,
}

impl Season {
    /// The season of `body`: [`Season::new`], or for a tidally locked world one season all
    /// year round around its substellar point.
    #[must_use]
    pub fn for_body(a: &Atmo, body: &BodyContext, sign: f64, monsoon: f64) -> Self {
        let (rot, r, obl) = (body.exterior.rotation_h, body.exterior.radius_km, body.irradiance.obliquity_deg);
        match Substellar::for_body(body) {
            Some(s) => Self { cells: Cells::for_body(rot, r, 0.0, 0.0), sign: 0.0, itcz: vec![0.0; a.g.w], locked: Some(s), lon: (0..a.g.w).map(|x| a.g.lon(x)).collect() },
            None => Self::new(a, rot, r, obl, sign, monsoon),
        }
    }

    /// The ITCZ sits at the zonal latitude over open ocean and moves into the summer
    /// hemisphere over land, in proportion to the land (high land counts more: it heats the
    /// air higher up) on the summer side of the tropics. `monsoon` scales the shift
    /// (0 = zonal everywhere; 0.7 = Earth).
    #[must_use]
    pub fn new(a: &Atmo, rotation_h: f64, radius_km: f64, obliquity: f64, sign: f64, monsoon: f64) -> Self {
        let cells = Cells::for_body(rotation_h, radius_km, obliquity, sign);
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
                if w > 0.0 {
                    (s / w).min(1.2)
                } else {
                    0.0
                }
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
        Self { cells, sign, itcz, locked: None, lon: Vec::new() }
    }

    /// The zonal pattern alone (no land): for winds without a surface, and tests.
    #[must_use]
    pub fn zonal(cells: Cells, sign: f64, columns: usize) -> Self {
        Self { cells, sign, itcz: vec![cells.itcz; columns], locked: None, lon: Vec::new() }
    }

    #[must_use]
    pub fn wind(&self, x: usize, lat: f64) -> (f64, f64) {
        match &self.locked {
            Some(s) => s.wind(lat, self.lon[x]),
            None => self.cells.wind_at(lat, self.itcz[x]),
        }
    }

    #[must_use]
    pub fn convergence(&self, x: usize, lat: f64, land: f64) -> f64 {
        match &self.locked {
            Some(s) => s.convergence(lat, self.lon[x]),
            None => self.cells.convergence_at(lat, self.itcz[x], land),
        }
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
        let earth = Cells::for_body(24.0, 6371.0, 23.0, 0.0);
        let slow = Cells::for_body(240.0, 6371.0, 23.0, 0.0);
        assert_eq!(earth.count(), 3);
        assert!((earth.hadley - 30.0).abs() < 1e-9);
        assert!(slow.hadley > 60.0);
        // Trades in the tropics, westerlies in mid-latitudes.
        assert!(earth.wind(15.0).0 < -0.5 && earth.wind(-15.0).0 < -0.5);
        assert!(earth.wind(45.0).0 > 0.5 && earth.wind(-45.0).0 > 0.5);
        // Wet at the ITCZ, dry under the subtropical highs.
        assert!(earth.convergence(0.0) > 1.5 && earth.convergence(22.0) < 0.6 && earth.convergence(-22.0) < 0.6);
    }

    /// Cells by rotation and size: Earth three per hemisphere, Venus and Titan two (a wide
    /// Hadley cell and a polar cap), a world spinning four times faster seven, with
    /// alternating surface westerlies and easterlies, and a storm track under each jet.
    #[test]
    fn more_jets_on_faster_and_bigger_worlds() {
        assert_eq!(jets_for(24.6, 3390.0), 1, "Mars");
        assert_eq!(jets_for(5832.0, 6052.0), 0, "Venus");
        assert_eq!(jets_for(383.0, 2575.0), 0, "Titan");
        assert_eq!(jets_for(6.0, 6371.0), 3, "four times Earth's spin");
        assert_eq!(jets_for(24.0, 1.8 * 6371.0), 2, "a big super-Earth");
        let fast = Cells::for_body(6.0, 6371.0, 0.0, 0.0);
        assert_eq!(fast.count(), 7);
        let signs: Vec<bool> = (13..=89).map(|l| fast.wind(f64::from(l)).0 > 0.0).collect();
        let flips = signs.windows(2).filter(|w| w[0] != w[1]).count();
        assert!(flips >= 4, "alternating bands: {flips} flips");
        let peaks = (13..=88).filter(|&l| {
            let s = |x: i32| fast.storminess(f64::from(x));
            s(l) > s(l - 1) && s(l) >= s(l + 1)
        });
        assert_eq!(peaks.count(), 3, "a storm track per jet");
        let venus = Cells::for_body(5832.0, 6052.0, 0.0, 0.0);
        assert_eq!(venus.count(), 2, "Hadley cell to a polar cap");
        assert!(venus.wind(85.0).0 <= 0.0 && venus.storminess(50.0) == 0.0);
    }

    #[test]
    fn locked_worlds_flow_to_the_substellar_point() {
        let s = Substellar { lon: 0.0, superrotation: 0.0 };
        let east_of = s.wind(0.0, 60.0);
        let north_of = s.wind(40.0, 0.0);
        assert!(east_of.0 < -0.5 && north_of.1 < -0.3, "{east_of:?} {north_of:?}");
        assert!(s.wind(0.0, 179.0).0.abs() < 0.1, "calm at the antistellar point");
        assert!(s.convergence(0.0, 0.0) > 1.7 && s.convergence(0.0, 120.0) < 0.5);
    }

    /// A continent filling the northern subtropics of half the globe pulls July's ITCZ north
    /// over it, and the trades that cross the equator to reach it blow from the west.
    #[test]
    fn monsoon_itcz_follows_land() {
        let g = AtmoGrid::new(2.0);
        let land: Vec<f64> = (0..g.len()).map(|i| if i % g.w < g.w / 2 && g.lat(i / g.w) > 5.0 && g.lat(i / g.w) < 40.0 { 1.0 } else { 0.0 }).collect();
        let a = Atmo { g, wet: true, land, elev: vec![0.0; g.len()] };
        let jul = Season::new(&a, 24.0, 6371.0, 23.4, 1.0, 0.7);
        let (over_land, over_sea) = (jul.itcz[g.w / 4], jul.itcz[3 * g.w / 4]);
        assert!(over_land > 15.0 && (over_sea - 7.0).abs() < 0.5, "land {over_land} sea {over_sea}");
        let (u, v) = jul.wind(g.w / 4, 8.0);
        assert!(u > 0.2 && v > 0.0, "monsoon westerlies {u} {v}");
        let jan = Season::new(&a, 24.0, 6371.0, 23.4, -1.0, 0.7);
        assert!((jan.itcz[g.w / 4] + 7.0).abs() < 0.5, "no southern land, zonal ITCZ in January");
    }
}
