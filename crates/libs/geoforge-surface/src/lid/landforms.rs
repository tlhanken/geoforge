//! Landform shapes of rigid-lid worlds and the crater population, shared by the lid layers.
//!
//! Profiles are functions of `x = r / R` (distance from the centre over the feature's
//! radius) returning height in units of the feature's depth or height.

use geoforge_core::grid::{cells_within, Cell, Grid};
use geoforge_core::UnitVec3;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

/// Craters this wide or wider are mapped as impact basins (geology provinces); smaller ones
/// are drawn straight into the elevation.
pub const BASIN_MIN_KM: f64 = 400.0;

/// Cumulative crater density on a 4 Ga surface: craters wider than `D` km per km² is
/// `CRATER_K / D²` (about 2–5 thousand craters over 20 km on the Moon's highlands).
const CRATER_K: f64 = 0.02;

/// Heavily cratered surfaces stop gaining craters (new ones erase old ones).
const SATURATION: f64 = 3.0;

pub use crate::relief::relief_scale;

/// Craters accumulated on a surface of age `age_ma`, relative to a 4 Ga surface: Neukum's
/// production function (a steep early bombardment, then a slow constant rate), saturating.
#[must_use]
pub fn production(age_ma: f64) -> f64 {
    let n = |t_gyr: f64| 5.44e-14 * ((6.93 * t_gyr).exp() - 1.0) + 8.38e-4 * t_gyr;
    (n(age_ma.max(0.0) / 1000.0) / n(4.0)).min(SATURATION)
}

/// A crater or basin: centre and diameter.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Impact {
    pub at: UnitVec3,
    pub diameter_km: f64,
}

/// Craters between `d_min_km` and `d_max_km` across, scattered with the density their
/// ground's age implies (`surface_age_ma` per cell). Largest first.
#[must_use]
pub fn craters(grid: &dyn Grid, surface_age_ma: &[f32], radius_km: f64, d_min_km: f64, d_max_km: f64, seed: u64) -> Vec<Impact> {
    if d_max_km <= d_min_km {
        return Vec::new();
    }
    let km2 = radius_km * radius_km;
    // Expected craters per cell, then a cumulative table to draw cells from.
    let size_term = d_min_km.powi(-2) - d_max_km.powi(-2);
    let mut cum = Vec::with_capacity(grid.len());
    let mut total = 0.0;
    for c in grid.cells() {
        total += grid.area_sr(c) * km2 * CRATER_K * production(f64::from(surface_age_ma[c as usize])) * size_term;
        cum.push(total);
    }
    let mut rng = StdRng::seed_from_u64(seed);
    // Poisson count via its normal approximation (exact enough for the dozens to thousands
    // drawn here; a plain draw for small expectations).
    let count = if total < 30.0 {
        let mut k = 0;
        let mut p = rng.gen::<f64>();
        while p > (-total).exp() {
            k += 1;
            p *= rng.gen::<f64>();
        }
        k
    } else {
        let z: f64 = (0..12).map(|_| rng.gen::<f64>()).sum::<f64>() - 6.0;
        (total + z * total.sqrt()).round().max(0.0) as usize
    };
    let half = grid.spacing_rad() * 0.5;
    let mut out: Vec<Impact> = (0..count)
        .map(|_| {
            let u = rng.gen::<f64>() * total;
            let c = cum.partition_point(|&v| v < u).min(grid.len() - 1) as Cell;
            let at = grid.center(c).travel(rng.gen_range(0.0..360.0), rng.gen::<f64>() * half);
            let diameter_km = (d_min_km.powi(-2) - rng.gen::<f64>() * size_term).powf(-0.5);
            Impact { at, diameter_km }
        })
        .collect();
    out.sort_by(|a, b| b.diameter_km.total_cmp(&a.diameter_km));
    out
}

/// A height field (m) being drawn on.
pub struct Canvas<'a> {
    pub grid: &'a dyn Grid,
    /// The body's radius, km.
    pub radius_km: f64,
    pub height_m: Vec<f64>,
}

impl Canvas<'_> {
    /// Add `height_m · shape(r / radius)` to every cell within `reach × radius` of `at`.
    pub fn stamp(&mut self, at: UnitVec3, radius_km: f64, reach: f64, height_m: f64, shape: impl Fn(f64) -> f64) {
        let r_rad = radius_km / self.radius_km;
        let reach_rad = (r_rad * reach).max(self.grid.spacing_rad());
        for c in cells_within(self.grid, &[self.grid.locate(at)], reach_rad * 1.3, |_| true) {
            let x = self.grid.center(c).angle_to(at) / r_rad;
            if x <= reach {
                self.height_m[c as usize] += height_m * shape(x);
            }
        }
    }
}

/// Complex crater: flat floor, raised rim, ejecta apron, central peak in the larger ones.
#[must_use]
pub fn crater(x: f64, central_peak: bool) -> f64 {
    let floor = -(1.6 * (1.0 - x * x)).clamp(0.0, 1.0);
    let rim = 0.3 * (-((x - 1.0) / 0.18).powi(2)).exp();
    let ejecta = if x > 1.0 { 0.08 * x.powi(-3) } else { 0.0 };
    let peak = if central_peak { 0.25 * (-(x / 0.15).powi(2)).exp() } else { 0.0 };
    floor + rim + ejecta + peak
}

/// Depth of a complex crater, m, before gravity scaling.
#[must_use]
pub fn crater_depth_m(diameter_km: f64) -> f64 {
    250.0 * diameter_km.powf(0.4)
}

/// Multi-ring impact basin: deep flat floor, main rim, faint outer ring.
#[must_use]
pub fn basin(x: f64) -> f64 {
    let floor = -(1.0 - x.powi(4)).max(0.0);
    let rim = 0.2 * (-((x - 1.0) / 0.2).powi(2)).exp();
    let ring = 0.08 * (-((x - 1.5) / 0.15).powi(2)).exp();
    floor + rim + ring
}

/// Depth of an impact basin, m, before gravity scaling.
#[must_use]
pub fn basin_depth_m(diameter_km: f64) -> f64 {
    (1500.0 + 2.0 * diameter_km).min(8000.0)
}

/// Broad volcanic rise (Tharsis, Beta Regio).
#[must_use]
pub fn dome(x: f64) -> f64 {
    if x < 1.0 { (1.0 - x * x).powi(2) } else { 0.0 }
}

/// Shield volcano: concave flanks, a summit plateau and a caldera sunk into it.
#[must_use]
pub fn shield(x: f64) -> f64 {
    let cone = if x < 1.0 { (1.0 - x.max(0.12)).powf(1.3) } else { 0.0 };
    let caldera = if x < 0.1 { 0.12 * (1.0 - (x / 0.1).powi(2)) } else { 0.0 };
    cone - caldera
}

/// Corona: a raised ring around a sagging centre.
#[must_use]
pub fn corona(x: f64) -> f64 {
    let ring = (-((x - 0.85) / 0.15).powi(2)).exp();
    let sag = if x < 0.7 { -0.4 * (1.0 - (x / 0.7).powi(2)) } else { 0.0 };
    ring + sag
}

#[cfg(test)]
mod tests {
    use super::*;
    use geoforge_core::grid::{Grid, GridSpec};
    use geoforge_core::{Archetype, Key};

    #[test]
    fn production_is_front_loaded() {
        assert!((production(4000.0) - 1.0).abs() < 1e-9);
        assert!(production(3500.0) < 0.2, "{}", production(3500.0));
        assert!(production(500.0) < 0.01);
        assert_eq!(production(4500.0), SATURATION);
    }

    #[test]
    fn gravity_scales_relief() {
        let s = |a: Archetype| relief_scale(&a.context(Key(1), "x"));
        assert!((s(Archetype::Earth) - 1.0).abs() < 0.01);
        assert!((s(Archetype::Mars) - 2.63).abs() < 0.05);
        assert_eq!(s(Archetype::Moon), 3.0);
    }

    #[test]
    fn craters_follow_surface_age() {
        let grid = GridSpec::Equirect { width: 256, height: 128 }.build().unwrap();
        // Old north, young south.
        let age: Vec<f32> = grid.cells().map(|c| if grid.center(c).to_lat_lon().lat_deg > 0.0 { 4300.0 } else { 500.0 }).collect();
        let cs = craters(&grid, &age, 1737.4, 60.0, BASIN_MIN_KM, 7);
        let north = cs.iter().filter(|c| c.at.to_lat_lon().lat_deg > 0.0).count();
        assert!(north > 30 && north > 20 * (cs.len() - north), "north {north} of {}", cs.len());
        assert!(cs.iter().all(|c| (60.0..BASIN_MIN_KM).contains(&c.diameter_km)));
        assert!(cs.windows(2).all(|w| w[0].diameter_km >= w[1].diameter_km), "largest first");
        assert_eq!(cs, craters(&grid, &age, 1737.4, 60.0, BASIN_MIN_KM, 7), "deterministic");
    }

    #[test]
    fn profiles() {
        assert!(crater(0.0, false) < -0.9 && crater(1.0, false) > 0.2 && crater(3.0, false).abs() < 0.01);
        assert!(basin(0.0) < -0.99 && basin(1.0) > 0.15);
        assert!(dome(0.0) == 1.0 && dome(1.2) == 0.0);
        assert!(shield(0.3) > shield(0.6) && shield(0.0) < shield(0.15), "caldera at the summit");
        assert!(corona(0.85) > 0.9 && corona(0.0) < 0.0);
    }
}
