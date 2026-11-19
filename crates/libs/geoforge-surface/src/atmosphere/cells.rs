//! The general circulation by rotation: Hadley, Ferrel and polar cells, and monsoons.

use super::grid::Atmo;
use std::f64::consts::PI;

/// Hadley, Ferrel and polar cells of one season.
#[derive(Debug, Clone, Copy)]
pub struct Cells {
    /// Latitude of the ITCZ, degrees.
    pub itcz: f64,
    /// Width of each Hadley cell, degrees from the ITCZ.
    pub hadley: f64,
    /// Poleward edge of the Ferrel cells, degrees.
    pub ferrel: f64,
}

impl Cells {
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
        // The storm tracks follow the sun: poleward in summer, equatorward in winter, bringing
        // winter rain to the subtropical west coasts.
        let storm_lat = ((self.hadley + 18.0).min(self.ferrel) + 1.2 * self.itcz * lat.signum()).max(self.hadley);
        let storm = 0.8 * (-((lat.abs() - storm_lat) / 12.0).powi(2)).exp();
        let subpolar = if lat.abs() > storm_lat { 0.5 } else { 0.0 };
        (1.0 + itcz - high + storm + subpolar).max(0.2)
    }
}

/// Winds for one season: cell pattern plus monsoon flow toward (summer) or away from
/// (winter) large heated land masses.
pub fn winds(a: &Atmo, circ: &Cells, season: f64, monsoon: f64) -> Vec<(f64, f64)> {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hadley_cell_scales_with_rotation() {
        let earth = Cells::for_body(24.0, 23.0, 0.0);
        let slow = Cells::for_body(240.0, 23.0, 0.0);
        assert!((earth.hadley - 30.0).abs() < 1e-9);
        assert!(slow.hadley > 60.0);
        // Trades in the tropics, westerlies in mid-latitudes.
        assert!(earth.wind(15.0).0 < -0.5 && earth.wind(-15.0).0 < -0.5);
        assert!(earth.wind(45.0).0 > 0.5 && earth.wind(-45.0).0 > 0.5);
    }
}
