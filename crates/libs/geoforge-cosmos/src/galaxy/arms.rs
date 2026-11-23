//! Spiral arms: logarithmic spirals that gather a disk's young stars, gas and dust.
//!
//! Arm `k` of `m` lies at azimuth φ₀ + 2πk/m + ln(R/R₀)/tan p (Vallée 1995's model), so the
//! arms trail a disk turning clockwise seen from its north pole, as the Milky Way does. A
//! star's density is raised along the arms and lowered between them by a Gaussian of the
//! distance across the nearest arms, keeping each ring's mass: 1 + c (g/ḡ − 1), where g sums
//! the arms' Gaussians and ḡ is its mean around the ring.

use serde::{Deserialize, Serialize};
use std::f64::consts::{PI, TAU};

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Arms {
    pub count: u32,
    /// Pitch angle, degrees: small is tightly wound.
    pub pitch_deg: f64,
    /// Where the arms begin, kpc.
    pub start_kpc: f64,
    /// Azimuth of the first arm at its start, degrees (the x axis is 0).
    pub phase_deg: f64,
    /// An arm's half-width (the Gaussian's σ across it), kpc.
    pub width_kpc: f64,
}

/// Over this distance inside their start the arms fade out, kpc.
const FADE_KPC: f64 = 0.5;

impl Arms {
    /// How much the arms raise (above 1) or lower the density at (`x`, `y`) pc, for stars
    /// whose arm contrast is `contrast` (0 none, 1 all of a ring's mass in its arms).
    #[must_use]
    pub fn factor(&self, x: f64, y: f64, contrast: f64) -> f64 {
        if contrast <= 0.0 {
            return 1.0;
        }
        (1.0 + contrast * self.excess(x, y)).max(0.0)
    }

    /// The arms' excess at (`x`, `y`) pc for full contrast: the density factor is
    /// 1 + contrast × this, so one excess serves stars of every contrast.
    #[must_use]
    pub fn excess(&self, x: f64, y: f64) -> f64 {
        let r = x.hypot(y) / 1000.0;
        let fade = ((r - self.start_kpc + FADE_KPC) / FADE_KPC).clamp(0.0, 1.0);
        if fade <= 0.0 {
            return 0.0;
        }
        let (g, mean) = self.gather(r, y.atan2(x));
        fade * (g / mean - 1.0)
    }

    /// The arms' Gaussians summed at radius `r` kpc and azimuth `phi`, and their mean
    /// around that ring.
    fn gather(&self, r: f64, phi: f64) -> (f64, f64) {
        let m = f64::from(self.count.max(1));
        let sin_p = self.pitch_deg.to_radians().sin();
        let tan_p = self.pitch_deg.to_radians().tan();
        let base = self.phase_deg.to_radians() + (r.max(1e-3) / self.start_kpc).ln() / tan_p;
        let mut g = 0.0;
        for k in 0..self.count.max(1) {
            let dphi = (phi - base - TAU * f64::from(k) / m + PI).rem_euclid(TAU) - PI;
            let across = r * dphi * sin_p;
            g += (-0.5 * (across / self.width_kpc).powi(2)).exp();
        }
        // ∫ over the ring of m Gaussians of σ_φ = w/(R sin p), divided by 2π.
        let mean = (m * self.width_kpc * (TAU).sqrt() / (TAU * r.max(1e-3) * sin_p)).min(1.0);
        (g, mean)
    }

    /// Distance from (`x`, `y`) pc to the nearest arm, kpc.
    #[must_use]
    pub fn distance_kpc(&self, x: f64, y: f64) -> f64 {
        let r = x.hypot(y) / 1000.0;
        let m = f64::from(self.count.max(1));
        let base = self.phase_deg.to_radians() + (r.max(1e-3) / self.start_kpc).ln() / self.pitch_deg.to_radians().tan();
        (0..self.count.max(1)).map(|k| ((y.atan2(x) - base - TAU * f64::from(k) / m + PI).rem_euclid(TAU) - PI).abs() * r * self.pitch_deg.to_radians().sin()).fold(f64::MAX, f64::min)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn milky_way() -> Arms {
        Arms { count: 4, pitch_deg: 13.1, start_kpc: 2.2, phase_deg: -9.0, width_kpc: 0.4 }
    }

    /// A ring keeps its mass, the arms hold more of it than the gaps, and inside their start
    /// there are none.
    #[test]
    fn arms_move_mass_around_a_ring() {
        let a = milky_way();
        for r in [8200.0, 12_000.0] {
            let n = 3600;
            let f: Vec<f64> = (0..n)
                .map(|i| {
                    let phi = TAU * f64::from(i) / f64::from(n);
                    a.factor(r * phi.cos(), r * phi.sin(), 0.8)
                })
                .collect();
            let mean = f.iter().sum::<f64>() / f64::from(n);
            assert!((mean - 1.0).abs() < 0.03, "{r}: {mean}");
            let (lo, hi) = f.iter().fold((f64::MAX, f64::MIN), |(lo, hi), &v| (lo.min(v), hi.max(v)));
            assert!(hi > 1.8 && lo < 0.4, "{r}: {lo}–{hi}");
        }
        // Further in, the arms crowd together and blur.
        let inner: Vec<f64> = (0..360).map(|i| a.factor(4000.0 * f64::from(i).to_radians().cos(), 4000.0 * f64::from(i).to_radians().sin(), 0.8)).collect();
        assert!(inner.iter().copied().fold(0.0, f64::max) > 1.2);
        assert_eq!(a.factor(1000.0, 0.0, 0.8), 1.0);
        assert_eq!(a.factor(8200.0, 0.0, 0.0), 1.0);
    }

    /// The Milky Way's phase puts the Sun (8.2 kpc on the x axis) between two arms.
    #[test]
    fn the_sun_lies_between_arms() {
        let a = milky_way();
        let gap = PI / 4.0 * 8.2 * 13.1f64.to_radians().sin();
        assert!((a.distance_kpc(8200.0, 0.0) - gap).abs() < 0.15, "{}", a.distance_kpc(8200.0, 0.0));
    }
}
