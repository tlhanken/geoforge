//! Deterministic draws from a seed.

use geoforge_core::seed::{range, unit, Seed};

/// Deterministic draws from a seed: each call uses a fresh labelled child.
pub(crate) struct Draw {
    seed: Seed,
    n: u64,
}

impl Draw {
    pub(crate) fn new(seed: Seed) -> Self {
        Self { seed, n: 0 }
    }
    pub(crate) fn u(&mut self) -> f64 {
        self.n += 1;
        unit(self.seed.child_u64(b"draw", self.n))
    }
    pub(crate) fn range(&mut self, lo: f64, hi: f64) -> f64 {
        self.n += 1;
        range(self.seed.child_u64(b"draw", self.n), lo, hi)
    }
    /// A whole number from `lo` to `hi`, both included, equally likely.
    pub(crate) fn int(&mut self, lo: u32, hi: u32) -> u32 {
        lo + ((self.u() * f64::from(hi - lo + 1)) as u32).min(hi - lo)
    }
    pub(crate) fn normal(&mut self, mean: f64, sd: f64) -> f64 {
        let u1 = self.u().max(1e-12);
        let u2 = self.u();
        mean + sd * (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }

    /// A Rayleigh deviate with scale `sigma` (its mode).
    pub(crate) fn rayleigh(&mut self, sigma: f64) -> f64 {
        sigma * (-2.0 * self.u().max(1e-12).ln()).sqrt()
    }

    /// A Poisson count with mean `mean`.
    pub(crate) fn poisson(&mut self, mean: f64) -> usize {
        if mean <= 0.0 {
            return 0;
        }
        if mean > 50.0 {
            return self.normal(mean, mean.sqrt()).round().max(0.0) as usize;
        }
        let (limit, mut p, mut k) = ((-mean).exp(), 1.0, 0);
        loop {
            p *= self.u();
            if p <= limit {
                return k;
            }
            k += 1;
        }
    }

    /// A direction drawn evenly over the sphere.
    pub(crate) fn direction(&mut self) -> [f64; 3] {
        let z = 2.0 * self.u() - 1.0;
        let phi = std::f64::consts::TAU * self.u();
        let s = (1.0 - z * z).sqrt();
        [s * phi.cos(), s * phi.sin(), z]
    }

    /// A rotation drawn evenly over all orientations (Shoemake 1992), as a matrix.
    pub(crate) fn rotation(&mut self) -> [[f64; 3]; 3] {
        let (u1, u2, u3) = (self.u(), self.u(), self.u());
        let tau = std::f64::consts::TAU;
        let (a, b) = ((1.0 - u1).sqrt(), u1.sqrt());
        let (x, y, z, w) = (a * (tau * u2).sin(), a * (tau * u2).cos(), b * (tau * u3).sin(), b * (tau * u3).cos());
        [
            [1.0 - 2.0 * (y * y + z * z), 2.0 * (x * y - z * w), 2.0 * (x * z + y * w)],
            [2.0 * (x * y + z * w), 1.0 - 2.0 * (x * x + z * z), 2.0 * (y * z - x * w)],
            [2.0 * (x * z - y * w), 2.0 * (y * z + x * w), 1.0 - 2.0 * (x * x + y * y)],
        ]
    }

    /// A standard normal deviate.
    pub(crate) fn gauss(&mut self) -> f64 {
        self.normal(0.0, 1.0)
    }
}
