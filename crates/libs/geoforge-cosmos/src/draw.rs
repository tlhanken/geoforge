//! Deterministic draws and names from a seed.

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
    pub(crate) fn normal(&mut self, mean: f64, sd: f64) -> f64 {
        let u1 = self.u().max(1e-12);
        let u2 = self.u();
        mean + sd * (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
}

/// A pronounceable name from a seed.
#[must_use]
pub fn name(seed: Seed) -> String {
    const ON: [&str; 16] = ["k", "t", "v", "s", "r", "m", "n", "l", "d", "th", "z", "b", "h", "sh", "g", "p"];
    const NU: [&str; 8] = ["a", "e", "i", "o", "u", "ae", "ai", "or"];
    let mut d = Draw::new(seed.child(b"name"));
    let syll = 2 + (d.u() * 2.0) as usize;
    let mut s = String::new();
    for _ in 0..syll {
        s.push_str(ON[(d.u() * ON.len() as f64) as usize % ON.len()]);
        s.push_str(NU[(d.u() * NU.len() as f64) as usize % NU.len()]);
    }
    let mut c = s.chars();
    c.next().map(|f| f.to_uppercase().collect::<String>() + c.as_str()).unwrap_or_default()
}
