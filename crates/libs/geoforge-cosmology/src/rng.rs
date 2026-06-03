//! Deterministic pseudo-random values from [`Seed`].

use geoforge_types::Seed;

/// Unit interval sample from a seed (deterministic).
#[must_use]
pub fn unit(seed: Seed) -> f64 {
    let mixed = geoforge_types::seed::mix(seed.value());
    (mixed as f64) / (u64::MAX as f64)
}

/// Sample in `[min, max)`.
#[must_use]
pub fn range(seed: Seed, min: f64, max: f64) -> f64 {
    min + unit(seed) * (max - min)
}

/// Derive a child seed for a labeled sub-draw.
#[must_use]
pub fn draw(parent: Seed, label: &[u8]) -> Seed {
    parent.child(label)
}
