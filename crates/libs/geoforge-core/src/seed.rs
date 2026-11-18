//! Deterministic seed derivation.
//!
//! The mixing function and derivation rules here are **frozen**: changing them changes every
//! generated world. Golden-value tests at the bottom of this file pin them.

use crate::key::Key;
use crate::product::ProductId;
use serde::{Deserialize, Serialize};

/// A 64-bit seed. Child seeds are derived with SplitMix64-style mixing.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(transparent)]
pub struct Seed(pub u64);

impl Seed {
    #[must_use]
    pub const fn new(value: u64) -> Self {
        Self(value)
    }

    #[must_use]
    pub const fn value(self) -> u64 {
        self.0
    }

    /// Child seed from a stable byte label (e.g. `b"plates"`).
    #[must_use]
    pub fn child(self, label: &[u8]) -> Self {
        let mut hash = self.0;
        for &byte in label {
            hash = mix(hash ^ u64::from(byte));
        }
        Self(mix(hash))
    }

    /// Child seed from a label and an integer (an index or a key).
    #[must_use]
    pub fn child_u64(self, label: &[u8], n: u64) -> Self {
        Self(mix(self.child(label).0 ^ mix(n)))
    }

    /// The seed for one product of one owner (a body, a system).
    ///
    /// Depends only on the root seed, the product's stable name and the owner's own key —
    /// never on a list index, a name, or the owner's ancestors. Attaching a thin world's
    /// system to a galaxy sector later therefore never changes any product seed.
    #[must_use]
    pub fn for_product(self, product: ProductId, owner: Key) -> Seed {
        self.child(product.name().as_bytes()).child_u64(b"owner", owner.0)
    }
}

/// SplitMix64 finalizer.
#[must_use]
pub fn mix(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = x;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Deterministic unit-interval sample in `[0, 1)`.
#[must_use]
pub fn unit(seed: Seed) -> f64 {
    (mix(seed.0) >> 11) as f64 / (1u64 << 53) as f64
}

/// Deterministic sample in `[min, max)`.
#[must_use]
pub fn range(seed: Seed, min: f64, max: f64) -> f64 {
    min + unit(seed) * (max - min)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn golden_mix() {
        // Frozen: changing these values is a breaking change to every world.
        assert_eq!(mix(0), 0xE220_A839_7B1D_CDAF);
        assert_eq!(mix(42), 0xBDD7_3226_2FEB_6E95);
    }

    #[test]
    fn golden_child_and_product() {
        let root = Seed(42);
        assert_eq!(root.child(b"plates"), root.child(b"plates"));
        assert_ne!(root.child(b"plates"), root.child(b"geology"));
        let a = root.for_product(ProductId::Plates, Key(1));
        let b = root.for_product(ProductId::Plates, Key(2));
        let c = root.for_product(ProductId::CrustType, Key(1));
        assert_ne!(a, b);
        assert_ne!(a, c);
        assert_eq!(a, Seed(42).for_product(ProductId::Plates, Key(1)));
    }

    #[test]
    fn unit_in_range() {
        for i in 0..1000 {
            let u = unit(Seed(i));
            assert!((0.0..1.0).contains(&u));
        }
    }
}
