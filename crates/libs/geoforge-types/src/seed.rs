//! Deterministic seed derivation for hierarchical procedural generation.

use serde::{Deserialize, Serialize};

/// Root or child seed for reproducible procedural generation.
///
/// Child seeds are derived with SplitMix64-style mixing so sibling branches
/// are independent but fully determined by the parent.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(transparent)]
pub struct Seed(pub u64);

impl Seed {
    /// Create a seed from an explicit value.
    #[must_use]
    pub const fn new(value: u64) -> Self {
        Self(value)
    }

    /// Raw seed value.
    #[must_use]
    pub const fn value(self) -> u64 {
        self.0
    }

    /// Derive a child seed from a stable byte label (e.g. `b"tectonics"`).
    #[must_use]
    pub fn child(self, label: &[u8]) -> Self {
        let mut hash = self.0;
        for &byte in label {
            hash = mix(hash ^ u64::from(byte));
        }
        Self(mix(hash))
    }

    /// Derive a child seed from a label and a coordinate-derived hash.
    #[must_use]
    pub fn child_at(self, label: &[u8], coordinate_hash: u64) -> Self {
        self.child(label).child_bytes(&coordinate_hash.to_le_bytes())
    }

    /// Derive a child seed from arbitrary bytes (e.g. serialized coordinates).
    #[must_use]
    pub fn child_bytes(self, bytes: &[u8]) -> Self {
        let mut hash = self.0;
        for chunk in bytes.chunks(8) {
            let mut buf = [0u8; 8];
            buf[..chunk.len()].copy_from_slice(chunk);
            hash = mix(hash ^ u64::from_le_bytes(buf));
        }
        Self(mix(hash))
    }
}

/// SplitMix64 finalizer — fast, deterministic avalanche.
#[must_use]
pub(crate) fn mix(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = x;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn same_parent_and_label_yield_same_child() {
        let a = Seed::new(42).child(b"geology");
        let b = Seed::new(42).child(b"geology");
        assert_eq!(a, b);
    }

    #[test]
    fn different_labels_yield_different_children() {
        let a = Seed::new(42).child(b"tectonics");
        let b = Seed::new(42).child(b"geology");
        assert_ne!(a, b);
    }

    #[test]
    fn different_parents_yield_different_children() {
        let a = Seed::new(1).child(b"stage");
        let b = Seed::new(2).child(b"stage");
        assert_ne!(a, b);
    }

    #[test]
    fn child_at_depends_on_coordinate_hash() {
        let base = Seed::new(99).child(b"pixel");
        let a = base.child_at(b"pixel", 100);
        let b = base.child_at(b"pixel", 200);
        assert_ne!(a, b);
    }
}
