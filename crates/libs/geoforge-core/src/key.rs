//! Stable identity: keys (addresses are in [`crate::address`]).

use serde::{Deserialize, Serialize};
use std::fmt;

/// Stable identity of an entity. Assigned at creation, never reused, never derived from a
/// list position or a name.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(transparent)]
pub struct Key(pub u64);

impl fmt::Display for Key {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// Hands out fresh keys for a container. Generated entities can instead take a key from
/// their procedural slot so that regeneration reproduces them.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct KeyAllocator {
    next: u64,
}

impl KeyAllocator {
    #[must_use]
    pub fn starting_at(next: u64) -> Self {
        Self { next }
    }

    pub fn fresh(&mut self) -> Key {
        self.next += 1;
        Key(self.next)
    }

    /// Make sure later fresh keys never collide with `key`.
    pub fn reserve(&mut self, key: Key) {
        self.next = self.next.max(key.0);
    }
}
