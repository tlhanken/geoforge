//! Stable identity: keys and addresses.

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

/// What kind of entity a key in an address names.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EntityKind {
    Universe,
    Galaxy,
    Sector,
    System,
    Body,
}

/// Where an entity lives: the chain of keys from the outermost container down.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Default, Serialize, Deserialize)]
#[serde(transparent)]
pub struct Address(pub Vec<(EntityKind, Key)>);

impl Address {
    #[must_use]
    pub fn root(kind: EntityKind, key: Key) -> Self {
        Self(vec![(kind, key)])
    }

    #[must_use]
    pub fn child(&self, kind: EntityKind, key: Key) -> Self {
        let mut v = self.0.clone();
        v.push((kind, key));
        Self(v)
    }

    #[must_use]
    pub fn last(&self) -> Option<(EntityKind, Key)> {
        self.0.last().copied()
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
