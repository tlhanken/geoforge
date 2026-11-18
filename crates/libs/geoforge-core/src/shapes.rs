//! The five data shapes every product takes. Storage, import, export, patching and
//! validation are written once per shape, not once per layer.

use crate::geom::UnitVec3;
use crate::grid::{Cell, GridSpec};
use crate::key::Key;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// Entities with identity and attributes, ordered by key.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(transparent)]
pub struct Table<E> {
    pub rows: BTreeMap<Key, E>,
}

impl<E> Default for Table<E> {
    fn default() -> Self {
        Self { rows: BTreeMap::new() }
    }
}

impl<E> Table<E> {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }
    pub fn insert(&mut self, key: Key, row: E) {
        self.rows.insert(key, row);
    }
    #[must_use]
    pub fn get(&self, key: Key) -> Option<&E> {
        self.rows.get(&key)
    }
    pub fn get_mut(&mut self, key: Key) -> Option<&mut E> {
        self.rows.get_mut(&key)
    }
    #[must_use]
    pub fn len(&self) -> usize {
        self.rows.len()
    }
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.rows.is_empty()
    }
    pub fn iter(&self) -> impl Iterator<Item = (Key, &E)> {
        self.rows.iter().map(|(k, v)| (*k, v))
    }
    #[must_use]
    pub fn max_key(&self) -> Key {
        self.rows.keys().next_back().copied().unwrap_or(Key(0))
    }
}

/// One value per grid cell.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Field<T> {
    pub grid: GridSpec,
    pub values: Vec<T>,
}

impl<T: Clone> Field<T> {
    #[must_use]
    pub fn filled(grid: GridSpec, value: T) -> Self {
        Self { grid, values: vec![value; grid.cell_count()] }
    }
}

impl<T> Field<T> {
    #[must_use]
    pub fn from_values(grid: GridSpec, values: Vec<T>) -> Self {
        assert_eq!(values.len(), grid.cell_count(), "field length must match its grid");
        Self { grid, values }
    }
    #[must_use]
    pub fn get(&self, c: Cell) -> &T {
        &self.values[c as usize]
    }
    pub fn set(&mut self, c: Cell, v: T) {
        self.values[c as usize] = v;
    }
    #[must_use]
    pub fn len(&self) -> usize {
        self.values.len()
    }
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }
}

/// An owner per cell plus the owners' table: plates, provinces, watersheds, kingdoms.
///
/// `owner` holds the owning row's key (keys of a partition must fit in `u32`); `0` means no
/// owner. Invariant, checked by [`Partition::validate`]: every non-zero owner exists in `table`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Partition<E> {
    pub owner: Field<u32>,
    pub table: Table<E>,
}

impl<E> Partition<E> {
    #[must_use]
    pub fn owner_of(&self, c: Cell) -> Option<Key> {
        match *self.owner.get(c) {
            0 => None,
            k => Some(Key(u64::from(k))),
        }
    }

    /// Cells per owner, ascending.
    #[must_use]
    pub fn cells_by_owner(&self) -> BTreeMap<Key, Vec<Cell>> {
        let mut m: BTreeMap<Key, Vec<Cell>> = BTreeMap::new();
        for (c, &o) in self.owner.values.iter().enumerate() {
            if o != 0 {
                m.entry(Key(u64::from(o))).or_default().push(c as Cell);
            }
        }
        m
    }

    /// # Errors
    /// Lists owners that are missing from the table.
    pub fn validate(&self) -> Result<(), crate::Error> {
        let mut missing = std::collections::BTreeSet::new();
        for &o in &self.owner.values {
            if o != 0 && !self.table.rows.contains_key(&Key(u64::from(o))) {
                missing.insert(o);
            }
        }
        if missing.is_empty() {
            Ok(())
        } else {
            Err(crate::Error::Invalid(format!("partition owners missing from table: {missing:?}")))
        }
    }

    /// Drop table rows that own no cells.
    pub fn prune(&mut self) {
        let used: std::collections::BTreeSet<u32> = self.owner.values.iter().copied().collect();
        self.table.rows.retain(|k, _| used.contains(&(k.0 as u32)));
    }
}

/// A graph over the surface: rivers, roads, plate boundaries.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Network<N, E> {
    pub nodes: Table<N>,
    pub edges: Vec<Edge<E>>,
}

impl<N, E> Default for Network<N, E> {
    fn default() -> Self {
        Self { nodes: Table::default(), edges: Vec::new() }
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Edge<E> {
    pub from: Key,
    pub to: Key,
    pub attrs: E,
}

/// Events in time, with causes.
#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
#[serde(transparent)]
pub struct Timeline {
    pub events: Table<Event>,
}

/// One event, `t_ma` million years before the world's present.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Event {
    pub t_ma: f64,
    pub kind: EventKind,
    pub at: Place,
    pub magnitude: f64,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub causes: Vec<Key>,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub note: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EventKind {
    CratonFormation,
    PaleoOrogeny,
    RiftOnset,
    OceanOpening,
    SubductionOnset,
    ContinentalCollision,
    LargeIgneousProvince,
    HotspotTrack,
    TransformInitiation,
}

/// Where an event happened.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum Place {
    /// A point on the body.
    Point { at: UnitVec3 },
    /// A row of another product (a plate, a province, a boundary segment).
    Entity { product: crate::product::ProductId, key: Key },
    /// Two plates' shared boundary.
    PlatePair { a: Key, b: Key, at: UnitVec3 },
}
