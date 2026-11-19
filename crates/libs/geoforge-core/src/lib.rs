//! Geoforge v2 core: identity, seeds, geometry, grids, the five data shapes, body and
//! product types. Data only: generation algorithms live in other crates.

pub mod body;
pub mod derived;
pub mod geom;
pub mod grid;
pub mod key;
pub mod product;
pub mod radiation;
pub mod seed;
pub mod shapes;
pub mod surface;

pub use body::{Archetype, BodyContext, BodyExterior, BodyInterior, Elements, SolarSystem, SystemRoot, TectonicRegime};
pub use geom::{LatLon, UnitVec3, Vec3};
pub use grid::{AnyGrid, Cell, EquirectGrid, Grid, GridSpec};
pub use key::{Address, EntityKind, Key, KeyAllocator};
pub use product::{Generator, LayerId, ProductData, ProductId, Products, RunContext, Source, Status};
pub use seed::Seed;
pub use shapes::{Edge, Event, EventKind, Field, Network, Partition, Place, Table, Timeline};

/// Current schema version of every stored document.
pub const SCHEMA_VERSION: u32 = 1;

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("invalid: {0}")]
    Invalid(String),
    #[error("missing input: {0}")]
    MissingInput(ProductId),
    #[error("not applicable: {0}")]
    NotApplicable(product::Precondition),
    #[error("unsupported: {0}")]
    Unsupported(String),
}
