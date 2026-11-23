//! Generators for ice shells floating on an ocean: thin, young, broken shells like Europa's
//! and thick, old, grooved ones like Ganymede's. They fill the tectonics, geology and
//! elevation layers with the same products as the plate and rigid-lid models, so later layers
//! work unchanged.
//!
//! | Module | Layer |
//! |---|---|
//! | [`tectonics`] | one shell; old plains and renewed ice; the surface's age |
//! | [`geology`] | ridged plains, chaos, bands and ridges, or dark cratered and bright grooved terrain; basins; lithology |
//! | [`elevation`] | low relief: chaos blocks, raised ridges and bands, grooves, shallow craters |
//!
//! Between a shell 25 km thick and one 60 km thick the style passes from the one to the
//! other. Titan's own surface (dunes, methane lakes and rivers) is not modelled: its shell is
//! drawn as a thick one.

pub mod elevation;
pub mod geology;
pub mod tectonics;

pub use elevation::{IceShellElevation, IceShellElevationParams};
pub use geology::{IceShellGeology, IceShellGeologyParams};
pub use tectonics::{IceShellTectonics, IceShellTectonicsParams, ShellStyle};
