//! Generators for rigid-lid worlds: stagnant (Mars, Mercury, the Moon), episodic (Venus)
//! and geologically dead bodies. They fill the tectonics, geology and elevation layers with
//! the same products as the plate-tectonic generators, so later layers work unchanged.
//!
//! | Module | Layer |
//! |---|---|
//! | [`tectonics`] | one plate, highland–lowland crust, surface age |
//! | [`geology`] | basins, volcanic rises and shields, rifts, coronae, tesserae; lithology |
//! | [`elevation`] | the step, the landforms, craters, roughness |
//! | [`landforms`] | shapes and the crater population shared by the above |

pub mod elevation;
pub mod geology;
pub mod landforms;
pub mod tectonics;

pub use elevation::{LidElevation, LidElevationParams};
pub use geology::{LidGeology, LidGeologyParams};
pub use tectonics::{LidTectonics, LidTectonicsParams};
