//! The atmosphere model shared by the circulation and climate layers.
//!
//! | Module | Concern |
//! |---|---|
//! | [`grid`] | the coarse latitude–longitude grid and the world's surface on it |
//! | [`transport`] | carrying quantities with the wind |
//! | [`cells`] | Hadley, Ferrel and polar cells; monsoons |
//! | [`radiation`] | energy balance by latitude, seasons |

pub mod cells;
pub mod grid;
pub mod radiation;
pub mod transport;

pub use cells::{winds, Cells};
pub use grid::{Atmo, AtmoGrid};
pub use radiation::Radiation;
