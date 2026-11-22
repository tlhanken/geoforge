//! Bodies, orbits and solar systems.
//!
//! A parent owns each child's *exterior* (what is visible from outside: mass, radius, orbit,
//! luminosity). The child owns its *interior* (the shell stack and how its lid behaves). The
//! surface pipeline sees only a [`BodyContext`].
//!
//! | Module | Concern |
//! |---|---|
//! | [`system`] | the orbital tree: stars, bodies, belts, rings |
//! | [`matter`] | species, compositions, phases |
//! | [`interior`] | the shell stack, tectonic regime and what they imply |
//! | [`context`] | the handoff to the surface pipeline |
//! | [`presets`] | real bodies (Earth, Mars, Venus, …) as starting points |

pub mod context;
pub mod interior;
pub mod matter;
pub mod presets;
pub mod system;

pub use context::{BodyContext, Irradiance, TidalPartner};
pub use interior::{BodyInterior, Shell, ShellKind, TectonicRegime};
pub use matter::{Composition, Phase, Species};
pub use presets::Archetype;
pub use system::{BodyClass, BodyExterior, Elements, Epoch, GalacticSite, NodeKind, OrbitalNode, Population, PowerLaw, Precession, SolarSystem, Star, StellarStage, SystemRoot};

pub const EARTH_MASS_KG: f64 = 5.972e24;
pub const EARTH_RADIUS_KM: f64 = 6371.0;
pub const SOLAR_MASS_KG: f64 = 1.989e30;
pub const SOLAR_CONSTANT_W_M2: f64 = 1361.0;
pub const AU_KM: f64 = 1.495_978_707e8;

/// Grid used for a body's surface products, recorded in the manifest.
pub type SurfaceGrid = crate::grid::GridSpec;
