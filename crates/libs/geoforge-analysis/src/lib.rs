//! What a generated surface means. The surface layers ([`geoforge_surface`]) make the
//! world; this crate reads it: the climate as Köppen classes, the ground as soils, timber
//! and ore deposits (with the giant-deposit catalog, the one stored product made here),
//! natural hazards, the oceans' ice, reefs and productivity, and sea-level pressure.
//!
//! | Module | Concern |
//! |---|---|
//! | [`appearance`] | true colour: a generated surface, a body without one, a star |
//! | [`ground`] | one read of everything the analyses need, per cell |
//! | [`koppen`] | Köppen–Geiger classes from the climate |
//! | [`resources`] | ore deposits (λ prospectivity, local draws, giants), soils, timber |
//! | [`hazards`] | earthquakes, tsunamis, eruptions, cyclones, tornadoes, floods, fire, landslides, blizzards, dust, fog, auroras, lightning |
//! | [`oceans`] | sea ice, coral reefs and atolls, productivity |
//! | [`pressure`] | sea-level pressure from the winds |

pub mod appearance;
pub mod ground;
pub mod hazards;
pub mod koppen;
pub mod oceans;
pub mod pressure;
pub mod resources;

pub use resources::GiantDeposits;
