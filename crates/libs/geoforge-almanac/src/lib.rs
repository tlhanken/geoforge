//! A world in time: where its bodies are on any day, its calendar, eclipses and night sky,
//! the storms crossing its map, the weather day by day at a place, and the hazards near it.
//! Everything here is derived from the stored world and seeded by date, so any day costs
//! the same and comes out the same however it is asked for.
//!
//! | Module | Concern |
//! |---|---|
//! | [`ephemeris`] | Kepler positions, calendars, the equatorial frame and seasons, eclipses |
//! | [`sky`] | the night sky on a date: neighbourhood stars, the galaxy's glow, planets, moons |
//! | [`storms`] | the global storm field: lows with fronts, tropical cyclones, steered by winds and mountains |
//! | [`weather`] | daily weather at a place, its phenomena, and the weather map |
//! | [`exposure`] | hazards near a place: yearly rates, chances over a window, quakes |

pub mod drought;
pub mod ephemeris;
pub mod exposure;
pub mod sky;
pub mod storms;
pub mod weather;
pub mod year;

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("{0}")]
    Store(#[from] geoforge_store::StoreError),
    #[error("{0}")]
    Core(#[from] geoforge_core::Error),
    #[error("{0}")]
    Invalid(String),
}

pub type Result<T> = std::result::Result<T, Error>;
