//! A place's almanac for one of its years: each local day's weather, sunrise and sunset, how
//! lit the moons are and the day's events, with the year's months of rain, drought index and
//! soil water.

use crate::drought::{self, Droughts, Month};
use crate::ephemeris::{self, Calendar};
use crate::weather::{self, Day};
use crate::Result;
use geoforge_store::World;
use serde::Serialize;

/// Years of record the drought index of an almanac is read against.
pub const DROUGHT_YEARS: usize = 30;

/// One local day of an almanac.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct AlmanacDay {
    pub weather: Day,
    /// Sunrise and sunset in local hours from midnight; `None` in polar night (both 0 and
    /// the day's length under the midnight sun), and on a world without days.
    pub sun: Option<(f64, f64)>,
    /// Each moon's lit share at the day's start, as seen from the world.
    pub moons: Vec<(String, f64)>,
}

/// A place's year, day by day.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Almanac {
    pub lat: f64,
    pub lon: f64,
    /// The local year (0 the present one).
    pub year: i64,
    /// Local day length, hours (`None` for a world locked to its star).
    pub day_h: Option<f64>,
    pub days: Vec<AlmanacDay>,
    /// The year's whole months, from the days above: rain, SPI-3 read against the record
    /// before the year, and soil water.
    pub months: Vec<Month>,
    /// The place's droughts over the [`DROUGHT_YEARS`] before the year.
    pub droughts: Droughts,
}

/// The almanac of (`lat`, `lon`) for local year `year` (0 the present one).
///
/// # Errors
/// If the climate is not generated yet.
pub fn year(world: &World, lat: f64, lon: f64, year: i64) -> Result<Almanac> {
    let body = world.load_body()?;
    let calendar = Calendar::new(&body, &[]);
    let year_days = calendar.year_days;
    let step = calendar.solar_day_h.map_or(1.0, |h| h / 24.0);
    let start = year as f64 * year_days;
    let count = (year_days / step).round() as usize;
    // Three months before the year too, for the first months' three-month drought index,
    // and two days after, so its last month is whole.
    let lead = (3.0 * year_days / 12.0 / step).ceil() as usize;
    let tail = (2.0 / step).ceil() as usize;
    let forecast = weather::forecast(world, lat, lon, start - lead as f64 * step, lead + count + tail)?;
    let droughts = drought::droughts(world, lat, lon, start - 1.0, DROUGHT_YEARS)?;
    let drawn: Vec<(f64, f64, f64)> = forecast.days.iter().map(|d| (d.day + 0.5 * step, d.precipitation_mm, d.mean_c)).collect();
    let months = droughts.read(&drawn).into_iter().filter(|m| m.day >= start && m.day < start + year_days).collect();
    let system = world.load_system()?;
    let season = weather::season_of(system.clone(), world.manifest.focus, &body, year_days);
    let obliquity = body.irradiance.obliquity_deg;
    let days = forecast
        .days
        .into_iter()
        .filter(|w| w.year == year)
        .map(|w| {
            let t = w.day;
            let sun = calendar.solar_day_h.and_then(|h| drought::sun_times(lat, obliquity * season(t + 0.5 * step), h));
            let moons = system.as_ref().map_or_else(Vec::new, |sys| {
                let placed = ephemeris::positions(sys, t);
                let at = |k| placed.iter().find(|p| p.key == k).map(|p| p.position_au);
                let Some(here) = at(world.manifest.focus) else { return Vec::new() };
                sys.children(world.manifest.focus).iter().filter_map(|m| at(m.key).map(|pos| (m.name.clone(), ephemeris::lit_fraction(pos, here, [0.0; 3])))).collect()
            });
            AlmanacDay { weather: w, sun, moons }
        })
        .collect();
    Ok(Almanac { lat, lon, year, day_h: calendar.solar_day_h, days, months, droughts })
}
