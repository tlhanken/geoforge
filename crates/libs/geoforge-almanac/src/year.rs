//! A place's almanac for one of its years: each local day's date in a calendar, weather,
//! sunrise and sunset, how lit the moons are and the day's events, with the year's months of
//! rain, drought index and soil water. The year is the calendar's: the whole days that begin
//! in it. The drought months are twelfths of the year whatever the calendar, as the seasons
//! follow the orbit.

use crate::drought::{self, Droughts, Month};
use crate::ephemeris;
use crate::weather::{self, Day};
use crate::Result;
use geoforge_core::time::{Calendar, Date, Instant, Month as CalendarMonth};
use geoforge_store::World;
use serde::Serialize;

/// Years of record the drought index of an almanac is read against.
pub const DROUGHT_YEARS: usize = 20;

/// One local day of an almanac.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct AlmanacDay {
    pub date: Date,
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
    /// The local year (0 the present one), as the calendar writes it.
    pub year: i64,
    pub year_name: String,
    /// The calendar the dates are in, and its months this year.
    pub calendar: String,
    pub calendar_months: Vec<CalendarMonth>,
    /// Local day length, hours (`None` for a world locked to its star).
    pub day_h: Option<f64>,
    pub days: Vec<AlmanacDay>,
    /// The year's whole twelfths, from the days above: rain, SPI-3 read against the record
    /// before the year, and soil water.
    pub months: Vec<Month>,
    /// The place's droughts over the [`DROUGHT_YEARS`] before the year.
    pub droughts: Droughts,
}

/// The almanac of (`lat`, `lon`) for year `year` (0 the present one) of `calendar`.
///
/// # Errors
/// If the climate is not generated yet.
pub fn year(world: &World, lat: f64, lon: f64, year: i64, calendar: &dyn Calendar) -> Result<Almanac> {
    let body = world.load_body()?;
    let clock = *calendar.clock();
    let year_days = clock.year.earth_days();
    let step = clock.day_len().earth_days();
    let (first, last) = (clock.first_day_from(clock.year_start(year)), clock.first_day_from(clock.year_start(year + 1)));
    let start = clock.day_start(first).earth_days();
    let count = (last - first) as usize;
    // Three months before the year too, for the first months' three-month drought index,
    // and two days after, so its last month is whole.
    let lead = (3.0 * year_days / 12.0 / step).ceil() as usize;
    let tail = (2.0 / step).ceil() as usize;
    let forecast = weather::forecast(world, lat, lon, start - lead as f64 * step, lead + count + tail)?;
    let droughts = drought::droughts(world, lat, lon, start - 1.0, DROUGHT_YEARS)?;
    let drawn: Vec<drought::Drawn> = forecast.days.iter().map(|d| drought::Drawn::of(d, step)).collect();
    let (from, to) = (clock.year_start(year).earth_days(), clock.year_start(year + 1).earth_days());
    let months = droughts.read(&drawn).into_iter().filter(|m| m.day >= from && m.day < to).collect();
    let system = world.load_system()?;
    let season = weather::season_of(system.clone(), world.manifest.focus, &body, year_days);
    let obliquity = body.irradiance.obliquity_deg;
    let days = forecast
        .days
        .into_iter()
        .filter(|w| (first..last).contains(&clock.day_count(Instant::from_earth_days(w.day + 0.5 * step)).0))
        .map(|w| {
            let t = w.day;
            let sun = clock.day.and_then(|d| drought::sun_times(lat, obliquity * season(t + 0.5 * step), d.hours()));
            let moons = system.as_ref().map_or_else(Vec::new, |sys| {
                let placed = ephemeris::positions(sys, t);
                let at = |k| placed.iter().find(|p| p.key == k).map(|p| p.position_au);
                let Some(here) = at(world.manifest.focus) else { return Vec::new() };
                sys.children(world.manifest.focus).iter().filter_map(|m| at(m.key).map(|pos| (m.name.clone(), ephemeris::lit_fraction(pos, here, [0.0; 3])))).collect()
            });
            let date = calendar.date(Instant::from_earth_days(t + 0.5 * step));
            AlmanacDay { date, weather: w, sun, moons }
        })
        .collect();
    Ok(Almanac { lat, lon, year, year_name: calendar.year_name(year), calendar: calendar.name(), calendar_months: calendar.months(year), day_h: clock.day.map(|d| d.hours()), days, months, droughts })
}
