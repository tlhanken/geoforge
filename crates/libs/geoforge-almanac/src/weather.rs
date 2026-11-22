//! Weather day by day at a place: drawn, seeded, from its climate, the season of the date
//! and the storms passing over it ([`crate::storms`]), so that neighbours share their
//! weather and over the years it adds up to the climate the layers hold.
//!
//! - **Season.** The star's declination as the ephemeris has it, as a share of its yearly
//!   swing, lagging a month behind over land and a month and a half over sea: 1 at the
//!   height of summer, −1 in winter (the hemisphere's own). The day's normal temperature and
//!   rain are the climate's month for that point of the season, as the Köppen view rebuilds
//!   them (the warm half-year's share of rain kept).
//! - **Temperature.** The normal plus an anomaly that carries from day to day (lag-one
//!   correlation 0.7 a day, as for daily means at mid-latitudes; Wilks & Wilby 1999), its
//!   spread growing with continentality: 1 K where the seasons are flat, 5 K in a
//!   continental interior. The daily range is 9 K in wet climates to 18 K in deserts (rainforest
//!   8–9 K, temperate 8–12 K, deserts 15–20 K), a quarter of that over open sea, 70% on
//!   rainy days, scaled by the square root of the day's length.
//! - **Rain.** Wet days follow a two-state Markov chain (Richardson 1981) whose wet share
//!   comes from the day's normal rain, r/(r + 4 mm) (0.4 for 3 mm a day; Geng et al. 1986
//!   give 0.3–0.5 for temperate climates), with a lag-one persistence of 0.35; a wet day's
//!   amount is exponential about r divided by that share. Storms raise both where they
//!   pass and lower them elsewhere, scaled by the place's own long-run storminess so the
//!   year's rain is unchanged. Below 0.5 °C it falls as snow.
//! - **Storm air.** A warm sector adds warmth and the air behind a cold front takes it away
//!   (up to about one day-to-day spread), the rest of the spread staying local and random.
//! - **Wind.** The prevailing wind between January's and July's by the season plus the
//!   storms' winds, varying day to day.
//! - **Phenomena.** Each hazard near the place (within 25 km) happens on a day with a
//!   chance from its yearly rate in `place`, weighted toward the days that favour it and
//!   normalised so the year's count is the rate: thunderstorms (thunder days from the
//!   lightning tier) and tornadoes on warm days under fronts and warm sectors, floods and
//!   landslides under heavy storm rain, blizzards on cold stormy days, wildfires on hot dry
//!   windy days, dust storms on dry windy days, fog on calm days; earthquakes, eruptions,
//!   tsunamis and auroras on any day. Tropical cyclones are the storm field's own.
//!
//! Draws are keyed by the day itself, and a month of spin-up runs before the window, so a
//! date's weather does not depend on where the window starts (to within a thousandth of a
//! degree, the anomaly's memory of the start).

use crate::ephemeris::{self, Calendar, Equatorial};
use crate::storms::{self, Atmosphere, Effect, StormKind};
use crate::{Error, Result};
use geoforge_analysis::hazards::Hazard;
use geoforge_analysis::koppen;
use geoforge_core::grid::Grid;
use geoforge_core::seed::{unit, Seed};
use geoforge_core::UnitVec3;
use geoforge_store::World;
use serde::Serialize;

/// Lag-one correlation of daily mean temperature anomalies, per Earth day.
const TEMPERATURE_MEMORY: f64 = 0.7;
/// Lag-one correlation of wet and dry days.
const RAIN_PERSISTENCE: f64 = 0.35;
/// Daily rain at which half of days are wet, mm.
const HALF_WET_MM: f64 = 4.0;
/// Steps simulated before the window.
const SPIN_UP: i64 = 30;

/// One day's weather.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Day {
    /// Earth days from the present, at the step's start.
    pub day: f64,
    /// The local date: year and day of the year.
    pub year: i64,
    pub local_day: f64,
    pub mean_c: f64,
    pub high_c: f64,
    pub low_c: f64,
    /// Rain or the water of snow, mm.
    pub precipitation_mm: f64,
    pub snow: bool,
    pub wind_m_s: f64,
    /// The bearing the wind blows from, degrees.
    pub wind_from_deg: f64,
    /// Sea-level pressure from the storms, on a 1013 hPa background.
    pub pressure_hpa: f64,
    pub sky: &'static str,
    /// The storm passing, if any.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub storm: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub storm_kind: Option<StormKind>,
    /// What else the day brings (thunderstorms, tornadoes, quakes, auroras …).
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub events: Vec<Event>,
}

/// A phenomenon on a day.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub struct Event {
    pub kind: &'static str,
    pub what: &'static str,
}

/// Something a day can bring, at a rate a year near the place.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Phenomenon {
    Thunderstorm,
    Hazard(Hazard),
}

impl Phenomenon {
    fn event(self) -> Event {
        match self {
            Self::Thunderstorm => Event { kind: "thunderstorm", what: "a thunderstorm" },
            Self::Hazard(h) => Event { kind: h.slug(), what: h.event() },
        }
    }

    /// How much a moment favours it, given the storms and the day's normal temperature.
    fn likelihood(self, e: &Effect, normal_c: f64) -> f64 {
        let warm = |base: f64| ((normal_c - base) / 10.0).max(0.0);
        let wind = e.wind[0].hypot(e.wind[1]);
        match self {
            Self::Thunderstorm => warm(8.0) * (0.25 + 2.5 * e.rain),
            Self::Hazard(h) => match h {
                Hazard::Tornadoes => warm(12.0) * (0.05 + 3.0 * e.rain * (1.0 + e.warmth.max(0.0))),
                Hazard::Floods | Hazard::Landslides => 0.05 + 3.0 * e.rain * e.rain,
                Hazard::Wildfire => warm(5.0) * (1.0 - e.rain).powi(2) * (1.0 + wind / 10.0),
                Hazard::Blizzards => ((2.0 - normal_c) / 5.0).max(0.0) * (0.05 + e.rain) * (1.0 + wind / 10.0),
                Hazard::DustStorms => (1.0 - e.rain) * (1.0 + wind / 5.0),
                Hazard::Fog => (1.0 - 0.7 * e.rain) / (1.0 + wind / 4.0),
                // The storm field's own; never drawn here.
                Hazard::Cyclones => 0.0,
                Hazard::Earthquakes | Hazard::Tsunamis | Hazard::Volcanoes | Hazard::Auroras => 1.0,
            },
        }
    }
}

/// The weather at a place over a run of days.
#[derive(Debug, Clone, Serialize)]
pub struct Forecast {
    pub lat: f64,
    pub lon: f64,
    pub at_sea: bool,
    /// Earth days per step: the solar day, or one Earth day on a world without days.
    pub step_days: f64,
    pub days: Vec<Day>,
}

/// A place's climate, as the weather is drawn from it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Climate {
    /// Annual mean, °C; warmest-to-coldest month, °C; annual rain, mm; and its share in the
    /// warm half-year.
    pub mean_c: f64,
    pub range_c: f64,
    pub precipitation_mm: f64,
    pub warm_share: f64,
    pub at_sea: bool,
    /// Prevailing winds in January and July (east, north), m/s.
    pub wind_jan: [f64; 2],
    pub wind_jul: [f64; 2],
    /// Latitude, degrees: which hemisphere's season applies.
    pub lat: f64,
}

/// `count` days of weather at (`lat`, `lon`) from `from_day` (Earth days from the present).
///
/// # Errors
/// If the world's climate is not generated yet.
pub fn forecast(world: &World, lat: f64, lon: f64, from_day: f64, count: usize) -> Result<Forecast> {
    let body = world.load_body()?;
    let grid = world.manifest.grid.build()?;
    let p = world.load_products()?;
    let missing = || Error::Invalid("the climate is not generated yet: run generate".into());
    let c = grid.locate(UnitVec3::from_lat_lon(lat, lon)) as usize;
    let wind = |f: Option<&geoforge_core::shapes::Field<[f32; 2]>>| f.map_or([0.0; 2], |f| f.values[c].map(f64::from));
    let at_sea = body.interior.has_surface_liquid_water() && f64::from(p.elevation().ok_or_else(missing)?.values[c]) < 0.0;
    let climate = Climate {
        mean_c: f64::from(p.temperature().ok_or_else(missing)?.values[c]),
        range_c: f64::from(p.temperature_range().ok_or_else(missing)?.values[c]),
        precipitation_mm: f64::from(p.precipitation().ok_or_else(missing)?.values[c]),
        warm_share: f64::from(p.precipitation_warm_share().ok_or_else(missing)?.values[c]),
        at_sea,
        wind_jan: wind(p.wind_jan()),
        wind_jul: wind(p.wind_jul()),
        lat,
    };
    let calendar = Calendar::new(&body, &[]);
    let year = calendar.year_days;
    let step_days = calendar.solar_day_h.map_or(1.0, |h| h / 24.0);
    let system = world.load_system()?;
    let atm = atmosphere(world, &p, &grid, &body)?;
    let near = crate::exposure::hazards_near(world, lat, lon, NEAR_KM)?;
    let rates: Vec<(Phenomenon, f64)> = std::iter::once((Phenomenon::Thunderstorm, near.thunder_days))
        .chain(near.hazards.iter().filter_map(|h| Hazard::from_slug(h.hazard).map(|k| (Phenomenon::Hazard(k), h.per_year))))
        .filter(|(p, r)| *r > 0.0 && *p != Phenomenon::Hazard(Hazard::Cyclones))
        .collect();
    let to_day = from_day + count as f64 * step_days;
    let passing = atm.storms(from_day - step_days * SPIN_UP as f64, to_day);
    let place = Storminess::of(&atm, lat, lon, year);
    let storms_at = |t: f64| storms::effect(&passing, lat, lon, t, atm.radius_km, atm.omega);
    let seed = world.manifest.root_seed.child(b"weather").child_u64(b"cell", c as u64);
    let mut days = simulate(&climate, season_of(system, world.manifest.focus, &body, year), year, step_days, from_day, count, seed, &storms_at, &place, &rates);
    for d in &mut days {
        (d.year, d.local_day) = calendar.date(d.day);
    }
    Ok(Forecast { lat, lon, at_sea, step_days, days })
}

/// Phenomena are counted within this distance of the place, km.
pub const NEAR_KM: f64 = 25.0;

/// The storm field of a world's climate.
///
/// # Errors
/// If the climate is not generated yet.
pub fn atmosphere(world: &World, p: &geoforge_core::product::Products, grid: &dyn geoforge_core::grid::Grid, body: &geoforge_core::BodyContext) -> Result<Atmosphere> {
    let year = Calendar::new(body, &[]).year_days;
    let season = season_of(world.load_system()?, world.manifest.focus, body, year);
    Atmosphere::read(p, grid, body, year, season, world.manifest.root_seed).ok_or_else(|| Error::Invalid("the climate is not generated yet: run generate".into()))
}

/// A storm on the weather map.
#[derive(Debug, Clone, Serialize)]
pub struct StormNow {
    pub name: String,
    pub kind: StormKind,
    pub lat: f64,
    pub lon: f64,
    pub depth_hpa: f64,
    /// Days since it formed.
    pub age_days: f64,
}

/// Draw the weather map of day `day` (`width` wide) to `to` ([`weather_map_image`]).
/// Returns the storms on the map.
///
/// # Errors
/// If the climate is not generated, or on IO errors.
pub fn weather_map(world: &World, day: f64, width: u32, to: &std::path::Path) -> Result<Vec<StormNow>> {
    let (img, storms) = weather_map_image(world, day, width)?;
    img.save(to).map_err(geoforge_store::StoreError::from)?;
    Ok(storms)
}

/// The weather map of day `day`, `width` wide: the season's sea-level pressure plus the
/// storms, isobars every 4 hPa, rain, fronts and storm centres; and the storms on it.
///
/// # Errors
/// If the climate is not generated.
pub fn weather_map_image(world: &World, day: f64, width: u32) -> Result<(image::RgbImage, Vec<StormNow>)> {
    use rayon::prelude::*;
    let body = world.load_body()?;
    let grid = world.manifest.grid.build()?;
    let p = world.load_products()?;
    let missing = || Error::Invalid("the climate is not generated yet: run generate".into());
    let year = Calendar::new(&body, &[]).year_days;
    let system = world.load_system()?;
    let season = season_of(system.clone(), world.manifest.focus, &body, year);
    let atm = atmosphere(world, &p, &grid, &body)?;
    let storms = atm.storms(day, day);
    // The season's pressure: January's and July's blended, a month behind the star.
    let pressure = |july| geoforge_analysis::pressure::seasonal(&p, &grid, &body, july);
    let (jan, jul) = (pressure(false), pressure(true));
    let mix = 0.5 * (1.0 + season(day - year / 12.0));
    let elevation = &p.elevation().ok_or_else(missing)?.values;
    let sea_world = body.interior.has_surface_liquid_water();
    let height = width / 2;
    let pixels: Vec<([u8; 3], f32, f32)> = (0..(width * height) as usize)
        .into_par_iter()
        .map(|i| {
            let (lat, lon) = (90.0 - ((i / width as usize) as f64 + 0.5) * 180.0 / f64::from(height), -180.0 + ((i % width as usize) as f64 + 0.5) * 360.0 / f64::from(width));
            let c = grid.locate(UnitVec3::from_lat_lon(lat, lon)) as usize;
            let e = f64::from(elevation[c]);
            let base = if sea_world && e < 0.0 { [24, 52, 96] } else { geoforge_store::render::elevation_color(e.max(0.0)).map(|v| (f64::from(v) * 0.75) as u8) };
            let seasonal = match (&jan, &jul) {
                (Some(a), Some(b)) => a[c] * (1.0 - mix) + b[c] * mix,
                _ => 1013.0,
            };
            let fx = storms::effect(&storms, lat, lon, day, atm.radius_km, atm.omega);
            (base, (seasonal + fx.pressure_hpa) as f32, fx.rain as f32)
        })
        .collect();
    let now: Vec<(&storms::Storm, storms::Fix)> = storms.iter().filter_map(|s| s.at(day).map(|f| (s, f))).filter(|(_, f)| f.depth_hpa > 2.0).collect();
    let fronts: Vec<(Vec<(f64, f64)>, bool)> = now.iter().filter_map(|(s, f)| s.fronts(*f, atm.radius_km)).flat_map(|[warm, cold]| [(warm, true), (cold, false)]).collect();
    let centres: Vec<(f64, f64, bool)> = now.iter().map(|(s, f)| (f.lat, f.lon, s.kind == StormKind::TropicalCyclone)).collect();
    let base: Vec<[u8; 3]> = pixels.iter().map(|p| p.0).collect();
    let pressure: Vec<f32> = pixels.iter().map(|p| p.1).collect();
    let rain: Vec<f32> = pixels.iter().map(|p| p.2).collect();
    let img = geoforge_store::render::render_weather(width, &base, &pressure, &rain, &fronts, &centres);
    Ok((img, now.iter().map(|(s, f)| StormNow { name: s.name.clone(), kind: s.kind, lat: f.lat, lon: f.lon, depth_hpa: f.depth_hpa, age_days: day - s.born_day }).collect()))
}

/// The northern season (−1 to 1) on Earth day `t`: the star's declination as a share of its
/// swing; a plain sinusoid from the equinox for a world without a system; 0 without tilt.
#[must_use]
pub fn season_of(system: Option<geoforge_core::body::SolarSystem>, focus: geoforge_core::Key, body: &geoforge_core::BodyContext, year: f64) -> Box<dyn Fn(f64) -> f64> {
    let obliquity = body.exterior.obliquity_deg.to_radians().sin().abs();
    let frame = system.and_then(|s| {
        let (f, star) = (Equatorial::of(&s, focus)?, ephemeris::light_source(&s, focus)?);
        Some((s, f, star))
    });
    Box::new(move |t: f64| match &frame {
        _ if obliquity < 0.02 => 0.0,
        Some((s, f, star)) => f.star_declination(s, focus, *star, t).to_radians().sin() / obliquity,
        None => (std::f64::consts::TAU * t / year).sin(),
    })
}

/// How stormy a place is over a year, which scales the storms' share of its weather.
#[derive(Debug, Clone, PartialEq)]
pub struct Storminess {
    /// Wetness factors through the year, and their mean.
    wet: Vec<f64>,
    mean_wet: f64,
    /// Mean warmth the storms bring.
    warmth: f64,
    /// The storms' effects through the year: (Earth day, effect).
    samples: Vec<(f64, Effect)>,
    /// Whether storms reach the place at all.
    stormy: bool,
}

/// Wetness factor of a storm effect: the odds of rain relative to a calm day.
fn wetness(e: &Effect) -> f64 {
    0.15 + 4.0 * e.rain
}

impl Storminess {
    /// No storms: weather from the place's climate alone.
    #[must_use]
    pub fn calm() -> Self {
        let w = wetness(&Effect::default());
        Self { wet: vec![w], mean_wet: w, warmth: 0.0, stormy: false, samples: (0..24).map(|k| (f64::from(k) * 365.25 / 24.0, Effect::default())).collect() }
    }

    /// The chance of a wet day for a day whose normal wet share is `share`, under wetness `g`.
    fn chance(&self, share: f64, g: f64) -> f64 {
        (share * g / self.mean_wet).min(0.95)
    }

    /// The factor on a wet day's amount, so that the year's rain is the normal's:
    /// `share / E[chance · √g]`.
    fn amount_factor(&self, share: f64) -> f64 {
        let expected = self.wet.iter().map(|&g| self.chance(share, g) * g.sqrt()).sum::<f64>() / self.wet.len() as f64;
        share / expected.max(1e-12)
    }

    /// The place's storminess over the first year from the present (the same for every
    /// window asked for).
    #[must_use]
    pub fn of(atm: &Atmosphere, lat: f64, lon: f64, year: f64) -> Self {
        let span = year.min(400.0);
        let reference = atm.storms(0.0, span);
        let n = 240;
        let times: Vec<f64> = (0..n).map(|k| span * f64::from(k) / f64::from(n)).collect();
        let effects: Vec<Effect> = times.iter().map(|&t| storms::effect(&reference, lat, lon, t, atm.radius_km, atm.omega)).collect();
        let wet: Vec<f64> = effects.iter().map(wetness).collect();
        let mean_wet = wet.iter().sum::<f64>() / f64::from(n);
        let warmth = effects.iter().map(|e| e.warmth).sum::<f64>() / f64::from(n);
        let stormy = effects.iter().any(|e| e.rain > 0.05 || e.warmth.abs() > 0.05);
        Self { wet, mean_wet, warmth, stormy, samples: times.into_iter().zip(effects).collect() }
    }
}

/// Draw `count` steps of `step_days` from `from_day` for `climate`, with `season(t)` the
/// northern season (−1 to 1) on Earth day `t`, `year` the year in Earth days, `storms(t)`
/// what the storms do here, `place` how stormy it is and `rates` the phenomena near it a year.
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn simulate(
    climate: &Climate,
    season: impl Fn(f64) -> f64,
    year: f64,
    step_days: f64,
    from_day: f64,
    count: usize,
    seed: Seed,
    storms: &dyn Fn(f64) -> Effect,
    place: &Storminess,
    rates: &[(Phenomenon, f64)],
) -> Vec<Day> {
    let months = koppen::months(climate.mean_c, climate.range_c, climate.precipitation_mm, climate.warm_share);
    let lag = year / 12.0 * if climate.at_sea { 1.5 } else { 1.0 };
    let hemisphere = (climate.lat / 10.0).clamp(-1.0, 1.0);
    // The month for a point of the season: 1 is month 0 (the warmest), −1 month 6.
    let normal = |t: f64| {
        let s = (season(t - lag) * hemisphere).clamp(-1.0, 1.0);
        let k = 6.0 * s.acos() / std::f64::consts::PI;
        let (i, f) = (k.floor() as usize, k.fract());
        let (a, b) = (months[i.min(6)], months[(i + 1).min(6)]);
        (a.t + f * (b.t - a.t), (a.rain + f * (b.rain - a.rain)) * 12.0 / year, s)
    };
    let spread = (0.8 + 0.12 * climate.range_c).clamp(0.8, 6.0);
    // Where storms pass they carry most of the day-to-day change.
    let local = if place.stormy { 0.6 } else { 1.0 };
    let memory = TEMPERATURE_MEMORY.powf(step_days);
    let day_range = (9.0 + 9.0 * (-climate.precipitation_mm / 500.0).exp()) * if climate.at_sea { 0.25 } else { 1.0 } * step_days.sqrt().clamp(0.4, 3.0);
    // Each phenomenon's chance on a step: its rate spread over the year, weighted by how much
    // the moment favours it relative to the year's mean.
    let mean_likelihood: Vec<f64> = rates.iter().map(|(p, _)| place.samples.iter().map(|(t, e)| p.likelihood(e, normal(*t).0)).sum::<f64>() / place.samples.len() as f64).collect();
    let first = (from_day / step_days).floor() as i64;
    let (mut anomaly, mut wet) = (0.0, false);
    let mut out = Vec::with_capacity(count);
    for n in first - SPIN_UP..first + count as i64 {
        let t = n as f64 * step_days;
        let draw = |label: &[u8]| unit(seed.child_u64(label, n as u64));
        let gauss = |label: &[u8]| {
            let (u1, u2) = (draw(label).max(1e-12), unit(seed.child_u64(label, n as u64).child(b"2")));
            (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
        };
        let (normal_c, rain_mm_day, s) = normal(t);
        let e = storms(t);
        anomaly = memory * anomaly + local * spread * (1.0 - memory * memory).sqrt() * gauss(b"temperature");
        let rain = rain_mm_day * step_days;
        let share = rain / (rain + HALF_WET_MM * step_days.max(1.0));
        // Storms shift the odds and the amounts, keeping the long-run mean.
        let g = wetness(&e);
        let share_now = place.chance(share, g);
        // Storms carry their own persistence; local showers keep a little.
        let persistence = if place.stormy { 0.15 } else { RAIN_PERSISTENCE };
        let p_wet = if wet { share_now + persistence * (1.0 - share_now) } else { share_now * (1.0 - persistence) };
        wet = draw(b"wet") < p_wet;
        let mean_amount = rain / share.max(1e-9) * g.sqrt() * place.amount_factor(share);
        let amount = if wet { -mean_amount * draw(b"amount").max(1e-12).ln() } else { 0.0 };
        if n < first {
            continue;
        }
        let mean_c = normal_c + anomaly + spread * (e.warmth - place.warmth);
        let half = 0.5 * day_range * if wet { 0.7 } else { 1.0 };
        let mix = 0.5 * (1.0 + s * hemisphere.signum());
        let blowing = [0, 1].map(|k| climate.wind_jan[k] * (1.0 - mix) + climate.wind_jul[k] * mix + e.wind[k]);
        let rayleigh = (-draw(b"wind").max(1e-12).ln()).sqrt() * 2.0 / std::f64::consts::PI.sqrt();
        // Steady storm winds vary less than light breezes.
        let steady = blowing[0].hypot(blowing[1]);
        let wind = steady * (0.75 + 0.25 * rayleigh) + 2.0 * rayleigh;
        let from = ((-blowing[0]).atan2(-blowing[1]).to_degrees() + 40.0 / (1.0 + steady / 5.0) * gauss(b"bearing")).rem_euclid(360.0);
        let snow = wet && mean_c < 0.5;
        let sky = match () {
            _ if e.kind == Some(StormKind::TropicalCyclone) && wind > 33.0 => "tropical cyclone",
            _ if wind > 17.0 && wet => "storm",
            _ if snow => "snow",
            _ if amount > 20.0 => "heavy rain",
            _ if wet => "rain",
            _ if draw(b"cloud") < share => "cloudy",
            _ => "clear",
        };
        let mut events = Vec::new();
        for (k, (&(p, per_year), &mean)) in rates.iter().zip(&mean_likelihood).enumerate() {
            // Rates are per year of Earth days, whatever the world's year.
            let chance = if mean > 0.0 { per_year * step_days / crate::exposure::EARTH_YEAR_DAYS * p.likelihood(&e, normal_c) / mean } else { 0.0 };
            if unit(seed.child_u64(b"event", k as u64).child_u64(b"step", n as u64)) < chance {
                events.push(p.event());
            }
        }
        if e.kind == Some(StormKind::TropicalCyclone) && e.rain > 0.3 {
            events.push(Phenomenon::Hazard(Hazard::Cyclones).event());
        }
        let thunder = events.iter().any(|ev| ev.kind == "thunderstorm");
        let sky = if thunder && wet && sky != "tropical cyclone" { "thunderstorm" } else { sky };
        let storm = e.storm.clone().filter(|_| e.rain > 0.3 || wind > 12.0);
        let storm_kind = storm.as_ref().and(e.kind);
        out.push(Day {
            day: t,
            year: 0,
            local_day: 0.0,
            mean_c,
            high_c: mean_c + half,
            low_c: mean_c - half,
            precipitation_mm: amount,
            snow,
            wind_m_s: wind,
            wind_from_deg: from,
            pressure_hpa: 1013.0 + e.pressure_hpa,
            sky,
            storm,
            storm_kind,
            events,
        });
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn temperate() -> Climate {
        Climate { mean_c: 10.0, range_c: 20.0, precipitation_mm: 800.0, warm_share: 0.6, at_sea: false, wind_jan: [6.0, 0.0], wind_jul: [3.0, 0.0], lat: 50.0 }
    }

    fn earth_season(t: f64) -> f64 {
        (std::f64::consts::TAU * t / 365.25).sin()
    }

    /// Ten years of weather add back up to the climate: its mean, its rain, its seasons.
    #[test]
    fn weather_adds_up_to_climate() {
        let c = temperate();
        let days = simulate(&c, earth_season, 365.25, 1.0, 0.0, 3653, Seed(9), &|_| Effect::default(), &Storminess::calm(), &[]);
        let n = days.len() as f64;
        let mean = days.iter().map(|d| d.mean_c).sum::<f64>() / n;
        let rain = days.iter().map(|d| d.precipitation_mm).sum::<f64>() / 10.0;
        let wet = days.iter().filter(|d| d.precipitation_mm > 0.0).count() as f64 / n;
        assert!((mean - 10.0).abs() < 0.6, "mean {mean}");
        assert!((rain / 800.0 - 1.0).abs() < 0.1, "rain {rain}");
        assert!((0.25..0.5).contains(&wet), "wet share {wet}");
        // July (days 100–130 from the March equinox) is about the range warmer than January.
        let avg = |lo: usize, hi: usize| (0..10).flat_map(|y| days[y * 365 + lo..y * 365 + hi].iter()).map(|d| d.mean_c).sum::<f64>() / (10 * (hi - lo)) as f64;
        let swing = avg(110, 140) - avg(290, 320);
        assert!((15.0..21.0).contains(&swing), "swing {swing}");
        // Winter brings snow, summer none; westerlies blow from the west.
        assert!(days.iter().any(|d| d.snow) && days[100..160].iter().all(|d| !d.snow));
        let west = days.iter().filter(|d| (180.0..360.0).contains(&d.wind_from_deg)).count() as f64 / n;
        assert!(west > 0.8, "{west}");
    }

    /// A date's weather is the same whatever the window, and the same seed repeats it.
    #[test]
    fn a_date_has_one_weather() {
        let c = temperate();
        let a = simulate(&c, earth_season, 365.25, 1.0, 200.0, 10, Seed(3), &|_| Effect::default(), &Storminess::calm(), &[]);
        let b = simulate(&c, earth_season, 365.25, 1.0, 195.0, 15, Seed(3), &|_| Effect::default(), &Storminess::calm(), &[]);
        for (x, y) in a.iter().zip(&b[5..]) {
            assert!((x.mean_c - y.mean_c).abs() < 1e-3 && x.precipitation_mm == y.precipitation_mm && x.sky == y.sky, "{x:?} {y:?}");
        }
        assert_ne!(a, simulate(&c, earth_season, 365.25, 1.0, 200.0, 10, Seed(4), &|_| Effect::default(), &Storminess::calm(), &[]));
    }

    /// With storms passing, the years still add up to the climate, and two nearby places
    /// get rain on the same days while distant ones do not.
    #[test]
    fn storms_keep_the_climate_and_are_shared() {
        let atm = crate::storms::tests::earthlike(4);
        let c = temperate();
        let storms = atm.storms(-20.0, 3653.0);
        let run = |lat: f64, lon: f64, seed: u64| {
            let place = Storminess::of(&atm, lat, lon, 365.25);
            let at = |t: f64| storms::effect(&storms, lat, lon, t, atm.radius_km, atm.omega);
            simulate(&Climate { lat, ..c }, earth_season, 365.25, 1.0, 0.0, 3653, Seed(seed), &at, &place, &[])
        };
        let (a, near, far) = (run(50.0, -30.0, 1), run(50.5, -29.0, 2), run(50.0, 60.0, 3));
        let rain = a.iter().map(|d| d.precipitation_mm).sum::<f64>() / 10.0;
        let mean = a.iter().map(|d| d.mean_c).sum::<f64>() / a.len() as f64;
        assert!((rain / 800.0 - 1.0).abs() < 0.15, "rain {rain}");
        assert!((mean - 10.0).abs() < 0.8, "mean {mean}");
        let wet = |d: &[Day]| d.iter().map(|x| f64::from(u8::from(x.precipitation_mm > 0.0))).collect::<Vec<_>>();
        let (near_r, far_r) = (crate::storms::tests::correlation(&wet(&a), &wet(&near)), crate::storms::tests::correlation(&wet(&a), &wet(&far)));
        assert!(near_r > 0.3 && far_r < 0.1, "near {near_r}, far {far_r}");
        assert!(a.iter().any(|d| d.storm.is_some()));
    }

    /// Phenomena come at their yearly rates, on the days that favour them: thunder and
    /// tornadoes in the warm season and under storms, blizzards in winter.
    #[test]
    fn phenomena_keep_their_rates() {
        let atm = crate::storms::tests::earthlike(11);
        let c = Climate { lat: 45.0, ..temperate() };
        let storms = atm.storms(-20.0, 3653.0);
        let place = Storminess::of(&atm, 45.0, -30.0, 365.25);
        let at = |t: f64| storms::effect(&storms, 45.0, -30.0, t, atm.radius_km, atm.omega);
        let rates = [(Phenomenon::Thunderstorm, 25.0), (Phenomenon::Hazard(Hazard::Tornadoes), 2.0), (Phenomenon::Hazard(Hazard::Blizzards), 4.0), (Phenomenon::Hazard(Hazard::Earthquakes), 1.0)];
        let days = simulate(&c, earth_season, 365.25, 1.0, 0.0, 3653, Seed(2), &at, &place, &rates);
        let count = |kind: &str, d: &[Day]| d.iter().filter(|x| x.events.iter().any(|e| e.kind == kind)).count() as f64;
        for (kind, rate) in [("thunderstorm", 25.0), ("tornadoes", 2.0), ("blizzards", 4.0), ("earthquakes", 1.0)] {
            let per_year = count(kind, &days) / 10.0;
            assert!((per_year / rate - 1.0).abs() < 0.4, "{kind}: {per_year} a year against {rate}");
        }
        // Summer (days 60–150 from the March equinox) has most of the thunder, winter the blizzards.
        let summer: Vec<Day> = days.iter().filter(|d| (60.0..150.0).contains(&d.day.rem_euclid(365.25))).cloned().collect();
        assert!(count("thunderstorm", &summer) > 0.5 * count("thunderstorm", &days));
        assert_eq!(count("blizzards", &summer), 0.0);
        // Thunder comes with storms more than without.
        let stormy = days.iter().filter(|d| d.storm.is_some()).count() as f64;
        let thunder_stormy = days.iter().filter(|d| d.storm.is_some() && d.events.iter().any(|e| e.kind == "thunderstorm")).count() as f64;
        assert!(thunder_stormy / stormy > count("thunderstorm", &days) / days.len() as f64);
    }
}
