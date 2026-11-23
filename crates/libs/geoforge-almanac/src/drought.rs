//! Wet and dry years at a place: a slow rain anomaly, the drought index read from it, and the
//! soil's water.
//!
//! - **The slow anomaly.** The rain of any stretch of days is the climate's normal times a
//!   slow factor, whose logarithm is a month-to-month AR(1) (a dry spell keeps 60% of itself
//!   from one month to the next, so it outlasts a season) plus an ENSO-like cycle. The AR(1)
//!   draws are read off a 15° lattice, so neighbours share their wet and dry years. The cycle
//!   is one world-wide oscillation of two seeded periods, felt most in the tropics. Along the
//!   basin that sets it its warm phase wets the eastern end and its coast and dries the
//!   western end and the lands beyond it, and the land beyond the eastern coast, as El Niño
//!   wets coastal Ecuador and dries Indonesia, northern Australia and the Amazon when the
//!   Walker circulation's rising branch moves east and its sinking branches flank it
//!   (Ropelewski & Halpert 1987); elsewhere the sign of a region's response is a broad
//!   seeded pattern (the far teleconnections). Drier climates vary more: the logarithm's spread runs from 0.2 in wet
//!   climates to 0.5 in deserts.
//! - **The cycle's basin.** ENSO is a delayed oscillator (Suarez & Schopf 1988; Battisti &
//!   Hirst 1989): an equatorial Kelvin wave crosses the widest tropical basin eastward, and
//!   Rossby waves, a third as fast, carry the signal back to be reflected off its western
//!   shore, so the period grows with the basin's width (Kelvin waves cross the Pacific in
//!   about two months). The periods are Earth's (2.5–5 and 4–7 years, drawn) times the width
//!   of the world's widest equatorial sea over the Pacific's as measured the same way
//!   ([`EARTH_BASIN_KM`]), and the cycle weakens below 12,000 km and is gone under 4,000 km
//!   (the Atlantic's 5,000 km has only a weak Atlantic Niño). A world with no shore on the
//!   equator, or no sea, has no cycle. Earth-like thermoclines are assumed, so the wave
//!   speeds are Earth's.
//! - **SPI-3**, the standardized precipitation index over three months (McKee et al. 1993):
//!   each three-month total set against the same months in the other years, here by a
//!   log-normal fit; −1 and below is a moderate drought, −1.5 a severe one. Months whose
//!   three-month totals average under 10 mm are too dry to read. The record is the place's
//!   own daily weather over the years, so the index is set on the same storms and showers
//!   as the days it reads.
//! - **Soil water**: a 150 mm bucket filled by rain and emptied by evaporation, potential
//!   evaporation rising with temperature (5 mm a day at 25 °C, none below −3 °C) and cut as
//!   the soil dries.

use crate::ephemeris::Calendar;
use crate::weather::Climate;
use crate::Result;
use geoforge_core::seed::{unit, Seed};
use geoforge_store::World;
use serde::Serialize;
use std::f64::consts::{PI, TAU};

/// Share of a month's slow anomaly kept into the next.
const PERSISTENCE: f64 = 0.6;
/// Months of draws the AR(1) is summed over (0.6^24 ≈ 5e-6).
const MEMORY_MONTHS: i64 = 24;
/// Spacing of the lattice the regional draws are read off, degrees.
const REGION_DEG: f64 = 15.0;
/// Spacing of the lattice of the ENSO-like pattern, degrees.
const PATTERN_DEG: f64 = 45.0;
/// The largest share of a tropical place's slow anomaly that is the ENSO-like cycle.
const CYCLE_SHARE: f64 = 0.6;
/// The Pacific's width along the equator as [`equatorial_basin`] measures it on real Earth
/// (the 0.5° Earth world: 118°E to 80°W), km.
pub const EARTH_BASIN_KM: f64 = 18_000.0;
/// Basin widths over which the cycle grows from nothing to full, km.
const CYCLE_BASIN_KM: (f64, f64) = (4_000.0, 12_000.0);
/// How far the basin's own response reaches past its ends, degrees of longitude: the dry
/// lands beyond its western end (Indonesia, northern Australia), its wet eastern coast,
/// and the dry land beyond that (the Amazon); the half-width over which it fades beyond
/// those, and the latitudes it fills before fading over [`FADE_DEG`]'s half.
const WEST_DRY_DEG: f64 = 25.0;
const EAST_COAST_DEG: f64 = 5.0;
const EAST_DRY_DEG: f64 = 40.0;
const FADE_DEG: f64 = 20.0;
const RESPONSE_LAT_DEG: f64 = 15.0;

/// The widest stretch of sea along a world's equator.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct Basin {
    pub width_km: f64,
    /// Its western and eastern shores, degrees east.
    pub west_lon: f64,
    pub east_lon: f64,
}

impl Basin {
    /// The sign of the warm phase's rain at (`lat`, `lon`) from the basin's own shift (−1 at
    /// its western end to 1 at its eastern, −1 on the lands beyond either end but the eastern
    /// coast), and how much that sets the place's response (1 near the equator on the basin
    /// and the lands it reaches, fading away from them).
    #[must_use]
    pub fn response(&self, lat: f64, lon: f64) -> (f64, f64) {
        let span = (self.east_lon - self.west_lon).rem_euclid(360.0);
        let mut along = (lon - self.west_lon).rem_euclid(360.0);
        // Past the eastern end up to the far side's midpoint, before the western end after.
        if along > span + 0.5 * (360.0 - span) {
            along -= 360.0;
        }
        let fade = |d: f64, reach: f64| if d <= reach { 1.0 } else { (-((d - reach) / FADE_DEG).powi(2)).exp() };
        let (sign, near) = match along {
            a if a < 0.0 => (-1.0, fade(-a, WEST_DRY_DEG)),
            a if a <= span => (-(PI * a / span.max(1e-9)).cos(), 1.0),
            a if a - span <= EAST_COAST_DEG => (1.0, 1.0),
            a => (-1.0, fade(a - span, EAST_DRY_DEG)),
        };
        (sign, near * (-((lat.abs() - RESPONSE_LAT_DEG).max(0.0) / (0.5 * FADE_DEG)).powi(2)).exp())
    }
}

/// The widest basin along the equator of a world of radius `radius_km` where `sea(lat, lon)`:
/// a longitude (every quarter degree) is sea where most of 4°S, 2°S, 0°, 2°N and 4°N are, so
/// small islands do not split a basin. `None` if the equator is all sea or all land.
#[must_use]
pub fn equatorial_basin(sea: impl Fn(f64, f64) -> bool, radius_km: f64) -> Option<Basin> {
    const STEP: f64 = 0.25;
    let n = (360.0 / STEP) as usize;
    let open: Vec<bool> = (0..n).map(|i| [-4.0, -2.0, 0.0, 2.0, 4.0].iter().filter(|&&la| sea(la, -180.0 + STEP * i as f64)).count() >= 3).collect();
    let shore = open.iter().position(|&o| !o)?;
    // Runs of sea, counted from a shore so none wraps around.
    let (mut best, mut start, mut run) = ((0, 0), 0, 0);
    for k in 1..=n {
        if open[(shore + k) % n] {
            start = if run == 0 { k } else { start };
            run += 1;
            if run > best.1 {
                best = (start, run);
            }
        } else {
            run = 0;
        }
    }
    (best.1 > 0).then(|| {
        let (west, span) = (-180.0 + STEP * ((shore + best.0) % n) as f64, STEP * best.1 as f64);
        Basin { width_km: radius_km * span.to_radians(), west_lon: west, east_lon: (west + span + 180.0).rem_euclid(360.0) - 180.0 }
    })
}
/// Soil water capacity, mm.
pub(crate) const SOIL_MM: f64 = 150.0;
/// Snowmelt per degree above freezing, mm of water a day (degree-day factors for snow are
/// 2.5–6; Hock 2003).
const MELT_MM_PER_C_DAY: f64 = 3.0;

/// The soil's water (mm of [`SOIL_MM`]) and the snowpack (mm of water) after `days` days in
/// which `water_mm` falls (as snow if `snow`) at a mean of `mean_c` °C: snow melts at
/// [`MELT_MM_PER_C_DAY`] a degree, and the soil loses 0.18 (T + 3) mm a day when full, less
/// as it dries (a bucket; Thornthwaite-like in its temperature).
#[must_use]
pub(crate) fn bucket(soil: f64, pack: f64, water_mm: f64, snow: bool, mean_c: f64, days: f64) -> (f64, f64) {
    let pack = pack + if snow { water_mm } else { 0.0 };
    let melt = (MELT_MM_PER_C_DAY * mean_c.max(0.0) * days).min(pack);
    let dry = days * (0.18 * (mean_c + 3.0)).max(0.0) * soil / SOIL_MM;
    ((soil + if snow { 0.0 } else { water_mm } + melt - dry).clamp(0.0, SOIL_MM), pack - melt)
}

/// One drawn day as the drought record reads it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Drawn {
    /// Earth days from the present, at the day's middle.
    pub day: f64,
    pub rain_mm: f64,
    pub mean_c: f64,
    /// The soil's water at the day's end, a share of what it holds, and the snow lying, mm of
    /// water.
    pub soil: f64,
    pub snowpack_mm: f64,
}

impl Drawn {
    /// A day of weather `step_days` long.
    #[must_use]
    pub fn of(d: &crate::weather::Day, step_days: f64) -> Self {
        Self { day: d.day + 0.5 * step_days, rain_mm: d.precipitation_mm, mean_c: d.mean_c, soil: d.soil, snowpack_mm: d.snowpack_mm }
    }
}
/// SPI at and below which a month is in moderate and in severe drought.
pub const MODERATE_SPI: f64 = -1.0;
pub const SEVERE_SPI: f64 = -1.5;

/// The slow wetness of a place: a factor on its normal rain.
#[derive(Debug, Clone)]
pub struct Wetness {
    seed: Seed,
    lat: f64,
    lon: f64,
    /// Spread of the factor's logarithm.
    spread: f64,
    /// The ENSO-like cycle's share here, signed.
    cycle: f64,
    periods: [f64; 2],
    phases: [f64; 2],
    month_days: f64,
}

impl Wetness {
    /// The wetness of (`lat`, `lon`) on a world seeded `seed` whose year is `year_days` Earth
    /// days long and whose widest equatorial basin is `basin`, for a place with `annual_mm`
    /// of rain.
    #[must_use]
    pub fn new(seed: Seed, lat: f64, lon: f64, annual_mm: f64, year_days: f64, basin: Option<Basin>) -> Self {
        let draw = |label: &[u8]| unit(seed.child(b"cycle").child(label));
        let width = basin.map_or(0.0, |b| b.width_km);
        let earth_year = crate::exposure::EARTH_YEAR_DAYS * basin.map_or(1.0, |b| b.width_km / EARTH_BASIN_KM);
        let periods = [2.5 + 2.5 * draw(b"p1"), 4.0 + 3.0 * draw(b"p2")].map(|y| y * earth_year);
        let (along, weight) = basin.map_or((0.0, 0.0), |b| b.response(lat, lon));
        let pattern = weight * along + (1.0 - weight) * lattice(seed.child(b"pattern"), lat, lon, PATTERN_DEG);
        let (lo, hi) = CYCLE_BASIN_KM;
        let x = ((width - lo) / (hi - lo)).clamp(0.0, 1.0);
        let cycle = CYCLE_SHARE * x * x * (3.0 - 2.0 * x) * (-(lat / 30.0).powi(2)).exp() * (1.5 * pattern).tanh();
        Self { seed, lat, lon, spread: 0.2 + 0.3 * (-annual_mm.max(0.0) / 500.0).exp(), cycle, periods, phases: [TAU * draw(b"f1"), TAU * draw(b"f2")], month_days: year_days / 12.0 }
    }

    /// The wetness of a place on `world`, whose climate is `climate`, year `year_days` and
    /// widest equatorial basin `basin`.
    #[must_use]
    pub fn of(world: &World, lat: f64, lon: f64, climate: &Climate, year_days: f64, basin: Option<Basin>) -> Self {
        Self::new(world.manifest.root_seed.child(b"drought"), lat, lon, climate.precipitation_mm, year_days, basin)
    }

    /// The cycle's two periods, Earth years.
    #[must_use]
    pub fn periods_years(&self) -> [f64; 2] {
        self.periods.map(|p| p / crate::exposure::EARTH_YEAR_DAYS)
    }

    /// The share of the place's slow anomaly that is the cycle, signed: positive where its
    /// warm phase wets, negative where it dries.
    #[must_use]
    pub fn cycle_share(&self) -> f64 {
        self.cycle
    }

    /// The world-wide cycle on Earth day `t`, of unit variance.
    #[must_use]
    pub fn cycle(&self, t: f64) -> f64 {
        let [a, b] = [0, 1].map(|k| (TAU * t / self.periods[k] + self.phases[k]).sin());
        (a + 0.6 * b) / (0.5 * (1.0 + 0.36f64)).sqrt()
    }

    /// The slow anomaly on Earth day `t`, of about unit variance: read between months so it
    /// has no steps.
    #[must_use]
    pub fn index(&self, t: f64) -> f64 {
        let m = t / self.month_days - 0.5;
        let k = m.floor() as i64;
        let f = m - k as f64;
        let local = (1.0 - f) * self.ar(k) + f * self.ar(k + 1);
        (1.0 - self.cycle * self.cycle).sqrt() * local + self.cycle * self.cycle(t)
    }

    /// The factor on the normal rain on Earth day `t`; 1 over the long run.
    #[must_use]
    pub fn factor(&self, t: f64) -> f64 {
        (self.spread * self.index(t) - 0.5 * self.spread * self.spread).exp()
    }

    /// Month `k`'s AR(1) value, summed from its draws.
    fn ar(&self, k: i64) -> f64 {
        let scale = (1.0 - PERSISTENCE * PERSISTENCE).sqrt();
        (0..=MEMORY_MONTHS).map(|j| PERSISTENCE.powi(j as i32) * scale * lattice(self.seed.child_u64(b"month", (k - j) as u64), self.lat, self.lon, REGION_DEG)).sum()
    }
}

/// A unit-variance Gaussian field read off a lattice of `deg` spacing (bilinear, rescaled so
/// every point has unit variance), seeded by `seed`.
fn lattice(seed: Seed, lat: f64, lon: f64, deg: f64) -> f64 {
    let rows = (180.0 / deg).round() as i64;
    let cols = (360.0 / deg).round() as i64;
    let fy = ((lat + 90.0) / deg).clamp(0.0, rows as f64 - 1e-9);
    let fx = (lon + 180.0).rem_euclid(360.0) / deg;
    let (y0, x0) = (fy.floor() as i64, fx.floor() as i64);
    let (ty, tx) = (fy - y0 as f64, fx - x0 as f64);
    let node = |y: i64, x: i64| gauss(seed.child_u64(b"node", (y.min(rows) * cols + x.rem_euclid(cols)) as u64));
    let w = [(1.0 - ty) * (1.0 - tx), (1.0 - ty) * tx, ty * (1.0 - tx), ty * tx];
    let v = w[0] * node(y0, x0) + w[1] * node(y0, x0 + 1) + w[2] * node(y0 + 1, x0) + w[3] * node(y0 + 1, x0 + 1);
    v / w.iter().map(|w| w * w).sum::<f64>().sqrt()
}

fn gauss(seed: Seed) -> f64 {
    let (u1, u2) = (unit(seed.child(b"1")).max(1e-12), unit(seed.child(b"2")));
    (-2.0 * u1.ln()).sqrt() * (TAU * u2).cos()
}

/// One month of a place's simulated record.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct Month {
    /// Earth days from the present, at the month's middle.
    pub day: f64,
    pub rain_mm: f64,
    /// SPI-3 ending with this month; `None` where it is too dry to read or before the third.
    pub spi3: Option<f64>,
    /// Soil water at the month's end, as a share of what the soil holds.
    pub soil: f64,
}

/// A place's droughts over a run of years.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Droughts {
    pub years: usize,
    /// Share of the readable months in moderate (SPI-3 ≤ −1) and in severe (≤ −1.5) drought.
    pub moderate_share: f64,
    pub severe_share: f64,
    /// The longest run of months in moderate drought.
    pub longest_months: usize,
    /// Share of growing-season months (mean ≥ 5 °C) whose soil is under a quarter full;
    /// `None` without a growing season.
    pub dry_soil_share: Option<f64>,
    /// The record, month by month, ending with the month of `now`.
    pub months: Vec<Month>,
    /// The soil's water (a share) and the snow lying (mm of water) on the record's last day.
    pub soil_now: f64,
    pub snowpack_now_mm: f64,
    /// The world's ENSO-like cycle as the place feels it.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub enso: Option<crate::weather::Enso>,
    /// The months' index (Earth day ÷ month length, floored) of the record's first month.
    #[serde(skip)]
    first: i64,
    /// Each calendar month's fit of its three-month totals: mean and spread of ln(mm + 1).
    #[serde(skip)]
    fits: Vec<Option<(f64, f64)>>,
    #[serde(skip)]
    month_days: f64,
}

impl Droughts {
    /// The whole months of a run of days, in order: (the months' index, rain mm, mean °C,
    /// the soil's water at the month's end).
    fn totals(&self, days: &[Drawn]) -> Vec<(i64, f64, f64, f64)> {
        let index = |t: f64| (t / self.month_days).floor() as i64;
        let (Some(a), Some(b)) = (days.first(), days.last()) else { return Vec::new() };
        let (from, to) = (index(a.day - 1e-9) + 1, index(b.day + 1e-9) - 1);
        let mut sums: Vec<(f64, f64, usize, f64)> = vec![(0.0, 0.0, 0, 0.0); (to - from + 1).max(0) as usize];
        for d in days {
            let k = index(d.day) - from;
            if (0..sums.len() as i64).contains(&k) {
                let m = &mut sums[k as usize];
                *m = (m.0 + d.rain_mm, m.1 + d.mean_c, m.2 + 1, d.soil);
            }
        }
        sums.iter().zip(from..).map(|(m, k)| (k, m.0, m.1 / m.2.max(1) as f64, m.3)).collect()
    }

    /// The months of a run of days already drawn, in order, read against this record: each
    /// whole month's rain, its SPI-3 (once two months precede it) and the soil's water.
    #[must_use]
    pub fn read(&self, days: &[Drawn]) -> Vec<Month> {
        let totals = self.totals(days);
        totals
            .iter()
            .enumerate()
            .map(|(i, &(k, mm, _, soil))| {
                let spi3 = (i >= 2).then(|| totals[i - 2].1 + totals[i - 1].1 + mm).and_then(|three| self.fits[k.rem_euclid(12) as usize].map(|(mean, sd)| ((three + 1.0).ln() - mean) / sd));
                Month { day: (k as f64 + 0.5) * self.month_days, rain_mm: mm, spi3, soil }
            })
            .collect()
    }
}

/// A place's droughts from `years` of its days already drawn, in order, months of
/// `month_days`: each whole month's rain, SPI-3 and soil water, and how often and how long
/// it was in drought.
#[must_use]
pub fn record(days: &[Drawn], month_days: f64, years: usize) -> Droughts {
    let (soil_now, snowpack_now_mm) = days.last().map_or((0.0, 0.0), |d| (d.soil, d.snowpack_mm));
    let blank = Droughts {
        years,
        moderate_share: 0.0,
        severe_share: 0.0,
        longest_months: 0,
        dry_soil_share: None,
        months: Vec::new(),
        soil_now,
        snowpack_now_mm,
        enso: None,
        first: 0,
        fits: vec![None; 12],
        month_days,
    };
    let totals = blank.totals(days);
    let Some(&(first, ..)) = totals.first() else { return blank };
    let count = totals.len();
    let three: Vec<f64> = (0..count).map(|i| if i < 2 { 0.0 } else { (i - 2..=i).map(|j| totals[j].1).sum() }).collect();
    let mut fits = vec![None; 12];
    for (cal, fit) in fits.iter_mut().enumerate() {
        let logs: Vec<f64> = (2..count).filter(|&i| (first + i as i64).rem_euclid(12) as usize == cal).map(|i| (three[i] + 1.0).ln()).collect();
        let mean_mm = (2..count).filter(|&i| (first + i as i64).rem_euclid(12) as usize == cal).map(|i| three[i]).sum::<f64>() / logs.len().max(1) as f64;
        if logs.len() < 3 || mean_mm < 10.0 {
            continue;
        }
        let mean = logs.iter().sum::<f64>() / logs.len() as f64;
        let sd = (logs.iter().map(|l| (l - mean).powi(2)).sum::<f64>() / (logs.len() - 1) as f64).sqrt().max(1e-6);
        *fit = Some((mean, sd));
    }
    let d = Droughts { fits, first, ..blank };
    let months = d.read(days);
    let read: Vec<f64> = months.iter().filter_map(|m| m.spi3).collect();
    let share = |lim: f64| if read.is_empty() { 0.0 } else { read.iter().filter(|&&s| s <= lim).count() as f64 / read.len() as f64 };
    let (mut longest, mut run) = (0, 0);
    for m in &months {
        run = if m.spi3.is_some_and(|s| s <= MODERATE_SPI) { run + 1 } else { 0 };
        longest = longest.max(run);
    }
    let growing: Vec<(&Month, f64)> = months.iter().zip(&totals).filter(|(_, t)| t.2 >= 5.0).map(|(m, t)| (m, t.2)).collect();
    let dry_soil_share = (!growing.is_empty()).then(|| growing.iter().filter(|(m, _)| m.soil < 0.25).count() as f64 / growing.len() as f64);
    Droughts { moderate_share: share(MODERATE_SPI), severe_share: share(SEVERE_SPI), longest_months: longest, dry_soil_share, months, ..d }
}

/// The droughts at (`lat`, `lon`) over `years` years ending with the month of Earth day
/// `now`.
///
/// # Errors
/// If the climate is not generated yet.
pub fn droughts(world: &World, lat: f64, lon: f64, now: f64, years: usize) -> Result<Droughts> {
    let body = world.load_body()?;
    let calendar = Calendar::new(&body, &[]);
    let step = calendar.solar_day_h.map_or(1.0, |h| h / 24.0);
    let span = years as f64 * calendar.year_days;
    // A month more at the start, so the first whole month is whole.
    let from = now - span - calendar.year_days / 12.0;
    let f = crate::weather::forecast(world, lat, lon, from, ((now - from) / step).ceil() as usize + 1)?;
    let days: Vec<Drawn> = f.days.iter().map(|d| Drawn::of(d, step)).collect();
    Ok(Droughts { enso: f.enso, ..record(&days, calendar.year_days / 12.0, years) })
}

/// Day length and sunrise at latitude `lat` with the star at declination `dec` (degrees), on
/// a world whose solar day is `day_h` hours: (sunrise, sunset) in local hours from midnight,
/// `None` in polar night, and (0, `day_h`) under the midnight sun. The sun's centre counts as
/// up 0.83° below the horizon (its radius and refraction).
#[must_use]
pub fn sun_times(lat: f64, dec: f64, day_h: f64) -> Option<(f64, f64)> {
    let (phi, d) = (lat.to_radians(), dec.to_radians());
    let cos_h = ((-0.833f64).to_radians().sin() - phi.sin() * d.sin()) / (phi.cos() * d.cos()).max(1e-12);
    if cos_h >= 1.0 {
        return None;
    }
    let h = cos_h.max(-1.0).acos() / PI;
    Some((0.5 * day_h * (1.0 - h), 0.5 * day_h * (1.0 + h)))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn earth_season(t: f64) -> f64 {
        (TAU * (t - 80.0) / 365.25).sin()
    }

    fn climate(lat: f64, precip: f64) -> Climate {
        Climate { mean_c: 15.0, range_c: 15.0, precipitation_mm: precip, warm_share: 0.5, at_sea: false, wind_jan: [0.0; 2], wind_jul: [0.0; 2], lat }
    }

    #[test]
    fn the_slow_factor_keeps_the_mean_and_has_wet_and_dry_years() {
        let w = Wetness::new(Seed(4), 10.0, 30.0, 800.0, 365.25, Some(PACIFIC));
        let f: Vec<f64> = (0..365 * 40).map(|d| w.factor(f64::from(d))).collect();
        let mean = f.iter().sum::<f64>() / f.len() as f64;
        assert!((mean - 1.0).abs() < 0.08, "{mean}");
        let years: Vec<f64> = f.chunks(365).map(|y| y.iter().sum::<f64>() / 365.0).collect();
        let (lo, hi) = years.iter().fold((f64::MAX, f64::MIN), |(lo, hi), &y| (lo.min(y), hi.max(y)));
        assert!(lo < 0.85 && hi > 1.15, "driest year {lo:.2}, wettest {hi:.2}");
        // Neighbours share their years; the far side of the world does not.
        let near = Wetness::new(Seed(4), 11.0, 31.0, 800.0, 365.25, Some(PACIFIC));
        let far = Wetness::new(Seed(4), -40.0, -150.0, 800.0, 365.25, Some(PACIFIC));
        let corr = |o: &Wetness| {
            let (a, b): (Vec<f64>, Vec<f64>) = (0..480).map(|m| (w.index(f64::from(m) * 30.4), o.index(f64::from(m) * 30.4))).unzip();
            let (ma, mb) = (a.iter().sum::<f64>() / 480.0, b.iter().sum::<f64>() / 480.0);
            let cov: f64 = a.iter().zip(&b).map(|(x, y)| (x - ma) * (y - mb)).sum();
            cov / (a.iter().map(|x| (x - ma).powi(2)).sum::<f64>() * b.iter().map(|y| (y - mb).powi(2)).sum::<f64>()).sqrt()
        };
        assert!(corr(&near) > 0.9 && corr(&far).abs() < 0.4, "near {:.2}, far {:.2}", corr(&near), corr(&far));
    }

    /// `years` of calm daily weather (no storms) ending at day 0, under the slow factor.
    fn days(c: &Climate, w: &Wetness, years: usize, seed: u64) -> Vec<Drawn> {
        let n = (years as f64 * 365.25 + 31.0) as usize;
        crate::weather::simulate(c, earth_season, 365.25, 1.0, -(n as f64), n, Seed(seed), &|_| crate::storms::Effect::default(), &crate::weather::Storminess::calm(), &[], &|t| w.factor(t))
            .iter()
            .map(|d| Drawn::of(d, 1.0))
            .collect()
    }

    #[test]
    fn drought_index_reads_about_a_sixth_of_months_as_drought() {
        let c = climate(40.0, 700.0);
        let w = Wetness::new(Seed(7), 40.0, 10.0, 700.0, 365.25, Some(PACIFIC));
        let d = record(&days(&c, &w, 40, 8), 365.25 / 12.0, 40);
        assert!(d.months.len() >= 479, "{}", d.months.len());
        // A standardized index: about 16% of months at −1 or below, 7% at −1.5.
        assert!(d.moderate_share > 0.1 && d.moderate_share < 0.22, "{}", d.moderate_share);
        assert!(d.severe_share > 0.03 && d.severe_share < d.moderate_share, "{}", d.severe_share);
        assert!(d.longest_months >= 3);
        assert!(d.dry_soil_share.is_some());
        // Days drawn later read against it: a dry year reads as drought, a wet one does not.
        let run = |mm: f64| {
            let mut soil = 0.5 * SOIL_MM;
            (0..400)
                .map(|t| {
                    soil = bucket(soil, 0.0, mm, false, 12.0, 1.0).0;
                    Drawn { day: f64::from(t), rain_mm: mm, mean_c: 12.0, soil: soil / SOIL_MM, snowpack_mm: 0.0 }
                })
                .collect::<Vec<_>>()
        };
        let (dry, wet) = (d.read(&run(0.2)), d.read(&run(4.0)));
        assert_eq!((dry.len(), dry.iter().filter(|m| m.spi3.is_some()).count()), (13, 11));
        assert!(dry.iter().filter_map(|m| m.spi3).all(|s| s < -2.0) && wet.iter().filter_map(|m| m.spi3).all(|s| s > 0.5));
        assert!(dry.last().unwrap().soil < wet.last().unwrap().soil);
    }

    const PACIFIC: Basin = Basin { width_km: EARTH_BASIN_KM, west_lon: 118.0, east_lon: -80.0 };

    #[test]
    fn the_widest_equatorial_sea() {
        // Land from 80°W east to 123°E, as the Americas, Africa and Asia; a small island at
        // 150°W; the rest sea: the basin is the Pacific's 157° of the equator.
        let land = |la: f64, lo: f64| (-80.0..123.0).contains(&lo) || ((-150.5..-149.5).contains(&lo) && la.abs() < 1.0);
        let b = equatorial_basin(|la, lo| !land(la, lo), 6371.0).unwrap();
        assert!((b.west_lon - 123.0).abs() < 0.3 && (b.east_lon + 80.0).abs() < 0.3, "{b:?}");
        assert!((b.width_km - 6371.0 * 157f64.to_radians()).abs() < 50.0, "{b:?}");
        // No shore or no sea: no basin.
        assert!(equatorial_basin(|_, _| true, 6371.0).is_none() && equatorial_basin(|_, _| false, 6371.0).is_none());
        // A narrow basin has a weak, quick cycle; a wide one a strong, slow one.
        let at = |w: f64| Wetness::new(Seed(4), 0.0, 0.0, 800.0, 365.25, Some(Basin { width_km: w, ..PACIFIC }));
        let (narrow, wide) = (at(5_000.0), at(EARTH_BASIN_KM));
        assert!(narrow.cycle_share().abs() < 0.1 * wide.cycle_share().abs());
        assert!((narrow.periods_years()[0] / wide.periods_years()[0] - 5_000.0 / EARTH_BASIN_KM).abs() < 1e-9);
        // The warm phase wets the basin's east and dries its west, whatever the seed.
        for seed in 0..20 {
            let here = |lat: f64, lon: f64| Wetness::new(Seed(seed), lat, lon, 800.0, 365.25, Some(PACIFIC)).cycle_share();
            assert!(here(-2.0, -79.9) > 0.3 && here(-6.2, 106.8) < -0.3, "Guayaquil {} Jakarta {}", here(-2.0, -79.9), here(-6.2, 106.8));
            assert!(here(-3.1, -60.0) < -0.3 && here(-12.5, 130.8) < -0.2, "Manaus {} Darwin {}", here(-3.1, -60.0), here(-12.5, 130.8));
        }
        assert_eq!(Wetness::new(Seed(4), 0.0, 0.0, 800.0, 365.25, None).cycle_share(), 0.0);
    }

    #[test]
    fn the_bucket_fills_dries_and_holds_snow() {
        // Rain fills it to capacity; warm dry days empty it, ever slower.
        assert_eq!(bucket(140.0, 0.0, 30.0, false, 10.0, 1.0).0, SOIL_MM);
        let (a, _) = bucket(SOIL_MM, 0.0, 0.0, false, 20.0, 1.0);
        let (b, _) = bucket(0.2 * SOIL_MM, 0.0, 0.0, false, 20.0, 1.0);
        assert!((SOIL_MM - a - 4.14).abs() < 1e-9 && 0.2 * SOIL_MM - b < 1.0);
        // Snow lies until it melts, 3 mm a degree a day, then soaks in.
        let (soil, pack) = bucket(50.0, 0.0, 20.0, true, -5.0, 1.0);
        assert_eq!((soil, pack), (50.0, 20.0));
        let (soil, pack) = bucket(soil, pack, 0.0, false, 2.0, 1.0);
        assert!((pack - 14.0).abs() < 1e-9 && soil > 55.0);
    }

    #[test]
    fn a_continental_winter_keeps_its_snow_until_spring() {
        let c = Climate { mean_c: 2.0, range_c: 30.0, ..climate(55.0, 500.0) };
        let d = crate::weather::simulate(&c, earth_season, 365.25, 1.0, 0.0, 365 * 3, Seed(5), &|_| crate::storms::Effect::default(), &crate::weather::Storminess::calm(), &[], &|_| 1.0);
        // `earth_season` has the solstice at day −11, so day 0 is New Year's: the snow of
        // February is gone by June, and the soil after the melt is wetter than at summer's end.
        let at = |day: usize| &d[365 + day];
        let mean = |f: fn(&crate::weather::Day) -> f64, from: usize, to: usize| (from..to).map(|k| f(at(k))).sum::<f64>() / (to - from) as f64;
        assert!(mean(|x| x.snowpack_mm, 30, 60) > 50.0, "February snow {}", mean(|x| x.snowpack_mm, 30, 60));
        assert!(mean(|x| x.snowpack_mm, 150, 250) < 1.0);
        assert!(mean(|x| x.soil, 105, 135) > mean(|x| x.soil, 240, 270) + 0.2, "spring {} vs late summer {}", mean(|x| x.soil, 105, 135), mean(|x| x.soil, 240, 270));
    }

    #[test]
    fn a_desert_is_too_dry_to_read() {
        let c = Climate { mean_c: 25.0, precipitation_mm: 5.0, ..climate(22.0, 5.0) };
        let w = Wetness::new(Seed(7), 22.0, 10.0, 5.0, 365.25, Some(PACIFIC));
        let d = record(&days(&c, &w, 20, 8), 365.25 / 12.0, 20);
        assert!(d.months.iter().all(|m| m.spi3.is_none()));
        assert!(d.dry_soil_share.is_some_and(|s| s > 0.9));
    }

    #[test]
    fn the_sun_rises_later_in_winter_and_not_at_all_in_polar_night() {
        let (rise, set) = sun_times(0.0, 0.0, 24.0).unwrap();
        assert!((set - rise - 12.1).abs() < 0.1, "{rise} {set}");
        let (s_rise, s_set) = sun_times(51.5, 23.44, 24.0).unwrap();
        let (w_rise, w_set) = sun_times(51.5, -23.44, 24.0).unwrap();
        assert!((s_set - s_rise - 16.5).abs() < 0.3 && (w_set - w_rise - 7.9).abs() < 0.3, "{s_rise} {s_set} {w_rise} {w_set}");
        assert!(sun_times(80.0, -20.0, 24.0).is_none());
        assert_eq!(sun_times(80.0, 20.0, 24.0), Some((0.0, 24.0)));
    }
}
