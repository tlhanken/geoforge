//! Wet and dry years at a place: a slow rain anomaly, the drought index read from it, and the
//! soil's water.
//!
//! - **The slow anomaly.** The rain of any stretch of days is the climate's normal times a
//!   slow factor, whose logarithm is a month-to-month AR(1) (a dry spell keeps 60% of itself
//!   from one month to the next, so it outlasts a season) plus an ENSO-like cycle. The AR(1)
//!   draws are read off a 15° lattice, so neighbours share their wet and dry years. The cycle
//!   is one world-wide oscillation of two seeded periods between 2.5 and 7 years, felt most
//!   in the tropics, with the sign of a region's response set by a broad seeded pattern (as
//!   El Niño dries Indonesia and wets Peru). Drier climates vary more: the logarithm's spread
//!   runs from 0.2 in wet climates to 0.5 in deserts.
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
/// Soil water capacity, mm.
const SOIL_MM: f64 = 150.0;
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
    /// days long, for a place with `annual_mm` of rain.
    #[must_use]
    pub fn new(seed: Seed, lat: f64, lon: f64, annual_mm: f64, year_days: f64) -> Self {
        let draw = |label: &[u8]| unit(seed.child(b"cycle").child(label));
        let earth_year = crate::exposure::EARTH_YEAR_DAYS;
        let periods = [2.5 + 2.5 * draw(b"p1"), 4.0 + 3.0 * draw(b"p2")].map(|y| y * earth_year);
        let pattern = lattice(seed.child(b"pattern"), lat, lon, PATTERN_DEG);
        let cycle = CYCLE_SHARE * (-(lat / 30.0).powi(2)).exp() * (1.5 * pattern).tanh();
        Self { seed, lat, lon, spread: 0.2 + 0.3 * (-annual_mm.max(0.0) / 500.0).exp(), cycle, periods, phases: [TAU * draw(b"f1"), TAU * draw(b"f2")], month_days: year_days / 12.0 }
    }

    /// The wetness of a place on `world`, whose climate is `climate` and year `year_days`.
    #[must_use]
    pub fn of(world: &World, lat: f64, lon: f64, climate: &Climate, year_days: f64) -> Self {
        Self::new(world.manifest.root_seed.child(b"drought"), lat, lon, climate.precipitation_mm, year_days)
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
    /// The whole months of a run of days (`(Earth day, rain mm, mean °C)` each, in order):
    /// (the months' index, rain mm, mean °C).
    fn totals(&self, days: &[(f64, f64, f64)]) -> Vec<(i64, f64, f64)> {
        let index = |t: f64| (t / self.month_days).floor() as i64;
        let (Some(a), Some(b)) = (days.first(), days.last()) else { return Vec::new() };
        let (from, to) = (index(a.0 - 1e-9) + 1, index(b.0 + 1e-9) - 1);
        let mut sums: Vec<(f64, f64, usize)> = vec![(0.0, 0.0, 0); (to - from + 1).max(0) as usize];
        for &(t, mm, c) in days {
            let k = index(t) - from;
            if (0..sums.len() as i64).contains(&k) {
                let m = &mut sums[k as usize];
                (m.0, m.1, m.2) = (m.0 + mm, m.1 + c, m.2 + 1);
            }
        }
        sums.iter().zip(from..).map(|(m, k)| (k, m.0, m.1 / m.2.max(1) as f64)).collect()
    }

    /// The months of a run of days already drawn (`(Earth day, rain mm, mean °C)` each, in
    /// order), read against this record: each whole month's rain, its SPI-3 (once two
    /// months precede it) and the soil's water, carried on from the record's.
    #[must_use]
    pub fn read(&self, days: &[(f64, f64, f64)]) -> Vec<Month> {
        let totals = self.totals(days);
        let Some(&(from, ..)) = totals.first() else { return Vec::new() };
        let mut soil = self.months.iter().zip(self.first..).find(|(_, k)| *k == from - 1).map_or(0.5, |(m, _)| m.soil) * SOIL_MM;
        totals
            .iter()
            .enumerate()
            .map(|(i, &(k, mm, temp))| {
                soil = (soil + mm - self.month_days * (0.18 * (temp + 3.0)).max(0.0) * soil / SOIL_MM).clamp(0.0, SOIL_MM);
                let spi3 = (i >= 2).then(|| totals[i - 2].1 + totals[i - 1].1 + mm).and_then(|three| self.fits[k.rem_euclid(12) as usize].map(|(mean, sd)| ((three + 1.0).ln() - mean) / sd));
                Month { day: (k as f64 + 0.5) * self.month_days, rain_mm: mm, spi3, soil: soil / SOIL_MM }
            })
            .collect()
    }
}

/// A place's droughts from `years` of its days already drawn (`(Earth day, rain mm, mean
/// °C)` each, in order), months of `month_days`: each whole month's rain, SPI-3 and soil
/// water, and how often and how long it was in drought.
#[must_use]
pub fn record(days: &[(f64, f64, f64)], month_days: f64, years: usize) -> Droughts {
    let blank = Droughts { years, moderate_share: 0.0, severe_share: 0.0, longest_months: 0, dry_soil_share: None, months: Vec::new(), first: 0, fits: vec![None; 12], month_days };
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
    let days: Vec<(f64, f64, f64)> = f.days.iter().map(|d| (d.day + 0.5 * step, d.precipitation_mm, d.mean_c)).collect();
    Ok(record(&days, calendar.year_days / 12.0, years))
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
        let w = Wetness::new(Seed(4), 10.0, 30.0, 800.0, 365.25);
        let f: Vec<f64> = (0..365 * 40).map(|d| w.factor(f64::from(d))).collect();
        let mean = f.iter().sum::<f64>() / f.len() as f64;
        assert!((mean - 1.0).abs() < 0.08, "{mean}");
        let years: Vec<f64> = f.chunks(365).map(|y| y.iter().sum::<f64>() / 365.0).collect();
        let (lo, hi) = years.iter().fold((f64::MAX, f64::MIN), |(lo, hi), &y| (lo.min(y), hi.max(y)));
        assert!(lo < 0.85 && hi > 1.15, "driest year {lo:.2}, wettest {hi:.2}");
        // Neighbours share their years; the far side of the world does not.
        let near = Wetness::new(Seed(4), 11.0, 31.0, 800.0, 365.25);
        let far = Wetness::new(Seed(4), -40.0, -150.0, 800.0, 365.25);
        let corr = |o: &Wetness| {
            let (a, b): (Vec<f64>, Vec<f64>) = (0..480).map(|m| (w.index(f64::from(m) * 30.4), o.index(f64::from(m) * 30.4))).unzip();
            let (ma, mb) = (a.iter().sum::<f64>() / 480.0, b.iter().sum::<f64>() / 480.0);
            let cov: f64 = a.iter().zip(&b).map(|(x, y)| (x - ma) * (y - mb)).sum();
            cov / (a.iter().map(|x| (x - ma).powi(2)).sum::<f64>() * b.iter().map(|y| (y - mb).powi(2)).sum::<f64>()).sqrt()
        };
        assert!(corr(&near) > 0.9 && corr(&far).abs() < 0.4, "near {:.2}, far {:.2}", corr(&near), corr(&far));
    }

    /// `years` of calm daily weather (no storms) ending at day 0, under the slow factor.
    fn days(c: &Climate, w: &Wetness, years: usize, seed: u64) -> Vec<(f64, f64, f64)> {
        let n = (years as f64 * 365.25 + 31.0) as usize;
        crate::weather::simulate(c, earth_season, 365.25, 1.0, -(n as f64), n, Seed(seed), &|_| crate::storms::Effect::default(), &crate::weather::Storminess::calm(), &[], &|t| w.factor(t))
            .iter()
            .map(|d| (d.day + 0.5, d.precipitation_mm, d.mean_c))
            .collect()
    }

    #[test]
    fn drought_index_reads_about_a_sixth_of_months_as_drought() {
        let c = climate(40.0, 700.0);
        let w = Wetness::new(Seed(7), 40.0, 10.0, 700.0, 365.25);
        let d = record(&days(&c, &w, 40, 8), 365.25 / 12.0, 40);
        assert!(d.months.len() >= 479, "{}", d.months.len());
        // A standardized index: about 16% of months at −1 or below, 7% at −1.5.
        assert!(d.moderate_share > 0.1 && d.moderate_share < 0.22, "{}", d.moderate_share);
        assert!(d.severe_share > 0.03 && d.severe_share < d.moderate_share, "{}", d.severe_share);
        assert!(d.longest_months >= 3);
        assert!(d.dry_soil_share.is_some());
        // Days drawn later read against it: a dry year reads as drought, a wet one does not.
        let run = |mm: f64| (0..400).map(|t| (f64::from(t), mm, 12.0)).collect::<Vec<_>>();
        let (dry, wet) = (d.read(&run(0.2)), d.read(&run(4.0)));
        assert_eq!((dry.len(), dry.iter().filter(|m| m.spi3.is_some()).count()), (13, 11));
        assert!(dry.iter().filter_map(|m| m.spi3).all(|s| s < -2.0) && wet.iter().filter_map(|m| m.spi3).all(|s| s > 0.5));
        assert!(dry.last().unwrap().soil < wet.last().unwrap().soil);
    }

    #[test]
    fn a_desert_is_too_dry_to_read() {
        let c = Climate { mean_c: 25.0, precipitation_mm: 5.0, ..climate(22.0, 5.0) };
        let w = Wetness::new(Seed(7), 22.0, 10.0, 5.0, 365.25);
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
