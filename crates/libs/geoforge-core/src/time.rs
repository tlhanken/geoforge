//! Time: one representation under every calendar.
//!
//! - [`Instant`] and [`Span`] hold SI seconds: an instant from the epoch (the focus world's
//!   northern spring equinox, "the present"), a span for lengths. `f64` seconds resolve
//!   better than a millisecond over ±10,000 years. Earth days and hours are conversions.
//! - A world's [`Clock`] holds its constants: the local day (the mean solar day; none for a
//!   world locked to its star, which counts Earth days) and the local year (one orbit).
//!   Year `n` starts `n` years from the epoch; year 0 is the present one.
//! - A [`Calendar`] only divides each year into months and names them: [`Twelfths`] (twelve
//!   equal months, what climate statistics use), [`Lunar`] (a moon's lunations, new moon to
//!   new moon) and [`Custom`] (named months of whole local days, as a world's
//!   `calendar.json` gives them). No calendar changes how long a second, a day or a year is.
//! - Dates are whole days counted from the epoch, the same everywhere (like UTC), so the
//!   time of day keeps step with the sun. A day belongs to the year and the month in which
//!   it begins: as the year is rarely a whole number of days, calendar years differ by a
//!   day (Earth's are 365 or 366 days long, its leap years), while the clock's year stays
//!   one orbit.
//!
//! Deep time (rock and crust ages, plate motion) stays in Ma.

use crate::body::BodyContext;
use serde::{Deserialize, Serialize};
use std::ops::{Add, AddAssign, Div, Mul, Neg, Sub};

/// Seconds in an Earth day.
pub const EARTH_DAY_S: f64 = 86_400.0;

/// A length of time, SI seconds.
#[derive(Debug, Clone, Copy, Default, PartialEq, PartialOrd, Serialize, Deserialize)]
#[serde(transparent)]
pub struct Span(f64);

impl Span {
    pub const ZERO: Self = Self(0.0);

    #[must_use]
    pub const fn from_seconds(s: f64) -> Self {
        Self(s)
    }
    #[must_use]
    pub fn from_hours(h: f64) -> Self {
        Self(h * 3600.0)
    }
    #[must_use]
    pub fn from_earth_days(d: f64) -> Self {
        Self(d * EARTH_DAY_S)
    }
    #[must_use]
    pub const fn seconds(self) -> f64 {
        self.0
    }
    #[must_use]
    pub fn hours(self) -> f64 {
        self.0 / 3600.0
    }
    #[must_use]
    pub fn earth_days(self) -> f64 {
        self.0 / EARTH_DAY_S
    }
}

impl Add for Span {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self(self.0 + o.0)
    }
}
impl AddAssign for Span {
    fn add_assign(&mut self, o: Self) {
        self.0 += o.0;
    }
}
impl Sub for Span {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Self(self.0 - o.0)
    }
}
impl Neg for Span {
    type Output = Self;
    fn neg(self) -> Self {
        Self(-self.0)
    }
}
impl Mul<f64> for Span {
    type Output = Self;
    fn mul(self, k: f64) -> Self {
        Self(self.0 * k)
    }
}
impl Div<f64> for Span {
    type Output = Self;
    fn div(self, k: f64) -> Self {
        Self(self.0 / k)
    }
}
/// How many of one span fit in another.
impl Div for Span {
    type Output = f64;
    fn div(self, o: Self) -> f64 {
        self.0 / o.0
    }
}

/// A moment, SI seconds from the epoch (the focus world's northern spring equinox).
#[derive(Debug, Clone, Copy, Default, PartialEq, PartialOrd, Serialize, Deserialize)]
#[serde(transparent)]
pub struct Instant(f64);

impl Instant {
    /// The present: the focus world's northern spring equinox.
    pub const EPOCH: Self = Self(0.0);

    #[must_use]
    pub const fn from_seconds(s: f64) -> Self {
        Self(s)
    }
    #[must_use]
    pub fn from_earth_days(d: f64) -> Self {
        Self(d * EARTH_DAY_S)
    }
    /// Seconds from the epoch.
    #[must_use]
    pub const fn seconds(self) -> f64 {
        self.0
    }
    /// Earth days from the epoch.
    #[must_use]
    pub fn earth_days(self) -> f64 {
        self.0 / EARTH_DAY_S
    }
}

impl Add<Span> for Instant {
    type Output = Self;
    fn add(self, s: Span) -> Self {
        Self(self.0 + s.0)
    }
}
impl AddAssign<Span> for Instant {
    fn add_assign(&mut self, s: Span) {
        self.0 += s.0;
    }
}
impl Sub<Span> for Instant {
    type Output = Self;
    fn sub(self, s: Span) -> Self {
        Self(self.0 - s.0)
    }
}
impl Sub for Instant {
    type Output = Span;
    fn sub(self, o: Self) -> Span {
        Span(self.0 - o.0)
    }
}

/// A world's constant lengths of time.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Clock {
    /// The mean solar day; `None` for a world locked to its star.
    pub day: Option<Span>,
    /// One orbit.
    pub year: Span,
}

impl Clock {
    #[must_use]
    pub const fn new(day: Option<Span>, year: Span) -> Self {
        Self { day, year }
    }

    /// The clock of `body`: its mean solar day `1 / (1/spin ∓ 1/year)` (a retrograde spin,
    /// obliquity past 90°, adds) and its year.
    #[must_use]
    pub fn of_body(body: &BodyContext) -> Self {
        let year = body.irradiance.year_days;
        let spin = body.exterior.rotation_h.abs() / 24.0;
        let turn = if body.exterior.obliquity_deg > 90.0 { 1.0 } else { -1.0 };
        let day = (!body.is_tidally_locked()).then(|| Span::from_earth_days(1.0 / (1.0 / spin + turn / year)));
        Self { day, year: Span::from_earth_days(year) }
    }

    /// The day dates count in: the solar day, or an Earth day on a locked world.
    #[must_use]
    pub fn day_len(&self) -> Span {
        self.day.unwrap_or(Span::from_earth_days(1.0))
    }

    /// Days (of [`Self::day_len`]) in a year.
    #[must_use]
    pub fn days_per_year(&self) -> f64 {
        self.year / self.day_len()
    }

    /// The year `t` falls in (0 the present one) and the time since that year began.
    #[must_use]
    pub fn year_of(&self, t: Instant) -> (i64, Span) {
        let n = ((t - Instant::EPOCH) / self.year).floor();
        // Rounding can leave `t` a hair before its year's start; keep it in the year before.
        let into = t - self.year_start(n as i64);
        if into < Span::ZERO {
            (n as i64 - 1, into + self.year)
        } else {
            (n as i64, into)
        }
    }

    /// When year `n` begins.
    #[must_use]
    pub fn year_start(&self, n: i64) -> Instant {
        Instant::EPOCH + self.year * n as f64
    }

    /// The day `t` falls in, counted from the epoch, and the share of it gone.
    #[must_use]
    pub fn day_count(&self, t: Instant) -> (i64, f64) {
        let d = (t - Instant::EPOCH) / self.day_len();
        // A day's start, worked out as `k × day`, reads back as day `k`.
        let k = (d + 1e-9).floor();
        (k as i64, (d - k).max(0.0))
    }

    /// When day `k` (from the epoch) begins.
    #[must_use]
    pub fn day_start(&self, k: i64) -> Instant {
        Instant::EPOCH + self.day_len() * k as f64
    }

    /// The first day that begins at or after `t`.
    #[must_use]
    pub fn first_day_from(&self, t: Instant) -> i64 {
        ((t - Instant::EPOCH) / self.day_len() - 1e-9).ceil() as i64
    }

    /// The year day `k` (from the epoch) belongs to: the one it begins in.
    #[must_use]
    pub fn year_of_day(&self, k: i64) -> i64 {
        let n = self.year_of(self.day_start(k)).0;
        if self.first_day_from(self.year_start(n + 1)) <= k {
            n + 1
        } else if self.first_day_from(self.year_start(n)) > k {
            n - 1
        } else {
            n
        }
    }
}

/// A month of a year: its name, when it starts after the year's start, and how long it is.
/// Its days are those that begin inside it.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Month {
    pub name: String,
    pub start: Span,
    pub len: Span,
}

/// A date: the year (0 the present one), the month (index into the year's
/// [`Calendar::months`]), the day of the month and the share of that day gone, days of
/// [`Clock::day_len`] counted from 0. A day's year and month are those it begins in.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct Date {
    pub year: i64,
    pub month: usize,
    pub day: u32,
    pub of_day: f64,
}

/// Months on a world's clock.
pub trait Calendar: std::fmt::Debug + Send + Sync {
    /// The clock the calendar counts on.
    fn clock(&self) -> &Clock;

    /// What the calendar is called.
    fn name(&self) -> String;

    /// The months of year `year`, in order from its start, filling it exactly.
    fn months(&self, year: i64) -> Vec<Month>;

    /// The names of the days of the week, if the calendar keeps one.
    fn week(&self) -> &[String] {
        &[]
    }

    /// How a year is written.
    fn year_name(&self, year: i64) -> String {
        format!("year {year}")
    }

    /// The first day (counted from the epoch) of each month of year `year`, then the next
    /// year's: months and years are counted in the same whole days, so rounding in a month's
    /// start never moves a day between them.
    fn first_days(&self, year: i64) -> Vec<i64> {
        let clock = self.clock();
        let start = clock.year_start(year);
        self.months(year).iter().map(|m| clock.first_day_from(start + m.start)).chain(std::iter::once(clock.first_day_from(clock.year_start(year + 1)))).collect()
    }

    /// The date of `t`.
    fn date(&self, t: Instant) -> Date {
        let clock = self.clock();
        let (k, of_day) = clock.day_count(t);
        let year = clock.year_of_day(k);
        let firsts = self.first_days(year);
        let month = firsts[..firsts.len() - 1].iter().rposition(|&f| f <= k).unwrap_or(0);
        Date { year, month, day: (k - firsts[month]).max(0) as u32, of_day }
    }

    /// When `date` begins (its `of_day` ignored); `None` if its month or day is not in its year.
    fn start_of(&self, date: &Date) -> Option<Instant> {
        let firsts = self.first_days(date.year);
        let (&first, &next) = (firsts.get(date.month)?, firsts.get(date.month + 1)?);
        let k = first + i64::from(date.day);
        (k < next).then(|| self.clock().day_start(k))
    }

    /// The date of `t` written out: the weekday if the calendar keeps a week, the day and
    /// month, and the year.
    fn format_date(&self, t: Instant) -> String {
        let d = self.date(t);
        let month = self.months(d.year).into_iter().nth(d.month).map_or_else(String::new, |m| m.name);
        let week = self.week();
        let weekday = if week.is_empty() { String::new() } else { format!("{}, ", week[self.clock().day_count(t).0.rem_euclid(week.len() as i64) as usize]) };
        format!("{weekday}{} {month}, {}", d.day + 1, self.year_name(d.year))
    }

    /// `t` written out: its date and the time of day in local hours (24 to the day) at the
    /// prime meridian.
    fn format(&self, t: Instant) -> String {
        let minutes = (self.date(t).of_day * 1440.0).floor() as u32;
        format!("{}, {:02}:{:02}", self.format_date(t), minutes / 60, minutes % 60)
    }
}

/// Twelve equal months, each a twelfth of the year: what the seasons follow.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Twelfths(pub Clock);

impl Calendar for Twelfths {
    fn clock(&self) -> &Clock {
        &self.0
    }
    fn name(&self) -> String {
        "twelfths".into()
    }
    fn months(&self, _year: i64) -> Vec<Month> {
        twelfths(&self.0)
    }
}

fn twelfths(clock: &Clock) -> Vec<Month> {
    let len = clock.year / 12.0;
    (0..12).map(|i| Month { name: format!("month {}", i + 1), start: len * f64::from(i), len }).collect()
}

/// A moon's lunations, new moon to new moon, inside the clock's years: month `n` begins at
/// the year's `n`-th new moon, and the days before its first (if any) are month 0, so the
/// first and last months are cut by the year's ends. A moon whose month is longer than the
/// year or no longer than the day falls back to [`Twelfths`].
#[derive(Debug, Clone, PartialEq)]
pub struct Lunar {
    pub clock: Clock,
    pub moon: String,
    /// New moon to new moon.
    pub synodic: Span,
    /// Any new moon.
    pub new_moon: Instant,
}

impl Lunar {
    /// Does the moon's month fit between the day and the year?
    #[must_use]
    pub fn fits(&self) -> bool {
        self.synodic > self.clock.day_len() && self.synodic < self.clock.year
    }
}

impl Calendar for Lunar {
    fn clock(&self) -> &Clock {
        &self.clock
    }
    fn name(&self) -> String {
        if self.fits() {
            format!("lunar ({})", self.moon)
        } else {
            format!("twelfths ({}'s month does not fit the year)", self.moon)
        }
    }
    fn months(&self, year: i64) -> Vec<Month> {
        if !self.fits() {
            return twelfths(&self.clock);
        }
        let start = self.clock.year_start(year);
        // The year's new moons, from its start.
        let mut moons = Vec::new();
        let mut moon = self.new_moon + self.synodic * ((start - self.new_moon) / self.synodic).ceil();
        while moon - start < self.clock.year {
            moons.push(moon - start);
            moon += self.synodic;
        }
        // Month 0 holds the days that begin before the first new moon, if any do; otherwise
        // the first lunation's month takes that time too.
        let lead = self.clock.first_day_from(start + moons[0]) > self.clock.first_day_from(start);
        let starts: Vec<Span> = std::iter::once(Span::ZERO).chain(moons.into_iter().skip(usize::from(!lead))).collect();
        let ends = starts.iter().skip(1).copied().chain(std::iter::once(self.clock.year));
        starts.iter().zip(ends).enumerate().map(|(i, (&a, b))| Month { name: format!("{} {}", self.moon, i + usize::from(!lead)), start: a, len: b - a }).collect()
    }
}

/// A fantasy calendar as a world's `calendar.json` writes it.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CalendarDef {
    pub name: String,
    /// The months in order, lengths in local days (Earth days on a locked world). The last
    /// month takes whatever the others leave of the year; its own length is ignored.
    pub months: Vec<MonthDef>,
    /// Names of the week's days, counted from the epoch.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub week: Vec<String>,
    /// The number of the present year (default 0) and the era written after it.
    #[serde(default)]
    pub present_year: i64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub era: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MonthDef {
    pub name: String,
    #[serde(default)]
    pub days: f64,
}

/// A named calendar of months of whole local days ([`CalendarDef`]).
#[derive(Debug, Clone, PartialEq)]
pub struct Custom {
    pub clock: Clock,
    pub def: CalendarDef,
}

impl Custom {
    /// # Errors
    /// If there are no months, a length is not a whole number of days, or the months before the last fill
    /// the year (leaving the last less than half a day).
    pub fn new(clock: Clock, def: CalendarDef) -> Result<Self, String> {
        let Some((_, before)) = def.months.split_last() else { return Err(format!("calendar {}: no months", def.name)) };
        if let Some(m) = before.iter().find(|m| !(m.days >= 1.0 && m.days.fract() == 0.0)) {
            return Err(format!("calendar {}: month {} must be a whole number of days, not {}", def.name, m.name, m.days));
        }
        let left = clock.days_per_year() - before.iter().map(|m| m.days).sum::<f64>();
        if left < 0.5 {
            return Err(format!("calendar {}: the months before the last take {:.2} of the year's {:.2} days", def.name, clock.days_per_year() - left, clock.days_per_year()));
        }
        Ok(Self { clock, def })
    }
}

impl Calendar for Custom {
    fn clock(&self) -> &Clock {
        &self.clock
    }
    fn name(&self) -> String {
        self.def.name.clone()
    }
    fn months(&self, _year: i64) -> Vec<Month> {
        let day = self.clock.day_len();
        let mut start = Span::ZERO;
        let n = self.def.months.len();
        self.def
            .months
            .iter()
            .enumerate()
            .map(|(i, m)| {
                let len = if i + 1 == n { self.clock.year - start } else { day * m.days };
                let month = Month { name: m.name.clone(), start, len };
                start += len;
                month
            })
            .collect()
    }
    fn week(&self) -> &[String] {
        &self.def.week
    }
    fn year_name(&self, year: i64) -> String {
        let n = self.def.present_year + year;
        self.def.era.as_ref().map_or_else(|| format!("year {n}"), |era| format!("{n} {era}"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Archetype, Key};

    fn earth() -> Clock {
        Clock::of_body(&Archetype::Earth.context(Key(1), "Earth"))
    }

    #[test]
    fn earth_days_convert_exactly() {
        let t = Instant::from_earth_days(1.5);
        assert_eq!(t.seconds(), 129_600.0);
        assert_eq!(t.earth_days(), 1.5);
        assert_eq!((t + Span::from_hours(12.0)).earth_days(), 2.0);
        assert_eq!(Span::from_earth_days(2.0) / Span::from_hours(12.0), 4.0);
    }

    #[test]
    fn clocks_of_worlds() {
        let e = earth();
        assert!((e.day.unwrap().hours() - 24.0).abs() < 0.01 && (e.days_per_year() - 365.25).abs() < 0.02, "{e:?}");
        let mars = Clock::of_body(&Archetype::Mars.context(Key(1), "Mars"));
        assert!((mars.day.unwrap().hours() - 24.66).abs() < 0.01);
        let venus = Clock::of_body(&Archetype::Venus.context(Key(1), "Venus"));
        assert!((venus.day.unwrap().earth_days() - 116.75).abs() < 0.1, "retrograde");
        let (y, into) = e.year_of(Instant::from_earth_days(400.0));
        assert!(y == 1 && (into.earth_days() - 34.7).abs() < 0.1);
        assert_eq!(e.year_of(Instant::from_earth_days(-1.0)).0, -1);
        assert_eq!(e.year_of(e.year_start(7)).0, 7);
    }

    #[test]
    fn a_locked_world_counts_earth_days() {
        let mut b = Archetype::Earth.context(Key(1), "x");
        b.exterior.rotation_h = b.irradiance.year_days * 24.0;
        let c = Clock::of_body(&b);
        assert!(c.day.is_none());
        assert_eq!(c.day_len(), Span::from_earth_days(1.0));
        let d = Twelfths(c).date(Instant::from_earth_days(40.25));
        assert_eq!((d.month, d.day), (1, 9));
    }

    fn check(cal: &dyn Calendar, years: std::ops::Range<i64>) {
        let clock = *cal.clock();
        for y in years {
            let months = cal.months(y);
            let total: f64 = months.iter().map(|m| m.len.seconds()).sum();
            assert!((total - clock.year.seconds()).abs() < 1e-3, "{} year {y}: {total}", cal.name());
            assert!(months.iter().all(|m| m.len > Span::ZERO));
            assert!(months.windows(2).all(|w| (w[0].start + w[0].len - w[1].start).seconds().abs() < 1e-6));
        }
        // Round trips, and every calendar agrees with the clock on the year.
        for k in -400..400 {
            let t = Instant::from_earth_days(f64::from(k) * 9.37 + 0.3);
            let d = cal.date(t);
            assert_eq!(d.year, clock.year_of_day(clock.day_count(t).0), "a day is in the year it begins in");
            let start = cal.start_of(&d).unwrap();
            assert!(start <= t && t - start < clock.day_len(), "{} {d:?}", cal.name());
            assert!((((t - start) / clock.day_len()) - d.of_day).abs() < 1e-6);
        }
    }

    #[test]
    fn twelfths() {
        let c = Twelfths(earth());
        check(&c, -3..3);
        let months = c.months(0);
        assert_eq!(months.len(), 12);
        assert!((months[0].len.earth_days() - 30.44).abs() < 0.01);
        assert!(c.start_of(&Date { year: 0, month: 12, day: 0, of_day: 0.0 }).is_none());
        assert!(c.start_of(&Date { year: 0, month: 0, day: 31, of_day: 0.0 }).is_none());
        // Whole days: calendar years of 365 or 366 days on the constant orbit, a leap year about one in four.
        let clock = *c.clock();
        let lengths: Vec<i64> = (0..8).map(|y| clock.first_day_from(clock.year_start(y + 1)) - clock.first_day_from(clock.year_start(y))).collect();
        let total: i64 = lengths.iter().sum();
        assert!(lengths.iter().all(|&n| n == 365 || n == 366) && (total as f64 - 8.0 * clock.days_per_year()).abs() < 1.0, "{lengths:?}");
        let even = Twelfths(Clock::new(Some(Span::from_hours(24.0)), Span::from_earth_days(360.0)));
        assert_eq!(even.format(Instant::from_earth_days(31.5)), "2 month 2, year 0, 12:00");
    }

    #[test]
    fn lunar_months_follow_the_new_moons() {
        let clock = earth();
        let synodic = Span::from_earth_days(29.53);
        let new_moon = Instant::from_earth_days(10.0);
        let c = Lunar { clock, moon: "Moon".into(), synodic, new_moon };
        check(&c, -3..3);
        let months = c.months(0);
        // Ten days before the first new moon, twelve lunations, and the rest cut by the year.
        assert_eq!(months.len(), 14);
        assert_eq!(months[0].name, "Moon 0");
        assert!((months[0].len.earth_days() - 10.0).abs() < 1e-9 && (months[1].len - synodic).seconds().abs() < 1e-6);
        assert!(months.iter().skip(1).all(|m| m.len <= synodic + Span::from_seconds(1e-6)));
        // Every month after month 0 begins at a new moon.
        for m in &months[1..] {
            let phase = ((clock.year_start(0) + m.start - new_moon) / synodic).rem_euclid(1.0);
            assert!(phase.min(1.0 - phase) < 1e-9);
        }
        // A year that begins at a new moon has no month 0.
        let aligned = Lunar { new_moon: Instant::EPOCH, ..c.clone() };
        assert_eq!(aligned.months(0)[0].name, "Moon 1");
        // A month longer than the year falls back to twelfths.
        let slow = Lunar { synodic: Span::from_earth_days(400.0), ..c };
        assert!(!slow.fits() && slow.months(0).len() == 12);
    }

    #[test]
    fn a_fantasy_calendar() {
        let mars = Clock::of_body(&Archetype::Mars.context(Key(1), "Mars"));
        let json = r#"{"name": "Red reckoning", "months": [
            {"name": "Dust", "days": 167}, {"name": "Frost", "days": 167}, {"name": "Thaw", "days": 167}, {"name": "Bloom"}],
            "week": ["Sol", "Phobia", "Deima"], "present_year": 212, "era": "AL"}"#;
        let def: CalendarDef = serde_json::from_str(json).unwrap();
        let c = Custom::new(mars, def.clone()).unwrap();
        check(&c, -2..2);
        let months = c.months(0);
        assert!((months[3].len / mars.day_len() - (mars.days_per_year() - 501.0)).abs() < 1e-9, "the last takes the rest");
        let t = Instant::EPOCH + mars.day_len() * 170.25;
        let d = c.date(t);
        assert_eq!((d.month, d.day), (1, 3));
        assert_eq!(c.format(t), "Deima, 4 Frost, 212 AL, 06:00");
        assert_eq!(c.year_name(-1), "211 AL");
        let full = CalendarDef { months: vec![MonthDef { name: "All".into(), days: 700.0 }, MonthDef { name: "None".into(), days: 0.0 }], ..def.clone() };
        assert!(Custom::new(mars, full).is_err());
        let fraction = CalendarDef { months: vec![MonthDef { name: "Half".into(), days: 30.5 }, MonthDef { name: "Rest".into(), days: 0.0 }], ..def };
        assert!(Custom::new(mars, fraction).is_err());
    }

    #[test]
    fn whole_day_months_keep_their_lengths() {
        // Month starts summed from lengths and day starts multiplied out differ in the last
        // bit; neither may move a day into the wrong month.
        let months = (0..11)
            .map(|i| MonthDef { name: format!("m{i}"), days: [31.0, 28.0, 31.0, 30.0, 31.0, 30.0, 31.0, 31.0, 30.0, 31.0, 30.0][i] })
            .chain(std::iter::once(MonthDef { name: "m11".into(), days: 0.0 }));
        let def = CalendarDef { name: "g".into(), months: months.collect(), week: vec![], present_year: 0, era: None };
        for clock in [Clock::new(Some(Span::from_earth_days(1.0)), Span::from_earth_days(365.0)), earth(), Clock::new(Some(Span::from_hours(24.6)), Span::from_earth_days(686.98))] {
            let def =
                if clock.days_per_year() > 400.0 { CalendarDef { months: def.months.iter().map(|m| MonthDef { days: 2.0 * m.days, ..m.clone() }).collect(), ..def.clone() } } else { def.clone() };
            let c = Custom::new(clock, def.clone()).unwrap();
            for y in -3..4 {
                let (a, b) = (clock.first_day_from(clock.year_start(y)), clock.first_day_from(clock.year_start(y + 1)));
                let mut counts = vec![0; 12];
                for k in a..b {
                    let d = c.date(clock.day_start(k));
                    assert_eq!(d.year, y);
                    counts[d.month] += 1;
                    assert_eq!(c.start_of(&d), Some(clock.day_start(k)));
                }
                assert!(counts[..11].iter().zip(&def.months).all(|(&n, m)| f64::from(n) == m.days), "year {y}: {counts:?}");
            }
        }
    }

    #[test]
    fn calendars_agree_on_the_year_and_day() {
        let clock = earth();
        let cals: Vec<Box<dyn Calendar>> = vec![
            Box::new(Twelfths(clock)),
            Box::new(Lunar { clock, moon: "Moon".into(), synodic: Span::from_earth_days(29.53), new_moon: Instant::from_earth_days(3.0) }),
            Box::new(
                Custom::new(
                    clock,
                    CalendarDef { name: "x".into(), months: vec![MonthDef { name: "a".into(), days: 100.0 }, MonthDef { name: "b".into(), days: 0.0 }], week: vec![], present_year: 0, era: None },
                )
                .unwrap(),
            ),
        ];
        for k in 0..50 {
            let t = Instant::from_earth_days(f64::from(k) * 37.1 - 600.0);
            let day_of_year = |c: &dyn Calendar| {
                let d = c.date(t);
                (d.year, ((c.start_of(&d).unwrap() - clock.year_start(d.year)) / clock.day_len()).round() as i64)
            };
            let first = day_of_year(cals[0].as_ref());
            assert!(cals.iter().all(|c| day_of_year(c.as_ref()) == first));
        }
    }
}
