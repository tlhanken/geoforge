//! The calendars a world's dates can be shown in: twelfths of its year, the lunations of one
//! of its moons (new moons found from the ephemeris), or its own (`calendar.json`).

use crate::ephemeris::{self, Calendar as Days};
use crate::{Error, Result};
use geoforge_core::body::{NodeKind, SolarSystem};
use geoforge_core::key::Key;
use geoforge_core::time::{Calendar, Clock, Custom, Instant, Lunar, Span, Twelfths};
use geoforge_store::World;
use std::str::FromStr;

/// Which calendar to show dates in.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub enum Choice {
    /// The world's own calendar if it has one, twelfths of the year otherwise.
    #[default]
    World,
    Twelfths,
    /// The lunations of the named moon, or of the heaviest.
    Lunar(Option<String>),
}

impl FromStr for Choice {
    type Err = String;
    fn from_str(s: &str) -> std::result::Result<Self, String> {
        match s.split_once(':') {
            None if s == "world" => Ok(Self::World),
            None if s == "twelfths" => Ok(Self::Twelfths),
            None if s == "lunar" => Ok(Self::Lunar(None)),
            Some(("lunar", moon)) => Ok(Self::Lunar(Some(moon.to_string()))),
            _ => Err(format!("unknown calendar {s}: world, twelfths, lunar or lunar:<moon>")),
        }
    }
}

/// The calendar `choice` names on `world`.
///
/// # Errors
/// If the world's `calendar.json` does not fit its year, or the world has no such moon.
pub fn of_world(world: &World, choice: &Choice) -> Result<Box<dyn Calendar>> {
    let clock = Days::new(&world.load_body()?, &[]).clock();
    Ok(match choice {
        Choice::Twelfths => Box::new(Twelfths(clock)),
        Choice::World => match world.load_calendar()? {
            Some(def) => Box::new(Custom::new(clock, def).map_err(Error::Invalid)?),
            None => Box::new(Twelfths(clock)),
        },
        Choice::Lunar(moon) => {
            let system = world.load_system()?.ok_or_else(|| Error::Invalid("the world has no system, so no moons".into()))?;
            Box::new(lunar(&system, world.manifest.focus, clock, moon.as_deref())?)
        }
    })
}

/// The lunar calendar of `world`'s moon `moon` (the heaviest if `None`) on `clock`.
///
/// # Errors
/// If the world has no moons, or none of that name.
pub fn lunar(system: &SolarSystem, world: Key, clock: Clock, moon: Option<&str>) -> Result<Lunar> {
    let moons: Vec<_> = system.children(world).into_iter().filter(|n| matches!(n.node, NodeKind::Body(_))).collect();
    let chosen = match moon {
        Some(name) => moons.iter().find(|n| n.name.eq_ignore_ascii_case(name)),
        None => moons.iter().max_by(|a, b| system.mass_of(a.key).total_cmp(&system.mass_of(b.key))),
    };
    let names = || moons.iter().map(|n| n.name.as_str()).collect::<Vec<_>>().join(", ");
    let m = chosen.ok_or_else(|| Error::Invalid(if moons.is_empty() { "the world has no moons".into() } else { format!("no moon {}: {}", moon.unwrap_or(""), names()) }))?;
    let period = ephemeris::period_days(system, m.key).ok_or_else(|| Error::Invalid(format!("{} has no orbit", m.name)))?;
    let synodic = Span::from_earth_days(1.0 / (1.0 / period - 1.0 / clock.year.earth_days()));
    Ok(Lunar { clock, moon: m.name.clone(), synodic, new_moon: new_moon(system, world, m.key, synodic, Instant::EPOCH) })
}

/// The first new moon of `moon` seen from `world` after `from`: where its lit share is least.
#[must_use]
pub fn new_moon(system: &SolarSystem, world: Key, moon: Key, synodic: Span, from: Instant) -> Instant {
    let lit = |t: Instant| {
        let placed = ephemeris::positions(system, t.earth_days());
        let at = |k| placed.iter().find(|p| p.key == k).map_or([0.0; 3], |p| p.position_au);
        ephemeris::lit_fraction(at(moon), at(world), [0.0; 3])
    };
    // The least of 96 samples over one month, then golden sections about it.
    let step = synodic / 96.0;
    let k = (0..96).min_by(|&a, &b| lit(from + step * f64::from(a)).total_cmp(&lit(from + step * f64::from(b)))).unwrap_or(0);
    let (mut a, mut b) = (from + step * (f64::from(k) - 1.0), from + step * (f64::from(k) + 1.0));
    let g = 0.5 * (5f64.sqrt() - 1.0);
    for _ in 0..40 {
        let (c, d) = (b - (b - a) * g, a + (b - a) * g);
        if lit(c) < lit(d) {
            b = d;
        } else {
            a = c;
        }
    }
    let t = a + (b - a) * 0.5;
    if t < from {
        t + synodic
    } else {
        t
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn choices_parse() {
        assert_eq!("world".parse(), Ok(Choice::World));
        assert_eq!("twelfths".parse(), Ok(Choice::Twelfths));
        assert_eq!("lunar".parse(), Ok(Choice::Lunar(None)));
        assert_eq!("lunar:Phobos".parse(), Ok(Choice::Lunar(Some("Phobos".into()))));
        assert!("gregorian".parse::<Choice>().is_err());
    }

    #[test]
    fn the_moons_months_follow_its_new_moons() {
        let (system, body) = geoforge_cosmos::sol::archetype(geoforge_core::Seed::new(1), geoforge_core::Archetype::Earth);
        let earth = system.by_name("Earth").unwrap().key;
        let clock = Clock::of_body(&body);
        let c = lunar(&system, earth, clock, None).unwrap();
        assert_eq!(c.moon, "Moon");
        assert!((c.synodic.earth_days() - 29.53).abs() < 0.01, "{:?}", c.synodic);
        // The present is the 2026 March equinox (20 March, 14:46 UTC); the next new moon was
        // 17 April, 11:52 UTC, 27.88 days on. Keplerian elements miss the Moon's perturbations
        // by up to half a day.
        assert!((c.new_moon.earth_days() - 27.88).abs() < 0.6, "{:?}", c.new_moon.earth_days());
        assert!(lunar(&system, earth, clock, Some("Phobos")).is_err());
        // Phobos's month is shorter than Mars's day: twelfths instead.
        let (mars_system, mars) = geoforge_cosmos::sol::archetype(geoforge_core::Seed::new(1), geoforge_core::Archetype::Mars);
        let key = mars_system.by_name("Mars").unwrap().key;
        let phobos = lunar(&mars_system, key, Clock::of_body(&mars), Some("phobos")).unwrap();
        assert!(!phobos.fits() && phobos.months(0).len() == 12);
    }
}
