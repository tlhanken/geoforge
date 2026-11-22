//! Natural hazards near a place: their yearly rates over an area, the chance of each over a
//! stretch of days, and the quakes. Tropical cyclones are counted from the storm field's
//! tracks over ten years — the storms the daily weather has — so the overview and the
//! forecast agree; the rest come from the hazard layer's rates, kept in the world's
//! derived cache. Rates are per year of Earth days; seasonal hazards are packed into their
//! half of the world's own year.

use crate::ephemeris::Calendar;
use crate::storms;
use crate::{Error, Result};
use geoforge_analysis::ground::Ground;
use geoforge_analysis::hazards::{self, Hazard, Measure, Season};
use geoforge_core::grid::{cells_within, Cell, Grid};
use geoforge_core::UnitVec3;
use geoforge_store::{cache, Snapshot, World};
use serde::Serialize;

#[derive(Debug, Clone, Serialize)]
pub struct HazardHere {
    pub hazard: &'static str,
    pub event: &'static str,
    /// Expected events a year in the area (for hazards counted in days, days a year at a spot
    /// in it, on average).
    pub per_year: f64,
    /// Expected events in the window, its season counted.
    pub expected: f64,
    /// Chance of at least one in the window.
    pub chance: f64,
}

/// Hazards over an area: rates a year, the chance in a window, and the quakes.
#[derive(Debug, Clone, Serialize)]
pub struct HazardsNear {
    pub hazards: Vec<HazardHere>,
    pub quake_7: f64,
    pub quake_8: f64,
    pub largest_quake: f64,
    /// Days a year with thunder at the centre.
    pub thunder_days: f64,
}

/// The per-cell ground of a snapshot.
///
/// # Errors
/// If the world's layers are not generated yet.
pub fn ground(w: &Snapshot) -> Result<Ground<'_>> {
    Ground::read(&w.products, &w.grid, w.body.exterior.radius_km, w.body.interior.has_surface_liquid_water()).ok_or_else(|| Error::Invalid("the world is not generated yet: run generate".into()))
}

/// Hazards within `radius_km` of (`lat`, `lon`) over `days` from `from_day`. Tropical
/// cyclones are counted from the storm field's tracks, the same storms the daily weather
/// has; the rest come from the hazard layer's rates.
#[allow(clippy::too_many_arguments)]
pub fn hazards_in(world: &World, w: &Snapshot, g: &Ground, lat: f64, lon: f64, radius_km: f64, from_day: f64, days: f64, count_cyclones: bool) -> Result<HazardsNear> {
    let (grid, body) = (&w.grid, &w.body);
    let centre = grid.locate(UnitVec3::from_lat_lon(lat, lon));
    let area: Vec<Cell> = cells_within(grid, &[centre], radius_km / body.exterior.radius_km, |_| true);
    let area_km2: f64 = area.iter().map(|&c| g.area_km2(c as usize)).sum::<f64>().max(1.0);
    let c = centre as usize;
    let sum = |field: &[f64]| area.iter().map(|&k| field[k as usize] * g.area_km2(k as usize) / 1e4).sum::<f64>();
    let setting = hazards::Setting::new(&w.products, hazards::aurora_of(body));
    let calendar = Calendar::new(body, &[]);
    let warm_share = season_share(&calendar, from_day, days, lat >= 0.0);
    // Rates are a year of Earth days; the seasons follow the world's own year.
    let window_years = days / EARTH_YEAR_DAYS;
    let cyclones_a_year = match crate::weather::atmosphere(world, &w.products, grid, body).ok().filter(|_| count_cyclones) {
        Some(atm) => {
            let storms = atm.storms(0.0, storms::CYCLONE_YEARS * calendar.year_days);
            let earth_years = storms::CYCLONE_YEARS * calendar.year_days / EARTH_YEAR_DAYS;
            storms::passes(&storms, storms::StormKind::TropicalCyclone, lat, lon, radius_km + storms::CYCLONE_REACH_KM, body.exterior.radius_km) as f64 / earth_years
        }
        None => 0.0,
    };
    let mut out = Vec::new();
    for h in Hazard::ALL {
        if g.sea[c] && !h.at_sea() {
            continue;
        }
        let per_year = match h {
            Hazard::Cyclones => cyclones_a_year,
            _ => {
                let rate = cache::field(world, &format!("hazard_{}", h.slug()), HAZARD_VERSION, || hazards::rate(g, &setting, h));
                // Days of fog or dust at a spot do not add up over the area: average them.
                if h.measure() == Measure::Days {
                    sum(&rate) * 1e4 / area_km2
                } else {
                    sum(&rate)
                }
            }
        };
        if per_year <= 0.0 {
            continue;
        }
        // A seasonal hazard's year of events falls in its half of the year.
        let in_season = match h.season() {
            Season::Any => window_years,
            Season::Warm => 2.0 * warm_share * window_years,
            Season::Cold => 2.0 * (1.0 - warm_share) * window_years,
        };
        let expected = per_year * in_season;
        out.push(HazardHere { hazard: h.slug(), event: h.event(), per_year, expected, chance: 1.0 - (-expected).exp() });
    }
    let rate6 = cache::field(world, "quake_rate6", HAZARD_VERSION, || setting.quakes(grid, body.exterior.radius_km).0.clone());
    let corner = cache::field(world, "quake_corner", HAZARD_VERSION, || setting.quakes(grid, body.exterior.radius_km).1.clone());
    let quake = |m: f64| 1.0 - (-area.iter().map(|&k| hazards::quake_rate_at_least(rate6[k as usize], corner[k as usize], m) * g.area_km2(k as usize) / 1e4).sum::<f64>() * window_years).exp();
    let thunder_days = cache::field(world, "thunder_days", HAZARD_VERSION, || hazards::lightning(g).iter().map(|l| l.map_or(0.0, hazards::Lightning::thunder_days)).collect())[c];
    Ok(HazardsNear { hazards: out, quake_7: quake(7.0), quake_8: quake(8.0), largest_quake: area.iter().map(|&k| corner[k as usize]).fold(0.0, f64::max), thunder_days })
}

/// The hazards near a place, for a forecast: rates a year within `radius_km`, without the
/// tropical cyclones (a forecast has the storms themselves).
///
/// # Errors
/// If the world's layers are not generated yet.
pub fn hazards_near(world: &World, lat: f64, lon: f64, radius_km: f64) -> Result<HazardsNear> {
    let w = world.snapshot()?;
    let g = ground(&w)?;
    let year = Calendar::new(&w.body, &[]).year_days;
    hazards_in(world, &w, &g, lat, lon, radius_km, 0.0, year, false)
}

/// Version of the cached hazard fields: bump when a hazard formula changes.
const HAZARD_VERSION: u32 = 1;

/// Hazard rates are per year of this many Earth days.
pub const EARTH_YEAR_DAYS: f64 = 365.25;

/// Share of the window in the warm half of the year for the hemisphere.
pub fn season_share(calendar: &Calendar, from_day: f64, days: f64, north: bool) -> f64 {
    if days <= 0.0 {
        return 0.5;
    }
    let year = calendar.year_days;
    let steps = 400;
    let warm = (0..steps)
        .filter(|&i| {
            let phase = ((from_day + days * (f64::from(i) + 0.5) / f64::from(steps)) / year).rem_euclid(1.0);
            (phase < 0.5) == north
        })
        .count();
    warm as f64 / f64::from(steps)
}

#[cfg(test)]
mod tests {
    use super::*;
    use geoforge_core::{Archetype, Key};

    #[test]
    fn seasons_split_the_year() {
        let c = Calendar::new(&Archetype::Earth.context(Key(1), "Earth"), &[]);
        let year = c.year_days;
        assert!((season_share(&c, 0.0, year, true) - 0.5).abs() < 0.01);
        assert!((season_share(&c, 0.0, year / 2.0, true) - 1.0).abs() < 0.01, "the first half is the northern summer");
        assert!(season_share(&c, 0.0, year / 2.0, false) < 0.01);
        assert!((season_share(&c, year / 4.0, year / 2.0, true) - 0.5).abs() < 0.01);
    }
}
