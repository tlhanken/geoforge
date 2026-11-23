//! What is at a place: its climate and soil, the deposits likely within reach (and how many a
//! local map would hold), and the chance of each hazard over a stretch of days.
//!
//! Hazards come from [`geoforge_almanac::exposure`], which the daily weather shares.
//!
//! Deposits are drawn cell by cell, each a Poisson draw seeded by the world, the kind and the
//! cell, so any region's count is the sum of its cells' and every zoom agrees.

use crate::Result;
use geoforge_almanac::exposure::{self, HazardHere};
use geoforge_analysis::resources::{self, Deposit};
use geoforge_core::grid::{cells_within, Cell, Grid};
use geoforge_core::UnitVec3;
use geoforge_store::World;
use serde::Serialize;

#[derive(Debug, Clone, Serialize)]
pub struct DepositHere {
    pub kind: &'static str,
    /// Expected deposits in the area.
    pub expected: f64,
    /// Chance the area holds at least one.
    pub chance: f64,
    /// How many a local map of the area holds (seeded).
    pub drawn: u32,
    /// The settings at the centre and their rates per 10,000 km².
    pub settings: Vec<(&'static str, f64)>,
}

#[derive(Debug, Clone, Serialize)]
pub struct PlaceReport {
    pub lat: f64,
    pub lon: f64,
    pub radius_km: f64,
    pub area_km2: f64,
    pub elevation_m: f64,
    pub at_sea: bool,
    pub koppen: Option<&'static str>,
    pub soil: Option<&'static str>,
    pub fertility: Option<f64>,
    pub timber: Option<f64>,
    pub from_day: f64,
    pub days: f64,
    pub deposits: Vec<DepositHere>,
    /// The giant deposits within the area, and how far each is, km.
    pub giants: Vec<(geoforge_core::surface::GiantDeposit, f64)>,
    pub hazards: Vec<HazardHere>,
    /// Chance of a quake of magnitude 7 and 8 or more in the window, and the largest likely.
    pub quake_7: f64,
    pub quake_8: f64,
    pub largest_quake: f64,
    /// Droughts at the centre over the 20 years up to the window's start, and the drought
    /// index then; `None` at sea.
    pub droughts: Option<DroughtsHere>,
}

/// A place's droughts, without the month-by-month record.
#[derive(Debug, Clone, Serialize)]
pub struct DroughtsHere {
    pub years: usize,
    pub moderate_share: f64,
    pub severe_share: f64,
    pub longest_months: usize,
    pub dry_soil_share: Option<f64>,
    /// SPI-3 of the month the window starts in.
    pub spi_now: Option<f64>,
}

/// Report on the area within `radius_km` of (`lat`, `lon`) over `days` Earth days from
/// `from_day` (days from the present).
///
/// # Errors
/// If the world's layers are not generated yet.
pub fn report(world: &World, lat: f64, lon: f64, radius_km: f64, from_day: f64, days: f64) -> Result<PlaceReport> {
    let w = world.snapshot()?;
    let g = exposure::ground(&w)?;
    let centre = w.grid.locate(UnitVec3::from_lat_lon(lat, lon));
    let area: Vec<Cell> = cells_within(&w.grid, &[centre], radius_km / w.body.exterior.radius_km, |_| true);
    let area_km2: f64 = area.iter().map(|&c| g.area_km2(c as usize)).sum();
    let c = centre as usize;
    let sum = |field: &[f64]| area.iter().map(|&k| field[k as usize] * g.area_km2(k as usize) / 1e4).sum::<f64>();

    let deposit_seed = world.manifest.root_seed.child(b"deposits");
    let deposits = Deposit::ALL
        .into_iter()
        .map(|d| {
            let lambda = resources::prospectivity(&g, d);
            let kind_seed = deposit_seed.child(d.slug().as_bytes());
            let drawn = area.iter().map(|&k| resources::draw_count(lambda[k as usize] * g.area_km2(k as usize) / 1e4, kind_seed.child_u64(b"cell", u64::from(k)))).sum();
            let expected = sum(&lambda);
            DepositHere { kind: d.slug(), expected, chance: 1.0 - (-expected).exp(), drawn, settings: resources::settings_at(&g, d, c) }
        })
        .filter(|d| d.expected > 0.0)
        .collect();
    let near = exposure::hazards_in(world, &w, &g, lat, lon, radius_km, from_day, days, true)?;
    let here = UnitVec3::from_lat_lon(lat, lon);
    let giants = w
        .products
        .deposits()
        .map(|t| t.rows.values().map(|d| (d.clone(), here.angle_to(UnitVec3::from_lat_lon(d.lat_deg, d.lon_deg)) * w.body.exterior.radius_km)).filter(|(_, km)| *km <= radius_km).collect())
        .unwrap_or_default();
    Ok(PlaceReport {
        lat,
        lon,
        radius_km,
        area_km2,
        elevation_m: g.elevation[c],
        at_sea: g.sea[c],
        koppen: g.koppen[c].map(|k| k.code()),
        soil: resources::soil(&g, c).map(|s| s.name()),
        fertility: resources::soil(&g, c).map(|s| s.fertility()),
        timber: resources::timber(&g, c),
        from_day,
        days,
        deposits,
        giants,
        hazards: near.hazards,
        quake_7: near.quake_7,
        quake_8: near.quake_8,
        largest_quake: near.largest_quake,
        droughts: if g.sea[c] {
            None
        } else {
            let d = geoforge_almanac::drought::droughts(world, lat, lon, from_day, geoforge_almanac::year::DROUGHT_YEARS)?;
            let spi_now = d.months.last().and_then(|m| m.spi3);
            Some(DroughtsHere { years: d.years, moderate_share: d.moderate_share, severe_share: d.severe_share, longest_months: d.longest_months, dry_soil_share: d.dry_soil_share, spi_now })
        },
    })
}
