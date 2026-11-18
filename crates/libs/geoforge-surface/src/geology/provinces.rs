//! Province characteristics and constants, ported from geoforge v1 (Stage 2).
//!
//! Widths are **total** widths in km unless noted. Where v1 expanded a width on each side
//! of a boundary, the total here keeps v1's effective result, except collision orogens:
//! v1's 1000–2500 km *per side* (2000–5000 km total) is reduced to 1000–2500 km total,
//! which spans the Alps-to-Tibet range on Earth.

use geoforge_core::surface::{GeologicProvince as G, Province};

pub const LIP_PROBABILITY: f64 = 0.10;
pub const CONTINENTAL_FLOOD_BASALT_AREA_FRACTION: f64 = 0.10;
pub const OCEANIC_PLATEAU_AREA_FRACTION: f64 = 0.08;
/// Largest LIPs on Earth (Ontong Java, Siberian Traps) are a few million km².
pub const LIP_MAX_AREA_KM2: f64 = 4.0e6;

pub const OCEAN_TRENCH_WIDTH_KM: f64 = 50.0;
pub const ACCRETIONARY_WEDGE_BASE_WIDTH_KM: f64 = 100.0;
pub const ACCRETIONARY_WEDGE_WIDTH_FACTOR: f64 = 0.15;
pub const ACCRETIONARY_WEDGE_MAX_MULTIPLIER: f64 = 2.0;
pub const FOREARC_BASIN_BASE_WIDTH_KM: f64 = 100.0;
pub const FOREARC_BASIN_WIDTH_FACTOR: f64 = 0.15;
pub const FOREARC_BASIN_MAX_MULTIPLIER: f64 = 2.0;
pub const VOLCANIC_ARC_BASE_WIDTH_KM: f64 = 50.0;
pub const VOLCANIC_ARC_WIDTH_FACTOR: f64 = 0.15;
pub const VOLCANIC_ARC_MAX_MULTIPLIER: f64 = 2.0;
pub const BACKARC_BASIN_BASE_WIDTH_KM: f64 = 200.0;
pub const BACKARC_BASIN_WIDTH_FACTOR: f64 = 0.15;
pub const BACKARC_BASIN_MAX_MULTIPLIER: f64 = 2.0;
pub const BACKARC_BASIN_MIN_PLATE_AREA_KM2: f64 = 500_000.0;
pub const MIN_CONVERGENCE_RATE_CM_PER_YEAR: f64 = 2.0;

pub const COLLISION_BASE_WIDTH_KM: f64 = 1000.0;
pub const COLLISION_WIDTH_FACTOR: f64 = 0.15;
pub const COLLISION_MAX_MULTIPLIER: f64 = 2.5;

/// Rift zone half-width range. v1 expanded 250–400 km on each side (500–800 km total),
/// which turned short divergent stretches into discs; the East African rift system is
/// ~300–500 km across, so 150–250 km each side.
pub const RIFT_MIN_HALF_WIDTH_KM: f64 = 150.0;
pub const RIFT_HALF_WIDTH_RANGE_KM: f64 = 100.0;

pub const SPREADING_RATE_FAST: f64 = 10.0;
pub const SPREADING_RATE_MEDIUM: f64 = 5.0;
pub const SPREADING_RATE_SLOW: f64 = 2.0;
/// Ridge half-widths (v1 expanded these on each side).
pub const RIDGE_HALF_WIDTH_FAST_KM: f64 = 60.0;
pub const RIDGE_HALF_WIDTH_MEDIUM_KM: f64 = 80.0;
pub const RIDGE_HALF_WIDTH_SLOW_KM: f64 = 120.0;
pub const RIDGE_HALF_WIDTH_ULTRASLOW_KM: f64 = 150.0;
pub const FRACTURE_ZONE_HALF_WIDTH_KM: f64 = 100.0;

pub const SHIELD_AREA_FRACTION: f64 = 0.25;
pub const INTRACRATONIC_BASIN_MIN_AREA_KM2: f64 = 3_000_000.0;
pub const INTRACRATONIC_BASIN_PROBABILITY: f64 = 0.15;
/// v1 used 8 px (~180 km at 1800 px); Earth's intracratonic basins (Michigan, Williston,
/// Congo) are 300–1200 km across, so 300 km radius.
pub const INTRACRATONIC_BASIN_RADIUS_KM: f64 = 300.0;

pub const PALEO_OROGEN_MIN_AREA_KM2: f64 = 3_000_000.0;
pub const PALEO_OROGEN_COUNT_LARGE: usize = 2;
pub const PALEO_OROGEN_COUNT_MEDIUM: usize = 1;
pub const PALEO_OROGEN_MIN_WIDTH_KM: f64 = 200.0;
pub const PALEO_OROGEN_WIDTH_RANGE_KM: f64 = 200.0;
/// v1: endpoints at least 20 px apart (~450 km at 1800 px).
pub const PALEO_OROGEN_MIN_LENGTH_KM: f64 = 450.0;

pub const PASSIVE_MARGIN_DIVERGENT_WIDTH_KM: f64 = 350.0;
pub const PASSIVE_MARGIN_OTHER_WIDTH_KM: f64 = 130.0;

pub const HOTSPOT_PROBABILITY_LARGE_PLATE: f64 = 0.15;
pub const HOTSPOT_PROBABILITY_MEDIUM_PLATE: f64 = 0.10;
pub const HOTSPOT_MIN_PLATE_AREA_KM2: f64 = 2_000_000.0;
pub const HOTSPOT_LARGE_PLATE_THRESHOLD_KM2: f64 = 10_000_000.0;
pub const HOTSPOT_OCEANIC_MAX_LENGTH_KM: f64 = 2_400.0;
pub const HOTSPOT_CONTINENTAL_MAX_LENGTH_KM: f64 = 400.0;
/// Half-widths (v1 expanded these on each side).
pub const HOTSPOT_OCEANIC_HALF_WIDTH_KM: f64 = 50.0;
pub const HOTSPOT_CONTINENTAL_HALF_WIDTH_KM: f64 = 100.0;
/// v1: 20 px / 10 px from the plate edge (~450 / ~220 km at 1800 px).
pub const HOTSPOT_MIN_INTERIOR_DISTANCE_LARGE_KM: f64 = 450.0;
pub const HOTSPOT_MIN_INTERIOR_DISTANCE_MEDIUM_KM: f64 = 220.0;
pub const HOTSPOT_OCEANIC_MIN_AGE_MA: f64 = 10.0;
pub const HOTSPOT_OCEANIC_AGE_RANGE_MA: f64 = 18.0;
pub const HOTSPOT_CONTINENTAL_MIN_AGE_MA: f64 = 2.0;
pub const HOTSPOT_CONTINENTAL_AGE_RANGE_MA: f64 = 3.0;
pub const HOTSPOT_MIN_VISIBLE_LENGTH_KM: f64 = 200.0;

/// `base × min(1 + factor × (rate − 2 cm/yr), max)`: faster convergence, wider feature.
#[must_use]
pub fn rate_scaled(base: f64, factor: f64, max_mult: f64, rate_cm_yr: f64) -> f64 {
    let above = (rate_cm_yr - MIN_CONVERGENCE_RATE_CM_PER_YEAR).max(0.0);
    base * (1.0 + factor * above).min(max_mult)
}

fn p(kind: G, elevation_intensity: f64, roughness: f64, width_km: f64, intensity: f64, rate: f64) -> Province {
    Province { kind, elevation_intensity, roughness, width_km, intensity, convergence_rate_cm_yr: rate, boundary: None, age_ma: 0.0, anchor: None }
}

fn norm(rate: f64) -> f64 {
    (rate / 10.0).min(1.0)
}

#[must_use]
pub fn collision_orogen(rate: f64, width_km: f64) -> Province {
    p(G::CollisionOrogen, 1.0, 0.9, width_km, norm(rate), rate)
}
#[must_use]
pub fn accretionary_wedge(rate: f64, width_km: f64) -> Province {
    p(G::AccretionaryWedge, 0.2, 0.7, width_km, norm(rate), rate)
}
#[must_use]
pub fn continental_flood_basalt(area_km2: f64) -> Province {
    p(G::ContinentalFloodBasalt, 0.3, 0.2, area_km2.sqrt(), 0.8, 0.0)
}
#[must_use]
pub fn oceanic_plateau(area_km2: f64) -> Province {
    p(G::OceanicPlateau, -0.3, 0.3, area_km2.sqrt(), 0.7, 0.0)
}
#[must_use]
pub fn volcanic_arc(rate: f64, length_km: f64) -> Province {
    p(G::VolcanicArc, 0.8, 0.9, length_km, norm(rate), rate)
}
#[must_use]
pub fn forearc_basin(width_km: f64) -> Province {
    p(G::ForearcBasin, -0.2, 0.1, width_km, 0.3, 0.0)
}
#[must_use]
pub fn backarc_basin(width_km: f64) -> Province {
    p(G::BackarcBasin, -0.4, 0.2, width_km, 0.5, 0.0)
}
#[must_use]
pub fn craton(area_km2: f64) -> Province {
    p(G::Craton, 0.2, 0.1, area_km2.sqrt(), 0.0, 0.0)
}
#[must_use]
pub fn platform(area_km2: f64) -> Province {
    p(G::Platform, 0.15, 0.05, area_km2.sqrt(), 0.0, 0.0)
}
#[must_use]
pub fn intracratonic_basin(area_km2: f64) -> Province {
    p(G::IntracratonicBasin, 0.0, 0.05, area_km2.sqrt(), 0.0, 0.0)
}
#[must_use]
pub fn continental_rift(rate: f64, length_km: f64) -> Province {
    p(G::ContinentalRift, 0.3, 0.7, length_km, norm(rate), -rate)
}
#[must_use]
pub fn extended_crust(length_km: f64) -> Province {
    p(G::ExtendedCrust, 0.1, 0.5, length_km, 0.3, 0.0)
}
#[must_use]
pub fn mid_ocean_ridge(rate: f64, length_km: f64) -> Province {
    p(G::MidOceanRidge, -0.5, 0.7, length_km, norm(rate), -rate)
}
#[must_use]
pub fn abyssal_plain(area_km2: f64) -> Province {
    p(G::AbyssalPlain, -1.0, 0.05, area_km2.sqrt(), 0.0, 0.0)
}
#[must_use]
pub fn ocean_trench(rate: f64, length_km: f64) -> Province {
    p(G::OceanTrench, -1.2, 0.4, length_km, norm(rate), rate)
}
#[must_use]
pub fn oceanic_fracture_zone(length_km: f64) -> Province {
    p(G::OceanicFractureZone, -0.7, 0.6, length_km, 0.3, 0.0)
}
#[must_use]
pub fn oceanic_hotspot_track(length_km: f64) -> Province {
    p(G::OceanicHotspotTrack, -0.3, 0.9, length_km, 0.7, 0.0)
}
#[must_use]
pub fn continental_hotspot_track(length_km: f64) -> Province {
    p(G::ContinentalHotspotTrack, 0.4, 0.8, length_km, 0.8, 0.0)
}
#[must_use]
pub fn paleo_orogen(width_km: f64) -> Province {
    p(G::PaleoOrogen, 0.6, 0.6, width_km, 0.3, 0.0)
}
