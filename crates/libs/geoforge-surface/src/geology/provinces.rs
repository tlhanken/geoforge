//! Province characteristics and constants, ported from geoforge v1 (Stage 2).
//!
//! Widths are **total** widths in km unless noted. Where v1 expanded a width on each side
//! of a boundary, the total here keeps v1's effective result. Continental collisions are not
//! v1's: they are built asymmetrically (see [`Collision`]).

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

/// Ocean–continent ("Andean") margins. The trench lies at the foot of the continental
/// slope, the coast ~100 km inland of it, and the volcanic front ~200–300 km from the trench
/// (Chile, Peru, Cascadia, Sumatra): a narrow wedge under the slope, an onshore forearc (the
/// coastal range and central valley), then a volcanic cordillera ~150–250 km wide. Behind
/// it the upper plate is shortened, not stretched, so a retroarc foreland basin (the Chaco,
/// the Beni) takes the place of an extensional backarc. Widths grow with convergence rate as
/// the ocean–ocean ones do.
pub const ANDEAN_WEDGE_BASE_WIDTH_KM: f64 = 60.0;
pub const ANDEAN_FOREARC_BASE_WIDTH_KM: f64 = 70.0;
pub const ANDEAN_ARC_BASE_WIDTH_KM: f64 = 90.0;
pub const RETROARC_FORELAND_BASE_WIDTH_KM: f64 = 150.0;
/// The Altiplano–Puna: where fast convergence meets an old, heavy slab, the slab anchors
/// the trench while the continent drives into it, and the crust behind the arc thickens
/// into a plateau (Capitanio et al. 2011, Nature 480: the central Andes sit where the
/// Nazca slab is oldest). ~300–400 km wide at ~3.8 km, three quarters of Tibet's height.
/// The age threshold sits between the slab under Peru (~30–40 Ma, no plateau) and under
/// the Altiplano (~50 Ma). Very old slabs sink steeply and roll back instead, stretching
/// the upper plate (Molnar & Atwater 1978: the Japan Sea, the Marianas), so beyond
/// ~100 Ma there is no plateau.
pub const ANDEAN_PLATEAU_BASE_WIDTH_KM: f64 = 220.0;
pub const ANDEAN_PLATEAU_MIN_SLAB_AGE_MA: f64 = 50.0;
pub const ANDEAN_PLATEAU_MAX_SLAB_AGE_MA: f64 = 100.0;
pub const ANDEAN_PLATEAU_MIN_RATE_CM_PER_YEAR: f64 = 4.0;
pub const ANDEAN_PLATEAU_RELIEF: f64 = 0.75;

/// Slowest convergence that still builds a collision belt: the Alps rise at under 1 cm/yr.
pub const COLLISION_MIN_RATE_CM_PER_YEAR: f64 = 0.5;

/// The shape of a continent–continent collision, across the boundary.
///
/// One plate goes under the other: its upper crust is scraped off and stacked into a
/// steep **front range** just across the boundary (the Himalaya, the Zagros, the Alps), its
/// lithosphere slides beneath the upper plate and thickens it into a **plateau** behind the
/// range (Tibet, the Iranian plateau), and the load of the range bends the down-going plate
/// into a **foreland basin** in front (the Ganges plain, Mesopotamia and the Gulf).
///
/// Fast convergence stacks more crust: the range rises higher and the plateau grows
/// wider and higher. Calibrated on Earth at one degree, with PB2002's normal convergence
/// rates: the Alps (~0.5 cm/yr) have a 2 km range and no plateau, and the Himalaya
/// (~3.7 cm/yr) a 5 km range before Tibet, 5 km high and ~650 km wide. Rate is all this
/// knows; how long a collision has run matters too, so the Iranian plateau behind the
/// Zagros (~3 cm/yr, but half as old) comes out ~40% above its real 1.3 km. Peaks rise above
/// these means (the terrain detail adds ridges). Heights are for Earth's gravity; multiply
/// by [`crate::relief::relief_scale`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Collision {
    /// Width of the front range on the upper plate, km.
    pub range_km: f64,
    /// Width of the plateau behind it, km (0: none).
    pub plateau_km: f64,
    /// Width of the foreland basin on the lower plate, km.
    pub foreland_km: f64,
    /// Crest of the front range above the continental base, m.
    pub range_m: f64,
    /// Plateau surface above the continental base (itself 0.3–0.9 km), m.
    pub plateau_m: f64,
    /// Deepest point of the foreland basin below the continental base, m (sediment-filled).
    pub foreland_m: f64,
}

impl Collision {
    #[must_use]
    pub fn at_rate(rate_cm_yr: f64) -> Self {
        let r = rate_cm_yr.max(0.0);
        let ramp = |lo: f64, hi: f64| ((r - lo) / (hi - lo)).clamp(0.0, 1.0);
        Self {
            range_km: (150.0 + 30.0 * r).min(350.0),
            plateau_km: if r < 2.0 { 0.0 } else { (250.0 + 110.0 * r).min(1000.0) },
            foreland_km: (150.0 + 30.0 * r).min(350.0),
            range_m: 2000.0 + 3000.0 * ramp(0.5, 4.0),
            plateau_m: 5000.0 * smooth(ramp(2.0, 4.0)),
            foreland_m: 250.0,
        }
    }
}

fn smooth(t: f64) -> f64 {
    t * t * (3.0 - 2.0 * t)
}

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
/// Belt widths before swelling and pinching: the Appalachians are 150–500 km across with
/// their foothills, the Urals' ranges 50–150 km within a wider deformed zone.
pub const PALEO_OROGEN_MIN_WIDTH_KM: f64 = 200.0;
pub const PALEO_OROGEN_WIDTH_RANGE_KM: f64 = 200.0;
/// Where the belt's axis runs, as a share of the way from the craton's edge to the coast:
/// the younger, outer belt and an older one closer in (Grenville lies inside the
/// Appalachians).
pub const PALEO_OROGEN_OUTER_SHARE: (f64, f64) = (0.45, 0.65);
pub const PALEO_OROGEN_INNER_SHARE: (f64, f64) = (0.2, 0.32);
/// Accreted ground narrower than this beyond the craton holds no belt; wider than the
/// maximum counts no more when choosing a side.
pub const PALEO_OROGEN_MIN_ROOM_KM: f64 = 350.0;
pub const PALEO_OROGEN_MAX_ROOM_KM: f64 = 1500.0;
/// A coast within this of a plate boundary is an active margin.
pub const PALEO_OROGEN_ACTIVE_MARGIN_KM: f64 = 300.0;
/// Half the arc of the craton's outline a belt wraps, degrees (a 90–150° arc).
pub const PALEO_OROGEN_MIN_HALF_ARC_DEG: f64 = 45.0;
pub const PALEO_OROGEN_HALF_ARC_RANGE_DEG: f64 = 30.0;

pub const PASSIVE_MARGIN_DIVERGENT_WIDTH_KM: f64 = 350.0;
pub const PASSIVE_MARGIN_OTHER_WIDTH_KM: f64 = 130.0;

pub const HOTSPOT_PROBABILITY_LARGE_PLATE: f64 = 0.15;
pub const HOTSPOT_PROBABILITY_MEDIUM_PLATE: f64 = 0.10;
pub const HOTSPOT_MIN_PLATE_AREA_KM2: f64 = 2_000_000.0;
pub const HOTSPOT_LARGE_PLATE_THRESHOLD_KM2: f64 = 10_000_000.0;
pub const HOTSPOT_OCEANIC_MAX_LENGTH_KM: f64 = 2_400.0;
/// Yellowstone's track, the Snake River Plain, runs ~800 km.
pub const HOTSPOT_CONTINENTAL_MAX_LENGTH_KM: f64 = 800.0;
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
pub fn orogenic_plateau(rate: f64, width_km: f64) -> Province {
    p(G::OrogenicPlateau, 1.0, 0.45, width_km, norm(rate), rate)
}
/// A plateau behind an Andean arc. Its `elevation_intensity` is its height as a share of
/// a collision plateau's at the same rate.
#[must_use]
pub fn andean_plateau(rate: f64, width_km: f64) -> Province {
    p(G::OrogenicPlateau, ANDEAN_PLATEAU_RELIEF, 0.45, width_km, norm(rate), rate)
}
#[must_use]
pub fn foreland_basin(rate: f64, width_km: f64) -> Province {
    p(G::ForelandBasin, -0.1, 0.2, width_km, norm(rate), rate)
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
pub fn basement_uplift(rate: f64, width_km: f64) -> Province {
    p(G::BasementUplift, 0.7, 0.8, width_km, norm(rate), rate)
}

pub fn paleo_orogen(width_km: f64) -> Province {
    p(G::PaleoOrogen, 0.6, 0.6, width_km, 0.3, 0.0)
}

// Flat slabs: a buoyant ridge or plateau entering a trench under a continent lifts the slab
// flat beneath it, and the basement breaks into block ranges far inland (the Sierras
// Pampeanas today, the Laramide Rockies 80–50 Ma).
/// How far out from the trench a ridge or plateau on the down-going plate counts as arriving, km.
pub const FLAT_SLAB_ARRIVAL_KM: f64 = 400.0;
/// Half the along-trench length of a flat slab, km (the Pampean flat slab is ~600 km long).
pub const FLAT_SLAB_HALF_LENGTH_KM: f64 = 350.0;
/// Inland band of block uplifts over a flat slab, km from the trench.
pub const FLAT_SLAB_NEAR_KM: f64 = 450.0;
pub const FLAT_SLAB_FAR_KM: f64 = 1000.0;
/// A long continental margin with no ridge arriving now may carry a past flat slab's
/// uplifts, farther inland (the Rockies stand 1000–1600 km from the old trench).
pub const FOSSIL_FLAT_SLAB_MIN_MARGIN_KM: f64 = 2500.0;
pub const FOSSIL_FLAT_SLAB_PROBABILITY: f64 = 0.5;
pub const FOSSIL_FLAT_SLAB_HALF_LENGTH_KM: (f64, f64) = (400.0, 800.0);
pub const FOSSIL_FLAT_SLAB_NEAR_KM: f64 = 700.0;
pub const FOSSIL_FLAT_SLAB_FAR_KM: f64 = 1500.0;
pub const FOSSIL_FLAT_SLAB_AGE_MA: (f64, f64) = (50.0, 80.0);

#[cfg(test)]
mod tests {
    use super::Collision;

    #[test]
    fn faster_collisions_build_bigger() {
        let (alps, zagros, himalaya) = (Collision::at_rate(0.5), Collision::at_rate(3.0), Collision::at_rate(3.7));
        assert_eq!(alps.plateau_km, 0.0);
        assert!(alps.range_m < zagros.range_m && zagros.range_m <= himalaya.range_m);
        assert!(zagros.plateau_m < himalaya.plateau_m && zagros.plateau_km < himalaya.plateau_km);
        assert!((4500.0..=5000.0).contains(&himalaya.plateau_m));
    }
}
