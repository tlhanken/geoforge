//! Geological constants for province generation
//!
//! This module contains all empirically-derived geological parameters used in province generation.
//! Values are based on Earth observations and are documented with real-world examples.

// ============================================================================
// Large Igneous Provinces
// ============================================================================

/// Probability of generating a Large Igneous Province per suitable plate
///
/// Earth has ~10-20 major LIPs over geological time, making them quite rare events.
pub const LIP_PROBABILITY: f64 = 0.10;

/// Fraction of plate area occupied by continental flood basalts
pub const CONTINENTAL_FLOOD_BASALT_AREA_FRACTION: f64 = 0.1;

/// Fraction of plate area occupied by oceanic plateaus
pub const OCEANIC_PLATEAU_AREA_FRACTION: f64 = 0.08;

// ============================================================================
// Subduction Zone Features
// ============================================================================

/// Ocean trench width in kilometers
///
/// Trenches are narrow features at subduction zones.
/// Example: Mariana Trench ~70 km wide (we use 50 km for visibility)
pub const OCEAN_TRENCH_WIDTH_KM: f64 = 50.0;

/// Accretionary wedge base width in kilometers
///
/// Sediments scraped off the subducting plate pile up on the overriding side.
/// Example: Barbados accretionary wedge ~100 km wide
pub const ACCRETIONARY_WEDGE_BASE_WIDTH_KM: f64 = 100.0;

/// Width scaling factor per cm/year above minimum convergence rate
///
/// Faster subduction scrapes off more sediment, creating wider wedges.
pub const ACCRETIONARY_WEDGE_WIDTH_FACTOR: f64 = 0.15;

/// Maximum width multiplier for accretionary wedges
pub const ACCRETIONARY_WEDGE_MAX_MULTIPLIER: f64 = 2.0;

/// Forearc basin base width in kilometers
///
/// Sedimentary basin between the trench and volcanic arc.
/// Example: Great Valley (California) ~100 km wide
pub const FOREARC_BASIN_BASE_WIDTH_KM: f64 = 100.0;

/// Width scaling factor per cm/year for forearc basins
pub const FOREARC_BASIN_WIDTH_FACTOR: f64 = 0.15;

/// Maximum width multiplier for forearc basins
pub const FOREARC_BASIN_MAX_MULTIPLIER: f64 = 2.0;

/// Volcanic arc base width in kilometers
///
/// Chain of volcanoes above the subducting slab.
/// Example: Cascade Range ~50 km wide
pub const VOLCANIC_ARC_BASE_WIDTH_KM: f64 = 50.0;

/// Width scaling factor per cm/year for volcanic arcs
pub const VOLCANIC_ARC_WIDTH_FACTOR: f64 = 0.15;

/// Maximum width multiplier for volcanic arcs
pub const VOLCANIC_ARC_MAX_MULTIPLIER: f64 = 2.0;

/// Backarc basin base width in kilometers
///
/// Extensional basin behind the volcanic arc.
/// Example: Sea of Japan ~200 km wide
pub const BACKARC_BASIN_BASE_WIDTH_KM: f64 = 200.0;

/// Width scaling factor per cm/year for backarc basins
pub const BACKARC_BASIN_WIDTH_FACTOR: f64 = 0.15;

/// Maximum width multiplier for backarc basins
pub const BACKARC_BASIN_MAX_MULTIPLIER: f64 = 2.0;

/// Minimum plate area (km²) required to generate backarc basins
///
/// Only large plates develop backarc basins.
pub const BACKARC_BASIN_MIN_PLATE_AREA_KM2: u64 = 500_000;

/// Minimum convergence rate (cm/year) for subduction features
pub const MIN_CONVERGENCE_RATE_CM_PER_YEAR: f64 = 2.0;

// ============================================================================
// Continental Rifts
// ============================================================================

/// Minimum width for continental rifts in kilometers
///
/// Active rifts are zones of extension and volcanism.
/// Example: East African Rift ~250 km wide
pub const RIFT_MIN_WIDTH_KM: f64 = 250.0;

/// Additional random width range for continental rifts
///
/// Combined with MIN gives 250-400 km range.
/// Example: Rio Grande Rift ~400 km wide
pub const RIFT_WIDTH_RANGE_KM: f64 = 150.0;

// ============================================================================
// Mid-Ocean Ridges
// ============================================================================

/// Spreading rate threshold for fast-spreading ridges (cm/year)
///
/// Fast-spreading ridges have smooth gentle rises without deep rift valleys.
/// Example: East Pacific Rise >10 cm/yr
pub const SPREADING_RATE_FAST: f64 = 10.0;

/// Spreading rate threshold for medium-spreading ridges (cm/year)
pub const SPREADING_RATE_MEDIUM: f64 = 5.0;

/// Spreading rate threshold for slow-spreading ridges (cm/year)
///
/// Slow-spreading ridges have deep central rift valleys.
/// Example: Mid-Atlantic Ridge 2-5 cm/yr
pub const SPREADING_RATE_SLOW: f64 = 2.0;

/// Width of fast-spreading ridges in kilometers
///
/// Smooth gentle rise, no deep rift valley.
/// Example: East Pacific Rise ~60 km wide
pub const RIDGE_WIDTH_FAST_KM: f64 = 60.0;

/// Width of medium-spreading ridges in kilometers
pub const RIDGE_WIDTH_MEDIUM_KM: f64 = 80.0;

/// Width of slow-spreading ridges in kilometers
///
/// Deep central rift valley 1-2 km deep.
/// Example: Mid-Atlantic Ridge ~120 km wide
pub const RIDGE_WIDTH_SLOW_KM: f64 = 120.0;

/// Width of ultra-slow spreading ridges in kilometers
///
/// Highly irregular, deepest rift valleys.
/// Example: Gakkel Ridge ~150 km wide
pub const RIDGE_WIDTH_ULTRASLOW_KM: f64 = 150.0;

/// Width of oceanic fracture zones in kilometers
///
/// Transform faults offsetting ridge segments.
/// Example: Romanche Fracture Zone ~100 km wide
pub const FRACTURE_ZONE_WIDTH_KM: f64 = 100.0;

// ============================================================================
// Stable Continental Regions
// ============================================================================

/// Fraction of continental plate area occupied by exposed shield (craton core)
///
/// Shields are exposed Precambrian basement rock.
/// Example: Canadian Shield occupies ~25% of North America
pub const SHIELD_AREA_FRACTION: f64 = 0.25;

/// Minimum continental plate area (km²) for generating intracratonic basins
///
/// Basins only form in large, stable cratons.
/// Example: Michigan Basin, Illinois Basin
pub const INTRACRATONIC_BASIN_MIN_AREA_KM2: u64 = 3_000_000;

/// Probability of generating an intracratonic basin on large plates
pub const INTRACRATONIC_BASIN_PROBABILITY: f64 = 0.15;

/// Minimum radius for intracratonic basins in pixels
pub const INTRACRATONIC_BASIN_MIN_RADIUS_PX: f64 = 8.0;

/// Minimum pixel count for a valid intracratonic basin
pub const INTRACRATONIC_BASIN_MIN_PIXELS: usize = 10;

/// Area estimate per basin pixel in km²
pub const INTRACRATONIC_BASIN_AREA_PER_PIXEL_KM2: f64 = 10_000.0;

/// Minimum continental plate area (km²) for generating paleo-orogens
///
/// Ancient mountain belts only visible in large continents.
pub const PALEO_OROGEN_MIN_AREA_KM2: u64 = 3_000_000;

/// Number of paleo-orogen belts for large continents (>10M km²)
pub const PALEO_OROGEN_COUNT_LARGE: usize = 2;

/// Number of paleo-orogen belts for medium continents
pub const PALEO_OROGEN_COUNT_MEDIUM: usize = 1;

/// Minimum width for paleo-orogens in kilometers
pub const PALEO_OROGEN_MIN_WIDTH_KM: f64 = 200.0;

/// Additional random width range for paleo-orogens
pub const PALEO_OROGEN_WIDTH_RANGE_KM: f64 = 200.0;

/// Minimum distance squared for paleo-orogen endpoints (pixels)
pub const PALEO_OROGEN_MIN_DISTANCE_SQ: i32 = 400;

/// Minimum pixel count for a valid paleo-orogen
pub const PALEO_OROGEN_MIN_PIXELS: usize = 20;

/// Width of passive continental margins (divergent boundaries) in kilometers
///
/// Wide margins at old rift zones.
/// Example: US Atlantic coast ~350 km
pub const PASSIVE_MARGIN_DIVERGENT_WIDTH_KM: f64 = 350.0;

/// Width of passive continental margins (convergent/transform) in kilometers
pub const PASSIVE_MARGIN_OTHER_WIDTH_KM: f64 = 130.0;

// ============================================================================
// Hotspot Tracks
// ============================================================================

/// Probability of generating a hotspot track on large plates (>10M km²)
///
/// Pacific-sized plates are more likely to have hotspots.
pub const HOTSPOT_PROBABILITY_LARGE_PLATE: f64 = 0.15;

/// Probability of generating a hotspot track on medium plates
pub const HOTSPOT_PROBABILITY_MEDIUM_PLATE: f64 = 0.10;

/// Minimum plate area (km²) required for hotspot generation
pub const HOTSPOT_MIN_PLATE_AREA_KM2: u64 = 2_000_000;

/// Plate area threshold (km²) for "large" plate classification
pub const HOTSPOT_LARGE_PLATE_THRESHOLD_KM2: u64 = 10_000_000;

/// Maximum visible length of oceanic hotspot tracks in kilometers
///
/// Only young, geologically active portions are terrain-significant.
/// Example: Hawaiian Islands + atolls (0-28 Ma) ~2,400 km
/// (Total Hawaiian-Emperor chain is 6,200 km, but old seamounts are eroded flat)
pub const HOTSPOT_OCEANIC_MAX_LENGTH_KM: f64 = 2_400.0;

/// Maximum visible length of continental hotspot tracks in kilometers
///
/// Continental hotspots have shorter visible tracks.
/// Example: Yellowstone calderas (0-5 Ma) ~400 km
/// (Total Yellowstone track is 800 km, but old calderas are buried)
pub const HOTSPOT_CONTINENTAL_MAX_LENGTH_KM: f64 = 400.0;

/// Width of oceanic hotspot tracks in kilometers
///
/// Narrow seamount chains.
pub const HOTSPOT_OCEANIC_WIDTH_KM: f64 = 50.0;

/// Width of continental hotspot tracks in kilometers
///
/// Wider volcanic fields on land.
pub const HOTSPOT_CONTINENTAL_WIDTH_KM: f64 = 100.0;

/// Minimum distance from plate boundary for hotspot generation (large plates)
///
/// Hotspots form in deep plate interiors, far from edges.
/// Large plates: 400+ km from edge (20 pixels at 20 km/pixel)
pub const HOTSPOT_MIN_INTERIOR_DISTANCE_LARGE_PX: usize = 20;

/// Minimum distance from plate boundary for hotspot generation (medium plates)
///
/// Medium plates: 200+ km from edge (10 pixels at 20 km/pixel)
pub const HOTSPOT_MIN_INTERIOR_DISTANCE_MEDIUM_PX: usize = 10;

/// Minimum age for oceanic hotspot tracks in Ma
pub const HOTSPOT_OCEANIC_MIN_AGE_MA: f64 = 10.0;

/// Additional random age range for oceanic hotspot tracks in Ma
///
/// Combined with MIN gives 10-28 Ma range.
pub const HOTSPOT_OCEANIC_AGE_RANGE_MA: f64 = 18.0;

/// Minimum age for continental hotspot tracks in Ma
pub const HOTSPOT_CONTINENTAL_MIN_AGE_MA: f64 = 2.0;

/// Additional random age range for continental hotspot tracks in Ma
///
/// Combined with MIN gives 2-5 Ma range.
pub const HOTSPOT_CONTINENTAL_AGE_RANGE_MA: f64 = 3.0;

/// Minimum visible hotspot track length in kilometers
///
/// Tracks shorter than this are filtered out as too small to be terrain-significant.
pub const HOTSPOT_MIN_VISIBLE_LENGTH_KM: f64 = 200.0;

/// Conversion factor from cm/yr × Ma to km
///
/// Units: cm/yr × Ma × (1,000,000 yr/Ma) ÷ (100,000 cm/km) = km
pub const HOTSPOT_LENGTH_CONVERSION_FACTOR: f64 = 10.0;

// ============================================================================
// Sample/Interior Detection
// ============================================================================

/// Fraction of plate interior to sample for LIP generation
pub const LIP_INTERIOR_SAMPLE_FRACTION: f64 = 0.1;
