#![allow(clippy::too_many_arguments)]
//! Geological province generator (Stage 2)
//!
//! This module implements the main geological province generator, coordinating the complete
//! Stage 2 pipeline to generate 18 distinct province types from tectonic plate data.
//! Provinces are layered from oldest/deepest to youngest/highest to create realistic
//! geological hierarchies.
//!
//! # Pipeline Order
//!
//! 1. **Foundation Layers** (Ancient/Stable):
//!    - Oceanic base layer (abyssal plains on all oceanic plates)
//!    - Stable continental regions (cratons, shields, platforms)
//!    - Passive margins (extended crust at old boundaries)
//!
//! 2. **Active Features** (Overlay on foundation):
//!    - Collision orogens (continent-continent mountain belts: Himalayas, Alps)
//!    - Subduction zone systems (trench → accretionary wedge → forearc → volcanic arc → backarc)
//!    - Oceanic features (mid-ocean ridges, fracture zones)
//!    - Active continental rifts (volcanic rift zones)
//!    - Large Igneous Provinces (flood basalts, oceanic plateaus)
//!
//! 3. **Final Overlays**:
//!    - Hotspot tracks (linear volcanic chains, generated last to sit on top)
//!
//! # Province Types (18 total - IntracratonicBasin not implemented, ContinentalRift deferred)
//!
//! **Collision Orogens (1)**: CollisionOrogen (continental-continental convergence)
//! **Subduction Systems (5)**: OceanTrench, AccretionaryWedge, ForearcBasin, VolcanicArc, BackarcBasin
//!   - Applies to BOTH oceanic-oceanic (island arcs) AND oceanic-continental (Andes-style) convergence
//!   - Same province types, different elevations (determined in Stage 3 based on crust type)
//! **LIPs (3)**: ContinentalFloodBasalt, OceanicPlateau, ContinentalHotspotTrack*
//! **Stable (3)**: Craton, Platform, ExtendedCrust
//! **Oceanic (4)**: AbyssalPlain, MidOceanRidge, OceanicFractureZone, OceanicHotspotTrack*
//! **Deferred**: ContinentalRift (will be implemented later)
//!
//! *Hotspot tracks are generated separately as final overlays

use crate::geology::constants as geo_const;
use crate::geology::provinces::{GeologicProvince, ProvinceCharacteristics, ProvinceRegion};
use crate::geology::orogenic::{OrogenicBeltGenerator, OrogenicConfig};
use crate::map::terrain::TerrainMap;
use crate::map::spherical::PlanetaryParams;
use crate::tectonics::boundary_analysis::BoundarySegment;
use crate::tectonics::plates::{PlateInteraction, PlateStats, PlateType, PlateSeed};
use std::collections::{HashMap, HashSet};
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;

/// Index mapping plate IDs to their pixel coordinates
///
/// Built once at the start of generation to avoid O(plates × width × height) scanning.
type PlatePixelIndex = HashMap<u16, Vec<(usize, usize)>>;

/// Configuration for geological province generation
#[derive(Debug, Clone)]
pub struct GeologyConfig {
    /// Configuration for orogenic belts
    pub orogenic_config: OrogenicConfig,

    /// Probability of generating Large Igneous Provinces (0.0-1.0)
    /// Earth has ~10-20 major LIPs, so quite rare
    pub lip_probability: f64,

    /// Enable generation of oceanic provinces (ridges, trenches, abyssal plains)
    pub generate_oceanic: bool,

    /// Enable generation of stable continental regions (cratons, platforms)
    pub generate_stable_regions: bool,

    /// Enable generation of arc and basin systems
    pub generate_arc_systems: bool,

    /// Enable generation of extensional zones
    pub generate_extensional: bool,
}

impl Default for GeologyConfig {
    fn default() -> Self {
        Self {
            orogenic_config: OrogenicConfig::default(),
            lip_probability: geo_const::LIP_PROBABILITY,
            generate_oceanic: true,
            generate_stable_regions: true,
            generate_arc_systems: true,
            generate_extensional: true,
        }
    }
}

/// Geological province generator
pub struct GeologyGenerator {
    config: GeologyConfig,
    seed: u64,
    /// Planetary parameters (for map scale calculations and other properties)
    planetary_params: PlanetaryParams,
}

impl GeologyGenerator {
    pub fn new(config: GeologyConfig, seed: u64, planetary_params: PlanetaryParams) -> Self {
        Self { config, seed, planetary_params }
    }

    /// Generate all geological provinces from tectonic data
    ///
    /// Executes the complete Stage 2 pipeline, generating 20 distinct province types
    /// in a layered approach from ancient/stable foundations to active overlays.
    ///
    /// Returns a vector of province regions sorted by generation order (oldest first).
    pub fn generate_all_provinces(
        &self,
        boundaries: &[BoundarySegment],
        plate_stats: &HashMap<u16, PlateStats>,
        _plate_seeds: &[PlateSeed],  // Currently unused, plate motion stored in PlateStats
        plate_map: &TerrainMap<u16>,
    ) -> Vec<ProvinceRegion> {
        let mut regions = Vec::new();
        let mut rng = StdRng::seed_from_u64(self.seed);

        // Build plate pixel index once for O(width × height) instead of O(plates × width × height)
        let plate_index = self.build_plate_pixel_index(plate_map);

        // ========== FOUNDATION LAYERS (Ancient/Stable) ==========

        // Stage 2.6a: Oceanic Base Layer - Fill all oceanic plates with abyssal plains
        // (Oldest oceanic crust, foundation for all oceanic features)
        if self.config.generate_oceanic {
            let oceanic_base = self.generate_oceanic_base_layer(&plate_index, plate_stats);
            regions.extend(oceanic_base);
        }

        // Stage 2.4: Stable Continental Regions - Shields, platforms, intracratonic basins
        // (Ancient Precambrian basement = Cratons, 1.5-4 Ga old, foundation for all continental features)
        // Shield = exposed craton core, Platform = sedimentary-covered craton
        if self.config.generate_stable_regions {
            let stable = self.generate_stable_continental_base(plate_stats, &plate_index, &mut rng);
            regions.extend(stable);
        }

        // Stage 2.5a: Extended Crust - Passive continental margins
        // (Thinned/stretched crust at continental plate edges, wider at old divergent boundaries)
        if self.config.generate_extensional {
            let passive_margins = self.generate_passive_margins(boundaries, plate_stats, plate_map);
            regions.extend(passive_margins);
        }

        // Stage 2.4b: Paleo-orogens (Ancient mountain belts)
        // Linear features that cross-cut the stable continental cores
        if self.config.generate_stable_regions {
            let paleo_orogens = self.generate_paleo_orogens(plate_stats, &plate_index, plate_map, &mut rng);
            regions.extend(paleo_orogens);
        }

        // ========== ACTIVE FEATURES (Overlay on foundation) ==========

        // Stage 2.1: Orogenic Belts - Mountain building from convergent boundaries
        // (Compress/fold edges of cratons and platforms)
        let orogenic_gen = OrogenicBeltGenerator::new(
            self.config.orogenic_config.clone(),
            self.planetary_params.clone()
        );
        let orogens = orogenic_gen.generate_orogens(boundaries, plate_stats, plate_map);
        regions.extend(orogens);

        // Stage 2.3: Arc and Basin Systems - Oceanic-oceanic convergence
        // (Trenches, volcanic arcs, forearc/backarc basins overlay oceanic base)
        if self.config.generate_arc_systems {
            let arcs = self.generate_arc_systems(boundaries, plate_stats, plate_map, &mut rng);
            regions.extend(arcs);
        }

        // Stage 2.6b: Oceanic Overlays - Ridges and fractures
        // (Active spreading and transform features on oceanic base)
        if self.config.generate_oceanic {
            let oceanic_overlays = self.generate_oceanic_overlays(boundaries, plate_map, plate_stats);
            regions.extend(oceanic_overlays);
        }

        // Stage 2.5b: Active Continental Rifts - Volcanic rift zones (East African Rift style)
        // (Active rifting with flood basalts and volcanism - igneous features)
        if self.config.generate_extensional {
            let active_rifts = self.generate_active_continental_rifts(boundaries, plate_stats, plate_map, &mut rng);
            regions.extend(active_rifts);
        }

        // Stage 2.2: Large Igneous Provinces - Rare volcanic events
        // (Flood basalts, hotspots overlay everything)
        let lips = self.generate_large_igneous_provinces(plate_stats, &plate_index, &regions, &mut rng);
        regions.extend(lips);

        // FINAL OVERLAY: Hotspot tracks - Linear volcanic chains on top of all other features
        // (These should sit on top of oceanic base, abyssal plains, and any other provinces)
        // Generated LAST so they're visible on top of everything else
        if self.config.generate_oceanic {
            let hotspot_tracks = self.generate_hotspot_tracks(plate_stats, plate_map, &plate_index, &mut rng);
            regions.extend(hotspot_tracks);
        }

        regions
    }

    /// Generate arc and basin systems from oceanic-oceanic convergence
    ///
    /// Creates a complete subduction zone system:
    /// [Subducting Plate] → [Trench] → [Forearc Basin] → [Volcanic Arc] → [Backarc Basin] → [Overriding Plate]
    fn generate_arc_systems(
        &self,
        boundaries: &[BoundarySegment],
        plate_stats: &HashMap<u16, PlateStats>,
        plate_map: &TerrainMap<u16>,
        rng: &mut StdRng,
    ) -> Vec<ProvinceRegion> {
        let mut regions = Vec::new();

        for (idx, boundary) in boundaries.iter().enumerate() {
            if boundary.interaction_type != PlateInteraction::Convergent {
                continue;
            }

            // Check for ANY subduction (oceanic-oceanic OR oceanic-continental)
            if let (Some(stats_a), Some(stats_b)) = (plate_stats.get(&boundary.plate_a), plate_stats.get(&boundary.plate_b)) {
                // Determine if this is a subduction zone and identify plates
                let (is_subduction, subducting_plate, overriding_plate) = match (stats_a.plate_type, stats_b.plate_type) {
                    // Oceanic-Oceanic: Older plate (lower ID) subducts
                    (PlateType::Oceanic, PlateType::Oceanic) => {
                        if boundary.plate_a < boundary.plate_b {
                            (true, boundary.plate_a, boundary.plate_b)
                        } else {
                            (true, boundary.plate_b, boundary.plate_a)
                        }
                    },
                    // Oceanic-Continental: Oceanic plate always subducts (denser)
                    (PlateType::Oceanic, PlateType::Continental) => (true, boundary.plate_a, boundary.plate_b),
                    (PlateType::Continental, PlateType::Oceanic) => (true, boundary.plate_b, boundary.plate_a),
                    // Continental-Continental: Handled by OrogenicBeltGenerator
                    _ => (false, 0, 0),
                };

                if is_subduction {

                    // Generate complete subduction zone transect sequentially:
                    // Oceanic plate → Trench → Accretionary Wedge → Forearc → Volcanic Arc → Backarc → Continental plate
                    self.create_ocean_trench(boundary, subducting_plate, idx, plate_map, &mut regions);

                    let wedge_width_km = self.create_accretionary_wedge(boundary, overriding_plate, idx, plate_map, rng, &mut regions);
                    let (forearc_offset_km, forearc_width_km) = self.create_forearc_basin(
                        boundary, overriding_plate, wedge_width_km, idx, plate_map, rng, &mut regions
                    );

                    let (arc_offset_km, arc_width_km) = self.create_volcanic_arc(
                        boundary, overriding_plate, forearc_offset_km, forearc_width_km, idx, plate_map, rng, &mut regions
                    );

                    self.create_backarc_basin(
                        boundary, overriding_plate, idx, arc_offset_km, arc_width_km,
                        stats_a, stats_b, plate_map, rng, &mut regions
                    );
                }
            }
        }

        regions
    }

    /// Create ocean trench at convergent boundary (on subducting plate side only)
    ///
    /// Ocean trenches are narrow features on the subducting plate
    /// Example: Mariana Trench ~70 km wide
    fn create_ocean_trench(
        &self,
        boundary: &BoundarySegment,
        subducting_plate: u16,
        boundary_idx: usize,
        plate_map: &TerrainMap<u16>,
        regions: &mut Vec<ProvinceRegion>,
    ) {
        let chars = ProvinceCharacteristics::ocean_trench(
            boundary.relative_velocity,
            boundary.length_km,
        );

        // Trenches: narrow feature right at the boundary, extending slightly onto subducting plate
        // Expand minimally to make visible on map
        let trench_width_km = geo_const::OCEAN_TRENCH_WIDTH_KM;

        let expanded = self.expand_boundary_spherical(&boundary.pixels, trench_width_km / 2.0, plate_map);

        // Filter to only subducting plate pixels
        let pixels: Vec<(usize, usize)> = expanded.iter()
            .filter(|&&(x, y)| {
                let plate_id = plate_map.data[y * plate_map.width + x];
                plate_id == subducting_plate
            })
            .copied()
            .collect();

        regions.push(ProvinceRegion::new(
            pixels,
            chars,
            Some(boundary_idx)
        ));
    }

    /// Create accretionary wedge (on overriding plate side, adjacent to trench)
    ///
    /// Accretionary wedges are sediments scraped off the subducting plate
    /// and piled up on the overriding plate side. Width: 100-200 km
    /// Example: Barbados accretionary wedge ~100 km wide
    ///
    /// Returns the width in km for positioning the next feature
    fn create_accretionary_wedge(
        &self,
        boundary: &BoundarySegment,
        overriding_plate: u16,
        boundary_idx: usize,
        plate_map: &TerrainMap<u16>,
        _rng: &mut StdRng,
        regions: &mut Vec<ProvinceRegion>,
    ) -> f64 {
        // Dynamic width based on convergence rate (faster subduction = more sediment scraped off)
        let base_width_km = geo_const::ACCRETIONARY_WEDGE_BASE_WIDTH_KM;
        let rate_above_min = (boundary.relative_velocity - geo_const::MIN_CONVERGENCE_RATE_CM_PER_YEAR).max(0.0);
        let multiplier = 1.0 + (geo_const::ACCRETIONARY_WEDGE_WIDTH_FACTOR * rate_above_min);
        let clamped_multiplier = multiplier.min(geo_const::ACCRETIONARY_WEDGE_MAX_MULTIPLIER);
        let width_km = base_width_km * clamped_multiplier;

        // Expand onto overriding plate (opposite side from trench)
        let pixels = self.expand_boundary_toward_plate_spherical(
            &boundary.pixels,
            overriding_plate,
            width_km,
            plate_map
        );

        let chars = ProvinceCharacteristics::accretionary_wedge(
            boundary.relative_velocity,
            width_km
        );
        regions.push(ProvinceRegion::new(pixels, chars, Some(boundary_idx)));

        width_km
    }

    /// Create forearc basin (behind accretionary wedge)
    ///
    /// Starts after the accretionary wedge and extends inland
    ///
    /// Returns (offset_km, width_km) for positioning the next feature
    fn create_forearc_basin(
        &self,
        boundary: &BoundarySegment,
        overriding_plate: u16,
        wedge_width_km: f64,
        boundary_idx: usize,
        plate_map: &TerrainMap<u16>,
        _rng: &mut StdRng,
        regions: &mut Vec<ProvinceRegion>,
    ) -> (f64, f64) {
        let offset_km = wedge_width_km; // Start immediately after wedge ends (no gap)

        // Dynamic width based on convergence rate (faster subduction = more deformation/subsidence)
        let base_width_km = geo_const::FOREARC_BASIN_BASE_WIDTH_KM;
        let rate_above_min = (boundary.relative_velocity - geo_const::MIN_CONVERGENCE_RATE_CM_PER_YEAR).max(0.0);
        let multiplier = 1.0 + (geo_const::FOREARC_BASIN_WIDTH_FACTOR * rate_above_min);
        let clamped_multiplier = multiplier.min(geo_const::FOREARC_BASIN_MAX_MULTIPLIER);
        let width_km = base_width_km * clamped_multiplier;

        // Expand unidirectionally: from offset to offset+width
        // This creates a "ring" that doesn't overlap with previous features
        let outer_edge = self.expand_boundary_toward_plate_spherical(
            &boundary.pixels,
            overriding_plate,
            offset_km + width_km,
            plate_map
        );

        let inner_edge = self.expand_boundary_toward_plate_spherical(
            &boundary.pixels,
            overriding_plate,
            offset_km,
            plate_map
        );

        // Forearc basin = outer edge minus inner edge
        let inner_set: std::collections::HashSet<_> = inner_edge.iter().copied().collect();
        let pixels: Vec<(usize, usize)> = outer_edge.iter()
            .filter(|p| !inner_set.contains(p))
            .copied()
            .collect();

        let chars = ProvinceCharacteristics::forearc_basin(boundary.length_km);
        regions.push(ProvinceRegion::new(pixels, chars, Some(boundary_idx)));

        (offset_km, width_km)
    }

    /// Create volcanic arc (behind forearc basin)
    ///
    /// Returns (arc_offset_km, arc_width_km) for backarc positioning
    fn create_volcanic_arc(
        &self,
        boundary: &BoundarySegment,
        overriding_plate: u16,
        forearc_offset_km: f64,
        forearc_width_km: f64,
        boundary_idx: usize,
        plate_map: &TerrainMap<u16>,
        _rng: &mut StdRng,
        regions: &mut Vec<ProvinceRegion>,
    ) -> (f64, f64) {
        let offset_km = forearc_offset_km + forearc_width_km; // Start immediately after forearc ends (no gap)

        // Dynamic width based on convergence rate (faster subduction = more vigorous magmatism)
        let base_width_km = geo_const::VOLCANIC_ARC_BASE_WIDTH_KM;
        let rate_above_min = (boundary.relative_velocity - geo_const::MIN_CONVERGENCE_RATE_CM_PER_YEAR).max(0.0);
        let multiplier = 1.0 + (geo_const::VOLCANIC_ARC_WIDTH_FACTOR * rate_above_min);
        let clamped_multiplier = multiplier.min(geo_const::VOLCANIC_ARC_MAX_MULTIPLIER);
        let width_km = base_width_km * clamped_multiplier;

        // Expand unidirectionally: from offset to offset+width
        // This creates a "ring" that doesn't overlap with previous features
        let outer_edge = self.expand_boundary_toward_plate_spherical(
            &boundary.pixels,
            overriding_plate,
            offset_km + width_km,
            plate_map
        );

        let inner_edge = self.expand_boundary_toward_plate_spherical(
            &boundary.pixels,
            overriding_plate,
            offset_km,
            plate_map
        );

        // Volcanic arc = outer edge minus inner edge
        let inner_set: std::collections::HashSet<_> = inner_edge.iter().copied().collect();
        let arc_pixels: Vec<(usize, usize)> = outer_edge.iter()
            .filter(|p| !inner_set.contains(p))
            .copied()
            .collect();

        let chars = ProvinceCharacteristics::volcanic_arc(
            boundary.relative_velocity,
            boundary.length_km,
        );
        regions.push(ProvinceRegion::new(arc_pixels, chars, Some(boundary_idx)));

        (offset_km, width_km)
    }

    /// Create backarc basin (behind volcanic arc, only for large plates)
    #[allow(clippy::too_many_arguments)]
    fn create_backarc_basin(
        &self,
        boundary: &BoundarySegment,
        overriding_plate: u16,
        boundary_idx: usize,
        arc_offset_km: f64,
        arc_width_km: f64,
        stats_a: &PlateStats,
        stats_b: &PlateStats,
        plate_map: &TerrainMap<u16>,
        _rng: &mut StdRng,
        regions: &mut Vec<ProvinceRegion>,
    ) {
        // Only create backarc basin for large plates
        if stats_a.area_km2 <= geo_const::BACKARC_BASIN_MIN_PLATE_AREA_KM2
            && stats_b.area_km2 <= geo_const::BACKARC_BASIN_MIN_PLATE_AREA_KM2 {
            return;
        }

        let offset_km = arc_offset_km + arc_width_km; // Start immediately after arc ends (no gap)

        // Dynamic width based on convergence rate (faster subduction = more backarc extension)
        let base_width_km = geo_const::BACKARC_BASIN_BASE_WIDTH_KM;
        let rate_above_min = (boundary.relative_velocity - geo_const::MIN_CONVERGENCE_RATE_CM_PER_YEAR).max(0.0);
        let multiplier = 1.0 + (geo_const::BACKARC_BASIN_WIDTH_FACTOR * rate_above_min);
        let clamped_multiplier = multiplier.min(geo_const::BACKARC_BASIN_MAX_MULTIPLIER);
        let width_km = base_width_km * clamped_multiplier;

        // Expand unidirectionally: from offset to offset+width
        // This creates a "ring" that doesn't overlap with previous features
        let outer_edge = self.expand_boundary_toward_plate_spherical(
            &boundary.pixels,
            overriding_plate,
            offset_km + width_km,
            plate_map
        );

        let inner_edge = self.expand_boundary_toward_plate_spherical(
            &boundary.pixels,
            overriding_plate,
            offset_km,
            plate_map
        );

        // Backarc basin = outer edge minus inner edge
        let inner_set: std::collections::HashSet<_> = inner_edge.iter().copied().collect();
        let backarc_pixels: Vec<(usize, usize)> = outer_edge.iter()
            .filter(|p| !inner_set.contains(p))
            .copied()
            .collect();

        let chars = ProvinceCharacteristics::backarc_basin(boundary.length_km);
        regions.push(ProvinceRegion::new(backarc_pixels, chars, Some(boundary_idx)));
    }

    /// Generate active continental rifts (volcanic rift zones like East African Rift)
    ///
    /// Active continental rifts feature extensive volcanism and flood basalt eruptions.
    /// These are ACTIVE divergent boundaries within continental plates.
    /// Only applies to continental-continental divergent boundaries.
    ///
    /// Examples: East African Rift, Rio Grande Rift, Baikal Rift
    fn generate_active_continental_rifts(
        &self,
        boundaries: &[BoundarySegment],
        plate_stats: &HashMap<u16, PlateStats>,
        plate_map: &TerrainMap<u16>,
        rng: &mut StdRng,
    ) -> Vec<ProvinceRegion> {
        let mut regions = Vec::new();

        for (idx, boundary) in boundaries.iter().enumerate() {
            if boundary.interaction_type != PlateInteraction::Divergent {
                continue;
            }

            let stats_a = plate_stats.get(&boundary.plate_a);
            let stats_b = plate_stats.get(&boundary.plate_b);

            if let (Some(a), Some(b)) = (stats_a, stats_b) {
                // ONLY continental-continental divergent boundaries
                // (Oceanic-oceanic divergence = mid-ocean ridge, handled separately)
                if a.plate_type == PlateType::Continental && b.plate_type == PlateType::Continental {
                    // Active continental rift (narrow linear zone of extension)
                    let characteristics = ProvinceCharacteristics::continental_rift(
                        boundary.relative_velocity,
                        boundary.length_km,
                    );

                    let width_km = geo_const::RIFT_MIN_WIDTH_KM + rng.gen::<f64>() * geo_const::RIFT_WIDTH_RANGE_KM;

                    // Use spherical-aware expansion
                    let pixels = self.expand_boundary_spherical(&boundary.pixels, width_km, plate_map);

                    regions.push(ProvinceRegion::new(pixels, characteristics, Some(idx)));
                }
                // Continental-oceanic and oceanic-oceanic divergence handled elsewhere
            }
        }

        regions
    }

    /// Generate oceanic base layer: Fill all oceanic plates with abyssal plains
    ///
    /// This creates the foundation layer for oceanic regions. All other oceanic features
    /// (ridges, trenches, arcs, etc.) will be layered on top of this base.
    fn generate_oceanic_base_layer(
        &self,
        plate_index: &PlatePixelIndex,
        plate_stats: &HashMap<u16, PlateStats>,
    ) -> Vec<ProvinceRegion> {
        let mut regions = Vec::new();

        // Sort plates by ID for deterministic iteration order
        let mut sorted_plates: Vec<_> = plate_stats.iter().collect();
        sorted_plates.sort_by_key(|(plate_id, _)| **plate_id);

        // Fill ALL oceanic plates with abyssal plains as background
        for (plate_id, stats) in sorted_plates {
            if stats.plate_type == PlateType::Oceanic {
                // Get pixels from pre-built index (O(1) lookup)
                if let Some(plate_pixels) = plate_index.get(plate_id) {
                    if !plate_pixels.is_empty() {
                        let chars = ProvinceCharacteristics::abyssal_plain(stats.area_km2 as f64);
                        regions.push(ProvinceRegion::new(plate_pixels.clone(), chars, None));
                    }
                }
            }
        }

        regions
    }

    /// Generate oceanic overlays: Add ridges and fractures on top of abyssal plains
    ///
    /// This creates features that overlay the abyssal plain base layer:
    /// 1. Mid-ocean ridges at divergent boundaries
    /// 2. Fracture zones at transform boundaries
    ///
    /// Note: Hotspot tracks are generated separately as final overlays in generate_all_provinces()
    fn generate_oceanic_overlays(
        &self,
        boundaries: &[BoundarySegment],
        plate_map: &TerrainMap<u16>,
        plate_stats: &HashMap<u16, PlateStats>,
    ) -> Vec<ProvinceRegion> {
        let mut regions = Vec::new();

        // First: Overlay mid-ocean ridges at divergent oceanic boundaries
        for (idx, boundary) in boundaries.iter().enumerate() {
            if boundary.interaction_type == PlateInteraction::Divergent {
                if let (Some(a), Some(b)) = (plate_stats.get(&boundary.plate_a), plate_stats.get(&boundary.plate_b)) {
                    if a.plate_type == PlateType::Oceanic && b.plate_type == PlateType::Oceanic {
                        let chars = ProvinceCharacteristics::mid_ocean_ridge(
                            boundary.relative_velocity,
                            boundary.length_km,
                        );

                        // Mid-ocean ridges: Width depends on spreading rate
                        let spreading_rate = boundary.relative_velocity; // cm/year
                        let width_km = if spreading_rate > geo_const::SPREADING_RATE_FAST {
                            geo_const::RIDGE_WIDTH_FAST_KM  // Fast-spreading: narrow, smooth
                        } else if spreading_rate > geo_const::SPREADING_RATE_MEDIUM {
                            geo_const::RIDGE_WIDTH_MEDIUM_KM  // Medium-spreading
                        } else if spreading_rate > geo_const::SPREADING_RATE_SLOW {
                            geo_const::RIDGE_WIDTH_SLOW_KM // Slow-spreading: wider with rift valley
                        } else {
                            geo_const::RIDGE_WIDTH_ULTRASLOW_KM // Ultra-slow: widest, most irregular
                        };

                        // Use spherical-aware expansion to account for latitude
                        let pixels = self.expand_boundary_spherical(&boundary.pixels, width_km, plate_map);

                        regions.push(ProvinceRegion::new(pixels, chars, Some(idx)));
                    }
                }
            }
        }

        // Second: Overlay fracture zones at transform boundaries (only for oceanic plates)
        for (idx, boundary) in boundaries.iter().enumerate() {
            if boundary.interaction_type == PlateInteraction::Transform {
                if let (Some(a), Some(b)) = (plate_stats.get(&boundary.plate_a), plate_stats.get(&boundary.plate_b)) {
                    if a.plate_type == PlateType::Oceanic || b.plate_type == PlateType::Oceanic {
                        let chars = ProvinceCharacteristics::oceanic_fracture_zone(boundary.length_km);

                        let width_km = geo_const::FRACTURE_ZONE_WIDTH_KM;

                        // Use spherical-aware expansion
                        let pixels = self.expand_boundary_spherical(&boundary.pixels, width_km, plate_map);

                        regions.push(ProvinceRegion::new(pixels, chars, Some(idx)));
                    }
                }
            }
        }

        regions
    }

    /// Generate stable continental base layer (platforms with scattered shields)
    ///
    /// Creates the ancient Precambrian basement (cratons) that forms the foundation of continents.
    /// This runs BEFORE orogens, so mountain building will compress/overlay these stable cores.
    ///
    /// Terminology: Craton = Shield + Platform (the entire stable basement)
    /// - Shield = Exposed Precambrian basement rock at surface
    /// - Platform = Sedimentary cover over the same cratonic basement
    ///
    /// Structure (Earth-like):
    /// 1. Platform (pink) fills ALL continental pixels (base layer)
    /// 2. Shields (orange) = 2-4 scattered fragments (not centered!) overlaying platform
    /// 3. Extended crust (yellow) at passive margins (plate edges)
    /// 4. Intracratonic basins (purple-grey) occasional subsided areas
    fn generate_stable_continental_base(
        &self,
        plate_stats: &HashMap<u16, PlateStats>,
        plate_index: &PlatePixelIndex,
        rng: &mut StdRng,
    ) -> Vec<ProvinceRegion> {
        let mut regions = Vec::new();

        // Sort plates by ID for deterministic iteration order
        let mut sorted_plates: Vec<_> = plate_stats.iter().collect();
        sorted_plates.sort_by_key(|(plate_id, _)| **plate_id);

        // Process each continental plate
        for (plate_id, stats) in sorted_plates {
            if stats.plate_type != PlateType::Continental {
                continue;
            }

            // Get pixels from pre-built index (O(1) lookup)
            let plate_pixels = match plate_index.get(plate_id) {
                Some(pixels) if !pixels.is_empty() => pixels,
                _ => continue,
            };

            // 1. PLATFORM (PINK) - Fill ALL continental pixels as base layer
            // This is the sedimentary-covered cratonic basement
            let platform_chars = ProvinceCharacteristics::platform(stats.area_km2 as f64);
            regions.push(ProvinceRegion::new(plate_pixels.clone(), platform_chars, None));

            // 2. CONTINENTAL CORE / SHIELD (ORANGE) - Single massive core
            self.generate_continental_shield(plate_pixels, &mut regions);

            // 3. INTRACRATONIC BASINS - Occasional subsided areas
            self.generate_intracratonic_basin(stats, plate_pixels, rng, &mut regions);
        }

        regions
    }

    /// Generate continental shield (exposed cratonic core)
    ///
    /// Real continents grow around a central cratonic nucleus (e.g., Canadian Shield).
    /// Generates one large cohesive core representing 20-30% of plate area.
    fn generate_continental_shield(
        &self,
        plate_pixels: &[(usize, usize)],
        regions: &mut Vec<ProvinceRegion>,
    ) {
        // Find rough centroid
        let min_x = plate_pixels.iter().map(|(x, _)| x).min().unwrap_or(&0);
        let max_x = plate_pixels.iter().map(|(x, _)| x).max().unwrap_or(&0);
        let min_y = plate_pixels.iter().map(|(_, y)| y).min().unwrap_or(&0);
        let max_y = plate_pixels.iter().map(|(_, y)| y).max().unwrap_or(&0);
        let plate_center_x = (min_x + max_x) / 2;
        let plate_center_y = (min_y + max_y) / 2;

        // Target size: fraction of plate area for the exposed shield
        let target_shield_area = (plate_pixels.len() as f64 * geo_const::SHIELD_AREA_FRACTION) as usize;
        let target_radius = (target_shield_area as f64 / std::f64::consts::PI).sqrt();

        // Generate the core using a noise-distorted distance field
        // We want a cohesive blob, not a perfect circle
        let shield_pixels: Vec<(usize, usize)> = plate_pixels.iter()
            .filter(|&&(x, y)| {
                let dx = (x as i32 - plate_center_x as i32) as f64;
                let dy = (y as i32 - plate_center_y as i32) as f64;
                let dist = (dx * dx + dy * dy).sqrt();

                // Simple noise variation to make edge irregular
                // (Using pseudo-randomness based on coords to be deterministic but varied)
                let angle = dy.atan2(dx);
                let noise = (angle * 3.0).sin() * 0.2 + (angle * 7.0).cos() * 0.1;
                let varying_radius = target_radius * (1.0 + noise);

                dist < varying_radius
            })
            .copied()
            .collect();

        if !shield_pixels.is_empty() {
            let chars = ProvinceCharacteristics::craton(shield_pixels.len() as f64 * 2500.0);
            regions.push(ProvinceRegion::new(shield_pixels, chars, None));
        }
    }

    /// Generate intracratonic basin (subsided area within stable craton)
    ///
    /// Occasional circular subsided areas found in large, stable cratons.
    /// Example: Michigan Basin, Illinois Basin
    fn generate_intracratonic_basin(
        &self,
        stats: &PlateStats,
        plate_pixels: &[(usize, usize)],
        rng: &mut StdRng,
        regions: &mut Vec<ProvinceRegion>,
    ) {
        // Only on large plates with low probability
        if stats.area_km2 <= geo_const::INTRACRATONIC_BASIN_MIN_AREA_KM2
            || rng.gen::<f64>() >= geo_const::INTRACRATONIC_BASIN_PROBABILITY {
            return;
        }

        // Pick random basin center
        let basin_x = plate_pixels[rng.gen_range(0..plate_pixels.len())].0;
        let basin_y = plate_pixels[rng.gen_range(0..plate_pixels.len())].1;

        // Generate circular basin
        let basin_pixels: Vec<(usize, usize)> = plate_pixels.iter()
            .filter(|&&(x, y)| {
                let dx = (x as i32 - basin_x as i32).abs() as f64;
                let dy = (y as i32 - basin_y as i32).abs() as f64;
                let dist = (dx * dx + dy * dy).sqrt();
                dist < geo_const::INTRACRATONIC_BASIN_MIN_RADIUS_PX
            })
            .copied()
            .collect();

        if basin_pixels.len() > geo_const::INTRACRATONIC_BASIN_MIN_PIXELS {
            let chars = ProvinceCharacteristics::intracratonic_basin(
                basin_pixels.len() as f64 * geo_const::INTRACRATONIC_BASIN_AREA_PER_PIXEL_KM2
            );
            regions.push(ProvinceRegion::new(basin_pixels, chars, None));
        }
    }

    /// Generate passive continental margins (extended crust)
    ///
    /// Creates extended/thinned crust at the edges of continental plates.
    /// These are passive margins - old rift scars where rifting stopped (e.g., Atlantic coast).
    /// WIDER at divergent boundaries (where rifting occurred), narrower elsewhere.
    ///
    /// Only applies to CONTINENTAL plates.
    fn generate_passive_margins(
        &self,
        boundaries: &[BoundarySegment],
        plate_stats: &HashMap<u16, PlateStats>,
        plate_map: &TerrainMap<u16>,
    ) -> Vec<ProvinceRegion> {
        let mut regions = Vec::new();

        for (idx, boundary) in boundaries.iter().enumerate() {
            let stats_a = plate_stats.get(&boundary.plate_a);
            let stats_b = plate_stats.get(&boundary.plate_b);

            if let (Some(a), Some(b)) = (stats_a, stats_b) {
                // Only create passive margins at continental-oceanic or continental-continental boundaries
                let is_continental_edge = a.plate_type == PlateType::Continental || b.plate_type == PlateType::Continental;

                if !is_continental_edge {
                    continue; // Skip oceanic-oceanic boundaries
                }

                // Determine width based on boundary type
                let width_km = if boundary.interaction_type == PlateInteraction::Divergent {
                    // WIDE margins at divergent boundaries (old rift zones like Atlantic)
                    geo_const::PASSIVE_MARGIN_DIVERGENT_WIDTH_KM
                } else {
                    // Narrower margins at convergent/transform boundaries
                    geo_const::PASSIVE_MARGIN_OTHER_WIDTH_KM
                };

                // Use spherical-aware expansion
                let pixels = self.expand_boundary_spherical(&boundary.pixels, width_km, plate_map);

                // Filter to only include continental pixels
                let continental_pixels: Vec<(usize, usize)> = pixels.iter()
                    .filter(|&&(x, y)| {
                        let idx = y * plate_map.width + x;
                        if idx >= plate_map.data.len() {
                            return false;
                        }
                        let plate_id = plate_map.data[idx];
                        if let Some(stats) = plate_stats.get(&plate_id) {
                            stats.plate_type == PlateType::Continental
                        } else {
                            false
                        }
                    })
                    .copied()
                    .collect();

                if !continental_pixels.is_empty() {
                    let chars = ProvinceCharacteristics::extended_crust(
                        boundary.length_km,
                    );
                    regions.push(ProvinceRegion::new(continental_pixels, chars, Some(idx)));
                }
            }
        }

        regions
    }

    /// Generate Large Igneous Provinces (rare volcanic provinces)
    fn generate_large_igneous_provinces(
        &self,
        plate_stats: &HashMap<u16, PlateStats>,
        plate_index: &PlatePixelIndex,
        existing_regions: &[ProvinceRegion],
        rng: &mut StdRng,
    ) -> Vec<ProvinceRegion> {
        let mut regions = Vec::new();

        // Collect assigned pixels
        let mut assigned_pixels: HashSet<(usize, usize)> = HashSet::new();
        for region in existing_regions {
            for &pixel in &region.pixels {
                assigned_pixels.insert(pixel);
            }
        }

        // CRITICAL: Sort plates by ID to ensure deterministic iteration order
        let mut sorted_plates: Vec<_> = plate_stats.iter().collect();
        sorted_plates.sort_by_key(|(plate_id, _)| **plate_id);

        // Rarely generate LIPs on plates
        for (plate_id, stats) in sorted_plates {
            if rng.gen::<f64>() < self.config.lip_probability {
                let (province_type, area) = if stats.plate_type == PlateType::Continental {
                    (GeologicProvince::ContinentalFloodBasalt, stats.area_km2 as f64 * geo_const::CONTINENTAL_FLOOD_BASALT_AREA_FRACTION)
                } else {
                    (GeologicProvince::OceanicPlateau, stats.area_km2 as f64 * geo_const::OCEANIC_PLATEAU_AREA_FRACTION)
                };

                let characteristics = match province_type {
                    GeologicProvince::ContinentalFloodBasalt =>
                        ProvinceCharacteristics::continental_flood_basalt(area),
                    GeologicProvince::OceanicPlateau =>
                        ProvinceCharacteristics::oceanic_plateau(area),
                    _ => continue,
                };

                let pixels_sample = self.sample_plate_interior(*plate_id, plate_index, geo_const::LIP_INTERIOR_SAMPLE_FRACTION, rng);
                let unassigned: Vec<(usize, usize)> = pixels_sample.into_iter()
                    .filter(|p| !assigned_pixels.contains(p))
                    .collect();

                if !unassigned.is_empty() {
                    regions.push(ProvinceRegion::new(unassigned, characteristics, None));
                }
            }
        }

        regions
    }

    /// Helper: Core flood-fill expansion algorithm
    ///
    /// Generic expansion that can be filtered by a predicate function.
    /// This is the common logic extracted from expand_boundary variants.
    fn expand_boundary_filtered<F>(
        &self,
        boundary_pixels: &[(usize, usize)],
        distance_pixels: usize,
        plate_map: &TerrainMap<u16>,
        mut filter: F,
    ) -> Vec<(usize, usize)>
    where
        F: FnMut(usize, usize, &TerrainMap<u16>) -> bool,
    {
        let mut result = HashSet::new();
        let mut current_layer: Vec<(usize, usize)> = boundary_pixels.to_vec();

        // Add initial boundary pixels
        for &pixel in boundary_pixels {
            result.insert(pixel);
        }

        // Iteratively expand outward
        for _ in 0..distance_pixels {
            let mut next_layer = Vec::new();

            for &(x, y) in &current_layer {
                let neighbors = plate_map.get_neighbors(x, y);
                for (nx, ny) in neighbors {
                    // Skip if already in result
                    if result.contains(&(nx, ny)) {
                        continue;
                    }

                    // Apply filter predicate
                    if filter(nx, ny, plate_map) {
                        result.insert((nx, ny));
                        next_layer.push((nx, ny));
                    }
                }
            }

            // Stop if no new pixels were added
            if next_layer.is_empty() {
                break;
            }

            current_layer = next_layer;
        }

        result.into_iter().collect()
    }

    /// Helper: expand boundary pixels in all directions
    fn expand_boundary(
        &self,
        boundary_pixels: &[(usize, usize)],
        width_pixels: usize,
        plate_map: &TerrainMap<u16>,
    ) -> Vec<(usize, usize)> {
        // Accept all pixels (no filtering)
        self.expand_boundary_filtered(boundary_pixels, width_pixels, plate_map, |_, _, _| true)
    }

    /// Helper: expand boundary pixels toward a specific plate
    ///
    /// This expands the boundary in the direction of the target plate only,
    /// allowing us to position features (arc, basins) on one side of the boundary.
    fn expand_boundary_toward_plate(
        &self,
        boundary_pixels: &[(usize, usize)],
        target_plate: u16,
        distance_pixels: usize,
        plate_map: &TerrainMap<u16>,
    ) -> Vec<(usize, usize)> {
        // Only accept pixels belonging to the target plate
        self.expand_boundary_filtered(boundary_pixels, distance_pixels, plate_map, move |x, y, map| {
            let idx = y * map.width + x;
            idx < map.data.len() && map.data[idx] == target_plate
        })
    }

    /// Helper: Spherical-aware expansion with latitude-dependent distance
    ///
    /// Unlike `expand_boundary()` which uses fixed pixel distance, this accounts for
    /// latitude: near poles, we need to expand MORE pixels to cover the same km distance
    /// because pixels are "narrower" east-west at high latitudes.
    ///
    /// Uses average latitude of boundary to determine expansion distance.
    ///
    /// # Arguments
    /// * `boundary_pixels` - Starting pixels
    /// * `target_width_km` - Desired width in kilometers
    /// * `plate_map` - The terrain map
    ///
    /// # Returns
    /// Expanded pixel set with latitude-corrected distances
    fn expand_boundary_spherical(
        &self,
        boundary_pixels: &[(usize, usize)],
        target_width_km: f64,
        plate_map: &TerrainMap<u16>,
    ) -> Vec<(usize, usize)> {
        if boundary_pixels.is_empty() {
            return Vec::new();
        }

        // Calculate average latitude of the boundary
        let avg_lat = {
            let mut sum = 0.0;
            for &(x, y) in boundary_pixels {
                let (lat, _) = plate_map.projection.pixel_to_coords(x, y);
                sum += lat;
            }
            sum / boundary_pixels.len() as f64
        };

        // Calculate effective km per pixel at this latitude
        let km_per_px_ns = self.km_per_pixel(plate_map);
        let lat_rad = avg_lat.to_radians();
        // At high latitudes, pixels are compressed east-west by cos(lat)
        // So we need MORE iterations to cover the same km
        let latitude_factor = lat_rad.cos().abs().max(0.1); // Avoid division by zero near poles
        let effective_km_per_px = km_per_px_ns * latitude_factor;

        // Calculate pixel distance needed at this latitude
        let distance_pixels = (target_width_km / effective_km_per_px).ceil() as usize;

        // Use standard expansion with corrected pixel distance
        self.expand_boundary(boundary_pixels, distance_pixels, plate_map)
    }

    /// Helper: Spherical-aware expansion toward a specific plate
    ///
    /// Like `expand_boundary_spherical()` but only expands onto the target plate.
    ///
    /// # Arguments
    /// * `boundary_pixels` - Starting pixels
    /// * `target_plate` - Plate ID to expand toward
    /// * `target_width_km` - Desired width in kilometers
    /// * `plate_map` - The terrain map
    ///
    /// # Returns
    /// Expanded pixel set with latitude-corrected distances, filtered to target plate
    fn expand_boundary_toward_plate_spherical(
        &self,
        boundary_pixels: &[(usize, usize)],
        target_plate: u16,
        target_width_km: f64,
        plate_map: &TerrainMap<u16>,
    ) -> Vec<(usize, usize)> {
        if boundary_pixels.is_empty() {
            return Vec::new();
        }

        // Calculate average latitude of the boundary
        let avg_lat = {
            let mut sum = 0.0;
            for &(x, y) in boundary_pixels {
                let (lat, _) = plate_map.projection.pixel_to_coords(x, y);
                sum += lat;
            }
            sum / boundary_pixels.len() as f64
        };

        // Calculate effective km per pixel at this latitude
        let km_per_px_ns = self.km_per_pixel(plate_map);
        let lat_rad = avg_lat.to_radians();
        let latitude_factor = lat_rad.cos().abs().max(0.1);
        let effective_km_per_px = km_per_px_ns * latitude_factor;

        // Calculate pixel distance needed at this latitude
        let distance_pixels = (target_width_km / effective_km_per_px).ceil() as usize;

        // Use standard expansion with corrected pixel distance
        self.expand_boundary_toward_plate(boundary_pixels, target_plate, distance_pixels, plate_map)
    }

    /// Helper: sample interior pixels from a plate
    ///
    /// Uses pre-built PlatePixelIndex for O(1) lookup instead of O(width × height) scanning.
    fn sample_plate_interior(
        &self,
        plate_id: u16,
        plate_index: &PlatePixelIndex,
        fraction: f64,
        rng: &mut StdRng,
    ) -> Vec<(usize, usize)> {
        // Get all pixels for this plate from the index (O(1) lookup)
        let plate_pixels = match plate_index.get(&plate_id) {
            Some(pixels) => pixels,
            None => return Vec::new(),
        };

        // Randomly sample the requested fraction
        plate_pixels.iter()
            .filter(|_| rng.gen::<f64>() < fraction)
            .copied()
            .collect()
    }

    /// Helper: find pixels in deep plate interior (far from boundaries)
    ///
    /// Filters for pixels that are at least `min_distance` pixels away from any plate boundary.
    /// This ensures features like hotspots spawn in plate centers, not near edges.
    ///
    /// Uses pre-built PlatePixelIndex for O(1) lookup instead of O(width × height) scanning.
    fn find_deep_interior_pixels(
        &self,
        plate_id: u16,
        plate_index: &PlatePixelIndex,
        plate_map: &TerrainMap<u16>,
        min_distance: usize,
    ) -> Vec<(usize, usize)> {
        // Get all pixels for this plate from the index (O(1) lookup)
        let plate_pixels = match plate_index.get(&plate_id) {
            Some(pixels) => pixels,
            None => return Vec::new(),
        };

        // Filter for deep interior pixels (far from any boundary)
        plate_pixels.iter()
            .filter(|&&(x, y)| {
                // Check if all neighbors within min_distance are same plate
                for dy in -(min_distance as i32)..=(min_distance as i32) {
                    for dx in -(min_distance as i32)..=(min_distance as i32) {
                        let nx = x as i32 + dx;
                        let ny = y as i32 + dy;

                        // Out of bounds = not interior
                        if nx < 0 || ny < 0 || nx >= plate_map.width as i32 || ny >= plate_map.height as i32 {
                            return false;
                        }

                        // Different plate = boundary nearby, not interior
                        let idx = (ny as usize) * plate_map.width + (nx as usize);
                        if idx < plate_map.data.len() && plate_map.data[idx] != plate_id {
                            return false;
                        }
                    }
                }
                true // All neighbors same plate = deep interior
            })
            .copied()
            .collect()
    }

    /// Helper: Calculate kilometers per pixel for the given map
    ///
    /// Uses the map's projection resolution and the planetary radius
    /// to determine the linear distance represented by one pixel.
    ///
    /// # Returns
    /// Kilometers per pixel (guaranteed to be positive and non-zero)
    fn km_per_pixel(&self, plate_map: &TerrainMap<u16>) -> f64 {
        let km_per_px = plate_map.projection.km_per_pixel(self.planetary_params.radius_km);
        // Guard against zero or negative values (should never happen with valid inputs)
        km_per_px.max(0.01)
    }

    /// Build an index mapping plate IDs to their pixel coordinates
    ///
    /// This scans the map once to build a lookup table, avoiding O(plates × pixels) complexity.
    /// Complexity: O(width × height)
    fn build_plate_pixel_index(&self, plate_map: &TerrainMap<u16>) -> PlatePixelIndex {
        let mut index: PlatePixelIndex = HashMap::new();

        for (y, row) in plate_map.data.chunks(plate_map.width).enumerate() {
            for (x, &plate_id) in row.iter().enumerate() {
                index.entry(plate_id)
                    .or_default()
                    .push((x, y));
            }
        }

        index
    }

    /// Generate rare hotspot tracks (linear volcanic chains)
    ///
    /// Creates age-progressive volcanic chains as plates move over stationary mantle hotspots.
    /// Earth has ~40-50 major hotspots globally, so these are very rare (2-5 per world).
    ///
    /// # Geological Process
    ///
    /// Hotspots are stationary plumes of hot mantle material that melt through moving plates,
    /// creating volcanic features. As the plate moves, new volcanoes form while older ones
    /// become extinct and erode, creating an age-progressive chain.
    ///
    /// # Chain Properties
    ///
    /// - **Direction**: Opposite to plate motion (oldest features in direction of motion)
    /// - **Length**: Based on plate velocity × time, representing only VISIBLE portions:
    ///   - Oceanic: 2,400 km max (Hawaiian Islands + atolls, 10-28 Ma)
    ///   - Continental: 400 km max (Yellowstone calderas, 2-5 Ma)
    /// - **Stopping**: Chains stop at plate boundaries where they get subducted/destroyed
    /// - **Location**: Deep in plate interiors (200-400+ km from boundaries)
    ///
    /// # Examples
    ///
    /// - **Oceanic**: Hawaiian-Emperor chain (6,200 km total, but only ~2,400 km visible)
    /// - **Continental**: Yellowstone hotspot track (Snake River Plain calderas)
    ///
    /// # Parameters
    ///
    /// - Only generated on large plates (>2M km²)
    /// - 10-15% probability per large plate
    /// - Results in 2-5 hotspots globally for typical 20-plate world
    fn generate_hotspot_tracks(
        &self,
        plate_stats: &HashMap<u16, PlateStats>,
        plate_map: &TerrainMap<u16>,
        plate_index: &PlatePixelIndex,
        rng: &mut StdRng,
    ) -> Vec<ProvinceRegion> {
        let mut regions = Vec::new();

        // CRITICAL: Sort plates by ID to ensure deterministic iteration order
        let mut sorted_plates: Vec<_> = plate_stats.iter().collect();
        sorted_plates.sort_by_key(|(plate_id, _)| **plate_id);

        for (plate_id, stats) in sorted_plates {
            // Only on large plates, very rare
            if stats.area_km2 < geo_const::HOTSPOT_MIN_PLATE_AREA_KM2 {
                continue;
            }

            // ~10-15% chance per large plate = 2-5 hotspots globally
            let probability = if stats.area_km2 > geo_const::HOTSPOT_LARGE_PLATE_THRESHOLD_KM2 {
                geo_const::HOTSPOT_PROBABILITY_LARGE_PLATE
            } else {
                geo_const::HOTSPOT_PROBABILITY_MEDIUM_PLATE
            };

            if rng.gen::<f64>() > probability {
                continue;
            }

            // Select hotspot location in deep plate interior
            let hotspot_location = match self.select_hotspot_location(
                *plate_id, stats.area_km2, plate_index, plate_map, rng
            ) {
                Some(loc) => loc,
                None => continue,
            };

            // Calculate chain direction (OPPOSITE to plate motion)
            let chain_azimuth = (stats.seed.motion_direction + 180.0) % 360.0;

            // Calculate chain parameters (length, width, characteristics)
            let (chain_length_km, width_km, chars) =
                self.calculate_hotspot_chain_params(stats, rng);

            if chain_length_km < geo_const::HOTSPOT_MIN_VISIBLE_LENGTH_KM {
                continue;
            }

            // Create and widen the chain
            let chain_pixels = self.create_linear_chain(
                hotspot_location,
                chain_azimuth,
                chain_length_km,
                *plate_id,
                plate_map,
            );

            if chain_pixels.is_empty() {
                continue;
            }

            let widened = self.expand_boundary_spherical(&chain_pixels, width_km, plate_map);
            regions.push(ProvinceRegion::new(widened, chars, None));
        }

        regions
    }

    /// Select hotspot location in deep plate interior (far from boundaries)
    ///
    /// Hotspots form in plate centers, not at edges.
    /// Returns None if no suitable interior location found.
    fn select_hotspot_location(
        &self,
        plate_id: u16,
        plate_area_km2: u64,
        plate_index: &PlatePixelIndex,
        plate_map: &TerrainMap<u16>,
        rng: &mut StdRng,
    ) -> Option<(usize, usize)> {
        let min_distance = if plate_area_km2 > geo_const::HOTSPOT_LARGE_PLATE_THRESHOLD_KM2 {
            geo_const::HOTSPOT_MIN_INTERIOR_DISTANCE_LARGE_PX
        } else {
            geo_const::HOTSPOT_MIN_INTERIOR_DISTANCE_MEDIUM_PX
        };

        let interior_pixels = self.find_deep_interior_pixels(
            plate_id, plate_index, plate_map, min_distance
        );

        if interior_pixels.is_empty() {
            return None;
        }

        let hotspot_idx = rng.gen_range(0..interior_pixels.len());
        Some(interior_pixels[hotspot_idx])
    }

    /// Calculate hotspot chain parameters (length, width, characteristics)
    ///
    /// Returns: (chain_length_km, width_km, characteristics)
    ///
    /// Calculates chain length based on VISIBLE, geologically active portions only:
    /// - Oceanic: Hawaiian Islands + atolls (0-28 Ma), ~2,400 km max
    /// - Continental: Yellowstone calderas (0-5 Ma), ~400 km max
    fn calculate_hotspot_chain_params(
        &self,
        stats: &PlateStats,
        rng: &mut StdRng,
    ) -> (f64, f64, ProvinceCharacteristics) {
        let (time_ma, max_length_km, width_km, chars_fn): (f64, f64, f64, fn(f64) -> ProvinceCharacteristics) =
            if stats.plate_type == PlateType::Oceanic {
                // Oceanic: young islands/atolls
                let time = geo_const::HOTSPOT_OCEANIC_MIN_AGE_MA
                    + rng.gen::<f64>() * geo_const::HOTSPOT_OCEANIC_AGE_RANGE_MA;
                (time,
                 geo_const::HOTSPOT_OCEANIC_MAX_LENGTH_KM,
                 geo_const::HOTSPOT_OCEANIC_WIDTH_KM,
                 ProvinceCharacteristics::oceanic_hotspot_track)
            } else {
                // Continental: recent calderas
                let time = geo_const::HOTSPOT_CONTINENTAL_MIN_AGE_MA
                    + rng.gen::<f64>() * geo_const::HOTSPOT_CONTINENTAL_AGE_RANGE_MA;
                (time,
                 geo_const::HOTSPOT_CONTINENTAL_MAX_LENGTH_KM,
                 geo_const::HOTSPOT_CONTINENTAL_WIDTH_KM,
                 ProvinceCharacteristics::continental_hotspot_track)
            };

        // Formula: length = velocity (cm/yr) × time (Ma) × conversion factor
        // Units: cm/yr × Ma × 1,000,000 yr/Ma ÷ 100,000 cm/km = km
        let chain_length_km = (stats.seed.motion_speed * time_ma
            * geo_const::HOTSPOT_LENGTH_CONVERSION_FACTOR).min(max_length_km);

        let chars = chars_fn(chain_length_km);

        (chain_length_km, width_km, chars)
    }

    /// Create a linear chain of pixels from a starting point along an azimuth
    ///
    /// Generates pixels in a line opposite to plate motion direction.
    /// Stops at plate boundaries (where chains get subducted/destroyed).
    /// Chain lengths represent geologically active/visible portions only.
    fn create_linear_chain(
        &self,
        start: (usize, usize),
        azimuth_degrees: f64,
        length_km: f64,
        plate_id: u16,
        plate_map: &TerrainMap<u16>,
    ) -> Vec<(usize, usize)> {
        let mut chain = Vec::new();
        chain.push(start);

        let km_per_px = self.km_per_pixel(plate_map);
        let num_steps = (length_km / km_per_px) as usize;

        // Convert azimuth to radians and calculate direction vector
        let azimuth_rad = azimuth_degrees.to_radians();
        let dx = azimuth_rad.sin();
        let dy = -azimuth_rad.cos(); // Negative because y increases downward

        let mut x = start.0 as f64;
        let mut y = start.1 as f64;

        for _ in 0..num_steps {
            x += dx;
            y += dy;

            let xi = x.round() as i32;
            let yi = y.round() as i32;

            // Check map bounds
            if xi < 0 || yi < 0 || xi >= plate_map.width as i32 || yi >= plate_map.height as i32 {
                break;
            }

            let (xu, yu) = (xi as usize, yi as usize);

            // Stop at plate boundaries (subduction destroys the chain)
            let idx = yu * plate_map.width + xu;
            if idx >= plate_map.data.len() || plate_map.data[idx] != plate_id {
                break;
            }

            chain.push((xu, yu));
        }

        chain
    }

    /// Generate paleo-orogens (ancient inactive mountain belts)
    ///
    /// "Ghosts" of past collisions, these linear features cross active cratons.
    /// They represent sutures from previous supercontinent cycles.
    ///
    /// Logic:
    /// - Generate 1-2 linear belts crossing continental interiors
    /// - Random orientation, but often traversing the entire continent
    /// - Width ~300 km (eroded roots)
    fn generate_paleo_orogens(
        &self,
        plate_stats: &HashMap<u16, PlateStats>,
        plate_index: &PlatePixelIndex,
        plate_map: &TerrainMap<u16>,
        rng: &mut StdRng,
    ) -> Vec<ProvinceRegion> {
        let mut regions = Vec::new();
        
        // Deterministic iteration
        let mut sorted_plates: Vec<_> = plate_stats.iter().collect();
        sorted_plates.sort_by_key(|(plate_id, _)| **plate_id);

        for (plate_id, stats) in sorted_plates {
            // Only continental plates
            if stats.plate_type != PlateType::Continental { continue; }

            // Only distinct on medium/large continents
            if stats.area_km2 < geo_const::PALEO_OROGEN_MIN_AREA_KM2 { continue; }

            // 1-2 belts per large continent
            let num_belts = if stats.area_km2 > geo_const::HOTSPOT_LARGE_PLATE_THRESHOLD_KM2 {
                geo_const::PALEO_OROGEN_COUNT_LARGE
            } else {
                geo_const::PALEO_OROGEN_COUNT_MEDIUM
            };
            
            if let Some(pixels) = plate_index.get(plate_id) {
                if pixels.len() < 100 { continue; }
                
                for _ in 0..num_belts {
                    // Pick start/end points
                    // Try a few times to get points far apart
                    let mut start = pixels[0];
                    let mut end = pixels[0];
                    let mut valid = false;
                    
                    for _ in 0..10 {
                        let p1 = pixels[rng.gen_range(0..pixels.len())];
                        let p2 = pixels[rng.gen_range(0..pixels.len())];
                        
                        let dx = (p1.0 as i32 - p2.0 as i32).pow(2);
                        let dy = (p1.1 as i32 - p2.1 as i32).pow(2);
                        let dist_sq = dx + dy;
                        
                        // Minimum length check (approximate spatial distance square)
                        if dist_sq > geo_const::PALEO_OROGEN_MIN_DISTANCE_SQ { // heuristic minimum distance
                            start = p1;
                            end = p2;
                            valid = true;
                            break;
                        }
                    }
                    
                    if !valid { continue; }
                    
                    // Rasterize a line between start and end using Bresenham's algorithm
                    let mut line_pixels = Vec::new();
                    let x0 = start.0 as i32;
                    let y0 = start.1 as i32;
                    let x1 = end.0 as i32;
                    let y1 = end.1 as i32;
                    
                    let dx = (x1 - x0).abs();
                    let dy = -(y1 - y0).abs();
                    let sx = if x0 < x1 { 1 } else { -1 };
                    let sy = if y0 < y1 { 1 } else { -1 };
                    let mut err = dx + dy;
                    
                    let mut cx = x0;
                    let mut cy = y0;
                    
                    loop {
                        line_pixels.push((cx as usize, cy as usize));
                        if cx == x1 && cy == y1 { break; }
                        let e2 = 2 * err;
                        if e2 >= dy { err += dy; cx += sx; }
                        if e2 <= dx { err += dx; cy += sy; }
                    }
                    
                    // Expand line to width
                    let width_km = geo_const::PALEO_OROGEN_MIN_WIDTH_KM + rng.gen::<f64>() * geo_const::PALEO_OROGEN_WIDTH_RANGE_KM;
                    let expanded = self.expand_boundary_spherical(&line_pixels, width_km / 2.0, plate_map);

                    // Filter to keep only on this continent
                    // (Paleo-orogens don't cross into the ocean usually, they are internal sutures)
                    let final_pixels: Vec<(usize, usize)> = expanded.iter()
                        .filter(|&&(x, y)| {
                            let idx = y * plate_map.width + x;
                            idx < plate_map.data.len() && plate_map.data[idx] == *plate_id
                        })
                        .copied()
                        .collect();

                    if final_pixels.len() > geo_const::PALEO_OROGEN_MIN_PIXELS {
                         let chars = ProvinceCharacteristics::paleo_orogen(width_km);
                         regions.push(ProvinceRegion::new(final_pixels, chars, None));
                    }
                }
            }
        }
        
        regions
    }
}
