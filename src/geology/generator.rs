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
//! **Collision Orogens (1)**: CollisionOrogen
//! **Subduction Systems (5)**: OceanTrench, AccretionaryWedge, ForearcBasin, VolcanicArc, BackarcBasin
//! **LIPs (3)**: ContinentalFloodBasalt, OceanicPlateau, ContinentalHotspotTrack*
//! **Stable (3)**: Craton, Platform, ExtendedCrust
//! **Oceanic (4)**: AbyssalPlain, MidOceanRidge, OceanicFractureZone, OceanicHotspotTrack*
//! **Deferred**: ContinentalRift (will be implemented later)
//!
//! *Hotspot tracks are generated separately as final overlays

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
            lip_probability: 0.1, // 10% chance per suitable location
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
        let lips = self.generate_large_igneous_provinces(plate_stats, plate_map, &regions, &mut rng);
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

            // Check if both plates are oceanic
            if let (Some(stats_a), Some(stats_b)) = (plate_stats.get(&boundary.plate_a), plate_stats.get(&boundary.plate_b)) {
                if stats_a.plate_type == PlateType::Oceanic && stats_b.plate_type == PlateType::Oceanic {
                    // Determine which plate subducts (older/denser - lower ID)
                    let (subducting_plate, overriding_plate) = if boundary.plate_a < boundary.plate_b {
                        (boundary.plate_a, boundary.plate_b)
                    } else {
                        (boundary.plate_b, boundary.plate_a)
                    };

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
        // Expand minimally (50 km) to make visible on map
        let trench_width_km = 50.0;

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
        // Base: 100 km, scales up to 200 km at high convergence rates
        let base_width_km = 100.0;
        let rate_above_min = (boundary.relative_velocity - 2.0).max(0.0); // 2 cm/yr minimum
        let multiplier = 1.0 + (0.15 * rate_above_min); // 15% wider per cm/year above minimum
        let clamped_multiplier = multiplier.min(2.0); // Max 2x base width (200 km)
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
        // Base: 100 km, scales up to 200 km at high convergence rates
        let base_width_km = 100.0;
        let rate_above_min = (boundary.relative_velocity - 2.0).max(0.0); // 2 cm/yr minimum
        let multiplier = 1.0 + (0.15 * rate_above_min); // 15% wider per cm/year above minimum
        let clamped_multiplier = multiplier.min(2.0); // Max 2x base width (200 km)
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
        // Base: 50 km, scales up to 100 km at high convergence rates
        let base_width_km = 50.0;
        let rate_above_min = (boundary.relative_velocity - 2.0).max(0.0); // 2 cm/yr minimum
        let multiplier = 1.0 + (0.15 * rate_above_min); // 15% wider per cm/year above minimum
        let clamped_multiplier = multiplier.min(2.0); // Max 2x base width (100 km)
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
        if stats_a.area_km2 <= 500_000 && stats_b.area_km2 <= 500_000 {
            return;
        }

        let offset_km = arc_offset_km + arc_width_km; // Start immediately after arc ends (no gap)

        // Dynamic width based on convergence rate (faster subduction = more backarc extension)
        // Base: 200 km, scales up to 400 km at high convergence rates
        let base_width_km = 200.0;
        let rate_above_min = (boundary.relative_velocity - 2.0).max(0.0); // 2 cm/yr minimum
        let multiplier = 1.0 + (0.15 * rate_above_min); // 15% wider per cm/year above minimum
        let clamped_multiplier = multiplier.min(2.0); // Max 2x base width (400 km)
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
                    // Active continental rift with flood basalts (igneous/volcanic)
                    let characteristics = ProvinceCharacteristics::continental_flood_basalt(
                        boundary.length_km,
                    );

                    let width_km = 250.0 + rng.gen::<f64>() * 150.0;

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
                        // Fast-spreading (>10 cm/yr): ~50-80 km wide, smooth gentle rise (East Pacific Rise)
                        // Slow-spreading (2-5 cm/yr): ~100-150 km wide, deep rift valley (Mid-Atlantic Ridge)
                        // Ultra-slow (<1 cm/yr): Can be wider and more irregular (Gakkel Ridge)
                        let spreading_rate = boundary.relative_velocity; // cm/year
                        let width_km = if spreading_rate > 10.0 {
                            60.0  // Fast-spreading: narrow, smooth
                        } else if spreading_rate > 5.0 {
                            80.0  // Medium-spreading
                        } else if spreading_rate > 2.0 {
                            120.0 // Slow-spreading: wider with rift valley
                        } else {
                            150.0 // Ultra-slow: widest, most irregular
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

                        let width_km = 100.0;

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

            // 2. SHIELDS (ORANGE) - Generate 2-4 scattered fragments (exposed basement)
            // Earth has ~35 cratons globally; large plates should have multiple fragments
            let num_shields = if stats.area_km2 > 5_000_000 {
                4 // Large plates get 4 shield fragments
            } else if stats.area_km2 > 2_000_000 {
                3 // Medium plates get 3
            } else if stats.area_km2 > 1_000_000 {
                2 // Smaller plates get 2
            } else {
                1 // Tiny plates get 1
            };

            // Find plate bounds for proper shield placement
            let min_x = plate_pixels.iter().map(|(x, _)| x).min().unwrap_or(&0);
            let max_x = plate_pixels.iter().map(|(x, _)| x).max().unwrap_or(&0);
            let min_y = plate_pixels.iter().map(|(_, y)| y).min().unwrap_or(&0);
            let max_y = plate_pixels.iter().map(|(_, y)| y).max().unwrap_or(&0);
            let plate_center_x = (min_x + max_x) / 2;
            let plate_center_y = (min_y + max_y) / 2;

            // Generate scattered shield centers (not all at centroid!)
            for _i in 0..num_shields {
                // Scatter shields across the plate interior (stay within 70% of plate radius)
                let offset_range = ((max_x - min_x) as f64 * 0.35) as i32;

                // Guard against single-pixel or very small plates
                let offset_x = if offset_range > 0 {
                    rng.gen_range(-offset_range..=offset_range)
                } else {
                    0
                };
                let offset_y = if offset_range > 0 {
                    rng.gen_range(-offset_range..=offset_range)
                } else {
                    0
                };

                let shield_center = (
                    (plate_center_x as i32 + offset_x).max(0) as usize,
                    (plate_center_y as i32 + offset_y).max(0) as usize,
                );

                let shield_size = if stats.area_km2 > 5_000_000 {
                    20.0 + rng.gen::<f64>() * 15.0 // 20-35 pixels radius (larger!)
                } else {
                    15.0 + rng.gen::<f64>() * 10.0  // 15-25 pixels radius
                };

                let shield_pixels: Vec<(usize, usize)> = plate_pixels.iter()
                    .filter(|&&(x, y)| {
                        let dx = (x as i32 - shield_center.0 as i32).abs() as f64;
                        let dy = (y as i32 - shield_center.1 as i32).abs() as f64;
                        let dist = (dx * dx + dy * dy).sqrt();
                        dist < shield_size
                    })
                    .copied()
                    .collect();

                if shield_pixels.len() > 30 {
                    let chars = ProvinceCharacteristics::craton(shield_pixels.len() as f64 * 5000.0);
                    regions.push(ProvinceRegion::new(shield_pixels, chars, None));
                }
            }

            // 3. INTRACRATONIC BASINS - Occasional subsided areas (10-15% chance on large plates)
            if stats.area_km2 > 3_000_000 && rng.gen::<f64>() < 0.15 {
                let basin_x = plate_pixels[rng.gen_range(0..plate_pixels.len())].0;
                let basin_y = plate_pixels[rng.gen_range(0..plate_pixels.len())].1;

                let basin_pixels: Vec<(usize, usize)> = plate_pixels.iter()
                    .filter(|&&(x, y)| {
                        let dx = (x as i32 - basin_x as i32).abs() as f64;
                        let dy = (y as i32 - basin_y as i32).abs() as f64;
                        let dist = (dx * dx + dy * dy).sqrt();
                        dist < 8.0
                    })
                    .copied()
                    .collect();

                if basin_pixels.len() > 10 {
                    let chars = ProvinceCharacteristics::intracratonic_basin(
                        basin_pixels.len() as f64 * 10_000.0
                    );
                    regions.push(ProvinceRegion::new(basin_pixels, chars, None));
                }
            }
        }

        regions
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
                    200.0 + 150.0 // 350 km
                } else {
                    // Narrower margins at convergent/transform boundaries
                    80.0 + 50.0 // 130 km
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
        plate_map: &TerrainMap<u16>,
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
                    (GeologicProvince::ContinentalFloodBasalt, stats.area_km2 as f64 * 0.1)
                } else {
                    (GeologicProvince::OceanicPlateau, stats.area_km2 as f64 * 0.08)
                };

                let characteristics = match province_type {
                    GeologicProvince::ContinentalFloodBasalt =>
                        ProvinceCharacteristics::continental_flood_basalt(area),
                    GeologicProvince::OceanicPlateau =>
                        ProvinceCharacteristics::oceanic_plateau(area),
                    _ => continue,
                };

                let pixels_sample = self.sample_plate_interior(*plate_id, plate_map, 0.1, rng);
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
    fn sample_plate_interior(
        &self,
        plate_id: u16,
        plate_map: &TerrainMap<u16>,
        fraction: f64,
        rng: &mut StdRng,
    ) -> Vec<(usize, usize)> {
        let mut pixels = Vec::new();

        for (y, row) in plate_map.data.chunks(plate_map.width).enumerate() {
            for (x, &pid) in row.iter().enumerate() {
                if pid == plate_id && rng.gen::<f64>() < fraction {
                    pixels.push((x, y));
                }
            }
        }

        pixels
    }

    /// Helper: find pixels in deep plate interior (far from boundaries)
    ///
    /// Filters for pixels that are at least `min_distance` pixels away from any plate boundary.
    /// This ensures features like hotspots spawn in plate centers, not near edges.
    fn find_deep_interior_pixels(
        &self,
        plate_id: u16,
        plate_map: &TerrainMap<u16>,
        min_distance: usize,
    ) -> Vec<(usize, usize)> {
        let mut interior = Vec::new();

        for (y, row) in plate_map.data.chunks(plate_map.width).enumerate() {
            for (x, &pid) in row.iter().enumerate() {
                if pid != plate_id {
                    continue;
                }

                // Check if all neighbors within min_distance are same plate
                let mut is_interior = true;
                'check: for dy in -(min_distance as i32)..=(min_distance as i32) {
                    for dx in -(min_distance as i32)..=(min_distance as i32) {
                        let nx = x as i32 + dx;
                        let ny = y as i32 + dy;

                        if nx < 0 || ny < 0 || nx >= plate_map.width as i32 || ny >= plate_map.height as i32 {
                            is_interior = false;
                            break 'check;
                        }

                        let idx = (ny as usize) * plate_map.width + (nx as usize);
                        if idx < plate_map.data.len() && plate_map.data[idx] != plate_id {
                            is_interior = false;
                            break 'check;
                        }
                    }
                }

                if is_interior {
                    interior.push((x, y));
                }
            }
        }

        interior
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
                    .or_insert_with(Vec::new)
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
        _plate_index: &PlatePixelIndex,
        rng: &mut StdRng,
    ) -> Vec<ProvinceRegion> {
        let mut regions = Vec::new();

        // CRITICAL: Sort plates by ID to ensure deterministic iteration order
        // HashMap iteration order is non-deterministic, which would cause RNG
        // to be called in different sequences on different runs
        let mut sorted_plates: Vec<_> = plate_stats.iter().collect();
        sorted_plates.sort_by_key(|(plate_id, _)| **plate_id);

        for (plate_id, stats) in sorted_plates {
            // Only on large plates (>2M km²), very rare
            if stats.area_km2 < 2_000_000 {
                continue;
            }

            // ~10-15% chance per large plate = 2-5 hotspots globally for typical world
            let probability = if stats.area_km2 > 10_000_000 {
                0.15 // Pacific-sized plates more likely
            } else {
                0.10
            };

            if rng.gen::<f64>() > probability {
                continue;
            }

            // Pick random hotspot location in DEEP PLATE INTERIOR (far from boundaries)
            // Hotspots form in plate centers, not at edges
            // Require at least 10-20 pixels from any boundary (200-400 km)
            let min_distance = if stats.area_km2 > 10_000_000 {
                20 // Large plates: 400+ km from edge
            } else {
                10 // Medium plates: 200+ km from edge
            };

            let interior_pixels = self.find_deep_interior_pixels(*plate_id, plate_map, min_distance);

            if interior_pixels.is_empty() {
                continue; // No suitable deep interior location found
            }

            let hotspot_idx = rng.gen_range(0..interior_pixels.len());
            let hotspot_location = interior_pixels[hotspot_idx];

            // Calculate chain direction (OPPOSITE to plate motion)
            let plate_azimuth = stats.seed.motion_direction;
            let chain_azimuth = (plate_azimuth + 180.0) % 360.0;

            // Calculate chain length based on VISIBLE, geologically active portions only
            //
            // We don't include heavily eroded/subsided portions that are about to be subducted.
            // This focuses on terrain-significant features.
            //
            // Oceanic example: Hawaiian chain
            //   - Total length: 6,200 km (0-70 Ma)
            //   - Visible portion: ~2,400 km (Hawaiian Islands + Northwestern Hawaiian atolls, 0-28 Ma)
            //   - Old seamounts (28-70 Ma): heavily subsided, eroded flat, not terrain-significant
            //
            // Continental example: Yellowstone hotspot track
            //   - Total length: ~800 km (0-16 Ma)
            //   - Visible portion: ~400 km (recent calderas 0-5 Ma)
            //   - Older calderas (5-16 Ma): buried under basalt flows, not visible
            //
            // Formula: length = velocity (cm/yr) × time (Ma) × 10,000 (conversion factor)
            let (time_ma, max_length_km) = if stats.plate_type == PlateType::Oceanic {
                // Oceanic: 10-28 Ma (young islands/atolls), max 2,400 km
                (10.0 + rng.gen::<f64>() * 18.0, 2400.0)
            } else {
                // Continental: 2-5 Ma (recent calderas), max 400 km
                (2.0 + rng.gen::<f64>() * 3.0, 400.0)
            };

            let chain_length_km = (stats.seed.motion_speed * time_ma * 10_000.0).min(max_length_km);

            if chain_length_km < 200.0 {
                continue; // Too short to be visible
            }

            // Create linear chain from hotspot location
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

            // Widen slightly for visibility
            // Real seamount chains: 50-200 km wide
            let width_km = if stats.plate_type == PlateType::Oceanic {
                50.0 // Oceanic: narrow seamount chain
            } else {
                100.0 // Continental: wider volcanic field
            };

            // Use spherical-aware expansion
            let widened = self.expand_boundary_spherical(&chain_pixels, width_km, plate_map);

            // Determine province type based on plate character
            let chars = if stats.plate_type == PlateType::Oceanic {
                ProvinceCharacteristics::oceanic_hotspot_track(chain_length_km)
            } else {
                ProvinceCharacteristics::continental_hotspot_track(chain_length_km)
            };

            regions.push(ProvinceRegion::new(widened, chars, None));
        }

        regions
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
}
