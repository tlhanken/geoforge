//! Geological province generation (Stage 2 complete)
//!
//! This module coordinates all geological province generation from tectonic data,
//! implementing the complete Stage 2 pipeline.

use crate::geology::provinces::{GeologicProvince, ProvinceCharacteristics, ProvinceRegion};
use crate::geology::orogenic::{OrogenicBeltGenerator, OrogenicConfig};
use crate::map::terrain::TerrainMap;
use crate::map::spherical::PlanetaryParams;
use crate::tectonics::boundary_analysis::BoundarySegment;
use crate::tectonics::plates::{PlateInteraction, PlateStats, PlateType, PlateSeed};
use std::collections::{HashMap, HashSet};
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;

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
    pub fn generate_all_provinces(
        &self,
        boundaries: &[BoundarySegment],
        plate_stats: &HashMap<u16, PlateStats>,
        _plate_seeds: &[PlateSeed],  // Reserved for future use (hotspot tracks)
        plate_map: &TerrainMap<u16>,
    ) -> Vec<ProvinceRegion> {
        let mut regions = Vec::new();
        let mut rng = StdRng::seed_from_u64(self.seed);

        // ========== FOUNDATION LAYERS (Ancient/Stable) ==========

        // Stage 2.6a: Oceanic Base Layer - Fill all oceanic plates with abyssal plains
        // (Oldest oceanic crust, foundation for all oceanic features)
        if self.config.generate_oceanic {
            let oceanic_base = self.generate_oceanic_base_layer(plate_map, plate_stats);
            regions.extend(oceanic_base);
        }

        // Stage 2.4: Stable Continental Regions - Shields, platforms, intracratonic basins
        // (Ancient Precambrian basement = Cratons, 1.5-4 Ga old, foundation for all continental features)
        // Shield = exposed craton core, Platform = sedimentary-covered craton
        if self.config.generate_stable_regions {
            let stable = self.generate_stable_continental_base(plate_stats, plate_map, &mut rng);
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

        // Stage 2.6b: Oceanic Overlays - Ridges, fractures, seamounts
        // (Active spreading and transform features on oceanic base)
        if self.config.generate_oceanic {
            let oceanic_overlays = self.generate_oceanic_overlays(boundaries, plate_map, plate_stats, &mut rng);
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
                    // Determine which plate subducts (older/denser)
                    // Without plate age data, use plate ID as deterministic proxy
                    // Lower ID = formed earlier in generation = "older"
                    let (_subducting_plate, overriding_plate) = if boundary.plate_a < boundary.plate_b {
                        (boundary.plate_a, boundary.plate_b)
                    } else {
                        (boundary.plate_b, boundary.plate_a)
                    };

                    // 1. OCEAN TRENCH - at the boundary (deepest point)
                    // Elevation: -2.5 to -3.0 (Mariana Trench = -11km)
                    let trench_chars = ProvinceCharacteristics::ocean_trench(
                        boundary.relative_velocity,
                        boundary.length_km,
                    );
                    regions.push(ProvinceRegion::new(
                        boundary.pixels.clone(),
                        trench_chars,
                        Some(idx)
                    ));

                    // 2. FOREARC BASIN - between trench and arc (50-100km from trench)
                    // Elevation: -0.8 to -1.2 (moderate depth, sediment accumulation)
                    let forearc_width_km = 50.0 + rng.gen::<f64>() * 50.0;
                    let km_per_px = self.km_per_pixel(plate_map);
                    let forearc_pixels = self.expand_boundary_toward_plate(
                        &boundary.pixels,
                        overriding_plate,
                        (forearc_width_km / km_per_px) as usize,
                        plate_map
                    );

                    let forearc_chars = ProvinceCharacteristics::forearc_basin(
                        boundary.length_km,
                    );
                    regions.push(ProvinceRegion::new(forearc_pixels, forearc_chars, Some(idx)));

                    // 3. VOLCANIC ARC - on overriding plate (150-300km from trench)
                    // Elevation: +1.5 to +3.0 (may emerge as islands depending on convergence rate)
                    let arc_offset_km = 150.0 + rng.gen::<f64>() * 150.0;
                    let arc_width_km = 100.0 + rng.gen::<f64>() * 100.0;

                    let arc_pixels = self.expand_boundary_toward_plate(
                        &boundary.pixels,
                        overriding_plate,
                        (arc_offset_km / km_per_px) as usize,
                        plate_map
                    );

                    // Further expand to create arc width
                    let arc_pixels_wide = self.expand_boundary(&arc_pixels,
                        (arc_width_km / km_per_px / 2.0) as usize,
                        plate_map
                    );

                    let arc_chars = ProvinceCharacteristics::volcanic_arc(
                        boundary.relative_velocity,
                        boundary.length_km,
                    );
                    regions.push(ProvinceRegion::new(arc_pixels_wide, arc_chars, Some(idx)));

                    // 4. BACKARC BASIN - behind the volcanic arc (extensional)
                    // Elevation: -0.5 to -0.8 (shallow, spreading behind arc)
                    // Only create if there's enough space (large plates)
                    if stats_a.area_km2 > 500_000 || stats_b.area_km2 > 500_000 {
                        let backarc_offset_km = arc_offset_km + arc_width_km + 50.0;
                        let backarc_width_km = 100.0 + rng.gen::<f64>() * 100.0;

                        let backarc_pixels = self.expand_boundary_toward_plate(
                            &boundary.pixels,
                            overriding_plate,
                            (backarc_offset_km / km_per_px) as usize,
                            plate_map
                        );

                        let backarc_pixels_wide = self.expand_boundary(&backarc_pixels,
                            (backarc_width_km / km_per_px / 2.0) as usize,
                            plate_map
                        );

                        let backarc_chars = ProvinceCharacteristics::backarc_basin(
                            boundary.length_km,
                        );
                        regions.push(ProvinceRegion::new(backarc_pixels_wide, backarc_chars, Some(idx)));
                    }
                }
            }
        }

        regions
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
                    let km_per_px = self.km_per_pixel(plate_map);
                    let pixels = self.expand_boundary(&boundary.pixels,
                        (width_km / km_per_px) as usize,
                        plate_map
                    );

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
        plate_map: &TerrainMap<u16>,
        plate_stats: &HashMap<u16, PlateStats>,
    ) -> Vec<ProvinceRegion> {
        let mut regions = Vec::new();

        // Fill ALL oceanic plates with abyssal plains as background
        for (plate_id, stats) in plate_stats {
            if stats.plate_type == PlateType::Oceanic {
                // Collect all pixels belonging to this oceanic plate
                let mut plate_pixels = Vec::new();
                for (y, row) in plate_map.data.chunks(plate_map.width).enumerate() {
                    for (x, &pid) in row.iter().enumerate() {
                        if pid == *plate_id {
                            plate_pixels.push((x, y));
                        }
                    }
                }

                if !plate_pixels.is_empty() {
                    let chars = ProvinceCharacteristics::abyssal_plain(stats.area_km2 as f64);
                    regions.push(ProvinceRegion::new(plate_pixels, chars, None));
                }
            }
        }

        regions
    }

    /// Generate oceanic overlays: Add ridges, fractures, seamounts on top of abyssal plains
    ///
    /// This creates features that overlay the abyssal plain base layer:
    /// 1. Mid-ocean ridges at divergent boundaries
    /// 2. Fracture zones at transform boundaries
    /// 3. Seamount fields (random placement)
    fn generate_oceanic_overlays(
        &self,
        boundaries: &[BoundarySegment],
        plate_map: &TerrainMap<u16>,
        plate_stats: &HashMap<u16, PlateStats>,
        rng: &mut StdRng,
    ) -> Vec<ProvinceRegion> {
        let mut regions = Vec::new();

        // First: Overlay mid-ocean ridges at divergent oceanic boundaries
        let km_per_px = self.km_per_pixel(plate_map);
        for (idx, boundary) in boundaries.iter().enumerate() {
            if boundary.interaction_type == PlateInteraction::Divergent {
                if let (Some(a), Some(b)) = (plate_stats.get(&boundary.plate_a), plate_stats.get(&boundary.plate_b)) {
                    if a.plate_type == PlateType::Oceanic && b.plate_type == PlateType::Oceanic {
                        let chars = ProvinceCharacteristics::mid_ocean_ridge(
                            boundary.relative_velocity,
                            boundary.length_km,
                        );

                        let width_km = 400.0;
                        let pixels = self.expand_boundary(&boundary.pixels,
                            (width_km / km_per_px) as usize,
                            plate_map
                        );

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
                        let pixels = self.expand_boundary(&boundary.pixels,
                            (width_km / km_per_px) as usize,
                            plate_map
                        );

                        regions.push(ProvinceRegion::new(pixels, chars, Some(idx)));
                    }
                }
            }
        }

        // Third: Add seamount fields randomly to some oceanic plates
        for (plate_id, stats) in plate_stats {
            if stats.plate_type == PlateType::Oceanic && stats.area_km2 > 300_000 {
                // 20% chance for seamount field on large oceanic plates
                if rng.gen::<f64>() < 0.2 {
                    let area = stats.area_km2 as f64 * 0.05; // Small portion of plate
                    let chars = ProvinceCharacteristics::seamount_field(area);

                    // Random sample of pixels from this plate
                    let pixels = self.sample_plate_interior(*plate_id, plate_map, 0.05, rng);

                    if !pixels.is_empty() {
                        regions.push(ProvinceRegion::new(pixels, chars, None));
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
        plate_map: &TerrainMap<u16>,
        rng: &mut StdRng,
    ) -> Vec<ProvinceRegion> {
        let mut regions = Vec::new();

        // Process each continental plate
        for (plate_id, stats) in plate_stats {
            if stats.plate_type != PlateType::Continental {
                continue;
            }

            // Collect all pixels belonging to this continental plate
            let mut plate_pixels = Vec::new();
            for (y, row) in plate_map.data.chunks(plate_map.width).enumerate() {
                for (x, &pid) in row.iter().enumerate() {
                    if pid == *plate_id {
                        plate_pixels.push((x, y));
                    }
                }
            }

            if plate_pixels.is_empty() {
                continue;
            }

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

                let km_per_px = self.km_per_pixel(plate_map);
                let pixels = self.expand_boundary(&boundary.pixels,
                    (width_km / km_per_px) as usize,
                    plate_map
                );

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

        // Rarely generate LIPs on plates
        for (plate_id, stats) in plate_stats {
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

    /// Helper: expand boundary pixels in all directions
    fn expand_boundary(
        &self,
        boundary_pixels: &[(usize, usize)],
        width_pixels: usize,
        plate_map: &TerrainMap<u16>,
    ) -> Vec<(usize, usize)> {
        let mut result = HashSet::new();
        let mut to_expand: Vec<(usize, usize)> = boundary_pixels.to_vec();

        for &pixel in boundary_pixels {
            result.insert(pixel);
        }

        for _ in 0..width_pixels {
            let mut next_layer = Vec::new();

            for &(x, y) in &to_expand {
                let neighbors = plate_map.get_neighbors(x, y);
                for neighbor in neighbors {
                    if !result.contains(&neighbor) {
                        result.insert(neighbor);
                        next_layer.push(neighbor);
                    }
                }
            }

            if next_layer.is_empty() {
                break;
            }

            to_expand = next_layer;
        }

        result.into_iter().collect()
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
        let mut result = HashSet::new();
        let mut current_layer: Vec<(usize, usize)> = boundary_pixels.to_vec();

        for _ in 0..distance_pixels {
            let mut next_layer = Vec::new();

            for &(x, y) in &current_layer {
                let neighbors = plate_map.get_neighbors(x, y);
                for (nx, ny) in neighbors {
                    // Only expand into pixels belonging to the target plate
                    let idx = ny * plate_map.width + nx;
                    if idx < plate_map.data.len() && plate_map.data[idx] == target_plate {
                        if !result.contains(&(nx, ny)) {
                            result.insert((nx, ny));
                            next_layer.push((nx, ny));
                        }
                    }
                }
            }

            if next_layer.is_empty() {
                break;
            }

            current_layer = next_layer;
        }

        result.into_iter().collect()
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
}
