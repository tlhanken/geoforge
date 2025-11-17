//! Comprehensive geological province generation (Stage 2 complete)
//!
//! This module coordinates all geological province generation from tecton data,
//! implementing the complete Stage 2 pipeline.

use crate::geology::provinces::{GeologicProvince, ProvinceCharacteristics, ProvinceRegion};
use crate::geology::orogenic::{OrogenicBeltGenerator, OrogenicConfig};
use crate::map::terrain::TerrainMap;
use crate::tectonics::boundary_analysis::BoundarySegment;
use crate::tectonics::plates::{PlateInteraction, PlateStats, PlateType, PlateSeed};
use std::collections::{HashMap, HashSet};
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;

/// Configuration for comprehensive geological province generation
#[derive(Debug, Clone)]
pub struct ComprehensiveGeologyConfig {
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

impl Default for ComprehensiveGeologyConfig {
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

/// Comprehensive geological province generator
pub struct ComprehensiveGeologyGenerator {
    config: ComprehensiveGeologyConfig,
    seed: u64,
}

impl ComprehensiveGeologyGenerator {
    pub fn new(config: ComprehensiveGeologyConfig, seed: u64) -> Self {
        Self { config, seed }
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

        // Stage 2.1: Orogenic Belts (from convergent boundaries)
        let orogenic_gen = OrogenicBeltGenerator::new(self.config.orogenic_config.clone());
        let orogens = orogenic_gen.generate_orogens(boundaries, plate_stats, plate_map);
        regions.extend(orogens);

        // Stage 2.3: Arc and Basin Systems (from oceanic-oceanic convergence)
        if self.config.generate_arc_systems {
            let arcs = self.generate_arc_systems(boundaries, plate_stats, plate_map, &mut rng);
            regions.extend(arcs);
        }

        // Stage 2.5: Extensional Zones (from divergent boundaries)
        if self.config.generate_extensional {
            let extensional = self.generate_extensional_zones(boundaries, plate_stats, plate_map, &mut rng);
            regions.extend(extensional);
        }

        // Stage 2.6: Oceanic Domains (mid-ocean ridges, trenches, abyssal plains)
        if self.config.generate_oceanic {
            let oceanic = self.generate_oceanic_domains(boundaries, plate_map, plate_stats, &mut rng);
            regions.extend(oceanic);
        }

        // Stage 2.4: Stable Continental Regions (cratons on large continental plates)
        if self.config.generate_stable_regions {
            let stable = self.generate_stable_regions(plate_stats, plate_map, &regions, &mut rng);
            regions.extend(stable);
        }

        // Stage 2.2: Large Igneous Provinces (rare, random placement)
        let lips = self.generate_large_igneous_provinces(plate_stats, plate_map, &regions, &mut rng);
        regions.extend(lips);

        regions
    }

    /// Generate arc and basin systems from oceanic-oceanic convergence
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
                    // Create volcanic arc
                    let arc_width_km = 200.0 + rng.gen::<f64>() * 100.0;
                    let characteristics = ProvinceCharacteristics::volcanic_arc(
                        boundary.relative_velocity,
                        boundary.length_km,
                    );

                    let pixels = self.expand_boundary(&boundary.pixels,
                        (arc_width_km / 71.0) as usize, // Rough km->pixels
                        plate_map
                    );

                    regions.push(ProvinceRegion::new(pixels, characteristics, Some(idx)));

                    // Create trench (deeper than arc)
                    let trench_chars = ProvinceCharacteristics::ocean_trench(
                        boundary.relative_velocity,
                        boundary.length_km,
                    );
                    regions.push(ProvinceRegion::new(
                        boundary.pixels.clone(),
                        trench_chars,
                        Some(idx)
                    ));
                }
            }
        }

        regions
    }

    /// Generate extensional zones from divergent boundaries
    fn generate_extensional_zones(
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
                let characteristics = if a.plate_type == PlateType::Continental || b.plate_type == PlateType::Continental {
                    // Continental rift
                    ProvinceCharacteristics::continental_rift(
                        boundary.relative_velocity,
                        boundary.length_km,
                    )
                } else {
                    // Mid-ocean ridge (handled in oceanic domains, skip here)
                    continue;
                };

                let width_km = 300.0 + rng.gen::<f64>() * 200.0;
                let pixels = self.expand_boundary(&boundary.pixels,
                    (width_km / 71.0) as usize,
                    plate_map
                );

                regions.push(ProvinceRegion::new(pixels, characteristics, Some(idx)));
            }
        }

        regions
    }

    /// Generate oceanic domains (ridges, abyssal plains, trenches, fracture zones, seamounts)
    ///
    /// Strategy:
    /// 1. Fill ALL oceanic plate pixels with abyssal plains (background)
    /// 2. Overlay mid-ocean ridges at divergent boundaries
    /// 3. Overlay fracture zones at transform boundaries
    /// 4. Overlay seamount fields randomly
    fn generate_oceanic_domains(
        &self,
        boundaries: &[BoundarySegment],
        plate_map: &TerrainMap<u16>,
        plate_stats: &HashMap<u16, PlateStats>,
        rng: &mut StdRng,
    ) -> Vec<ProvinceRegion> {
        let mut regions = Vec::new();

        // First: Fill ALL oceanic plates with abyssal plains as background
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

        // Second: Overlay mid-ocean ridges at divergent oceanic boundaries
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
                            (width_km / 71.0) as usize,
                            plate_map
                        );

                        regions.push(ProvinceRegion::new(pixels, chars, Some(idx)));
                    }
                }
            }
        }

        // Third: Overlay fracture zones at transform boundaries (only for oceanic plates)
        for (idx, boundary) in boundaries.iter().enumerate() {
            if boundary.interaction_type == PlateInteraction::Transform {
                if let (Some(a), Some(b)) = (plate_stats.get(&boundary.plate_a), plate_stats.get(&boundary.plate_b)) {
                    if a.plate_type == PlateType::Oceanic || b.plate_type == PlateType::Oceanic {
                        let chars = ProvinceCharacteristics::oceanic_fracture_zone(boundary.length_km);

                        let width_km = 100.0;
                        let pixels = self.expand_boundary(&boundary.pixels,
                            (width_km / 71.0) as usize,
                            plate_map
                        );

                        regions.push(ProvinceRegion::new(pixels, chars, Some(idx)));
                    }
                }
            }
        }

        // Fourth: Add seamount fields randomly to some oceanic plates
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

    /// Generate stable continental regions (cratons, platforms)
    fn generate_stable_regions(
        &self,
        plate_stats: &HashMap<u16, PlateStats>,
        plate_map: &TerrainMap<u16>,
        existing_regions: &[ProvinceRegion],
        rng: &mut StdRng,
    ) -> Vec<ProvinceRegion> {
        let mut regions = Vec::new();

        // Collect all pixels already assigned to geological provinces
        let mut assigned_pixels: HashSet<(usize, usize)> = HashSet::new();
        for region in existing_regions {
            for &pixel in &region.pixels {
                assigned_pixels.insert(pixel);
            }
        }

        // Generate cratons on large continental plates
        for (plate_id, stats) in plate_stats {
            if stats.plate_type == PlateType::Continental && stats.area_km2 > 1_000_000 {
                // Large continental plate - likely has a craton
                let chars = ProvinceCharacteristics::craton(stats.area_km2 as f64 * 0.4);

                // Find unassigned interior pixels
                let interior = self.sample_plate_interior(*plate_id, plate_map, 0.25, rng);
                let unassigned: Vec<(usize, usize)> = interior.into_iter()
                    .filter(|p| !assigned_pixels.contains(p))
                    .collect();

                if !unassigned.is_empty() {
                    regions.push(ProvinceRegion::new(unassigned, chars, None));
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

    /// Helper: expand boundary pixels
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
}
