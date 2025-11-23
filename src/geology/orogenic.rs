//! Orogenic belt generation from convergent plate boundaries
//!
//! This module implements Stage 2.1: Mountain-building zones created at convergent boundaries.

use crate::geology::provinces::{GeologicProvince, ProvinceCharacteristics, ProvinceRegion};
use crate::map::terrain::TerrainMap;
use crate::map::spherical::PlanetaryParams;
use crate::tectonics::boundary_analysis::BoundarySegment;
use crate::tectonics::plates::{PlateInteraction, PlateStats, PlateType};
use std::collections::{HashMap, HashSet};

/// Configuration for orogenic belt generation
///
/// Controls how convergent boundaries are converted into mountain-building zones.
#[derive(Debug, Clone)]
pub struct OrogenicConfig {
    /// Base width for collision orogens in km (continental-continental)
    ///
    /// This is the baseline width; actual width is scaled by convergence rate.
    /// Earth examples: Himalayas ~2000 km, Alps ~1000 km
    pub collision_base_width_km: f64,

    /// Base width for subduction orogens in km (oceanic-continental)
    ///
    /// Earth examples: Andes ~500 km, Cascades ~300 km
    pub subduction_base_width_km: f64,

    /// Minimum convergence rate to create an orogen (cm/year)
    ///
    /// Boundaries with slower convergence are considered inactive.
    /// Earth typical convergence: 2-10 cm/year
    pub min_convergence_rate: f64,

    /// Width scaling factor based on convergence rate
    ///
    /// Higher convergence → wider orogenic belt (more crustal shortening)
    /// Formula: width = base_width * (1 + rate_factor * (rate - min_rate))
    pub convergence_rate_width_factor: f64,

    /// Maximum width multiplier (prevents unrealistically wide belts)
    pub max_width_multiplier: f64,
}

impl Default for OrogenicConfig {
    fn default() -> Self {
        Self {
            // Base widths (km)
            collision_base_width_km: 1000.0,
            subduction_base_width_km: 400.0,

            // Activity thresholds
            min_convergence_rate: 2.0, // cm/year

            // Dynamic scaling
            convergence_rate_width_factor: 0.15, // 15% wider per cm/year above minimum
            max_width_multiplier: 2.5,           // Max 2.5x base width
        }
    }
}

impl OrogenicConfig {
    /// Create a new configuration with custom settings
    pub fn new(
        collision_base_width_km: f64,
        subduction_base_width_km: f64,
        min_convergence_rate: f64,
    ) -> Self {
        Self {
            collision_base_width_km,
            subduction_base_width_km,
            min_convergence_rate,
            ..Default::default()
        }
    }

    /// Calculate dynamic width based on convergence rate
    ///
    /// Higher convergence rates create wider orogenic belts due to
    /// increased crustal shortening and deformation.
    fn calculate_width(&self, base_width: f64, convergence_rate: f64) -> f64 {
        let rate_above_min = (convergence_rate - self.min_convergence_rate).max(0.0);
        let multiplier = 1.0 + (self.convergence_rate_width_factor * rate_above_min);
        let clamped_multiplier = multiplier.min(self.max_width_multiplier);
        base_width * clamped_multiplier
    }
}

/// Generator for orogenic belts from tectonic boundaries
pub struct OrogenicBeltGenerator {
    config: OrogenicConfig,
    /// Planetary parameters (for map scale calculations and other properties)
    planetary_params: PlanetaryParams,
}

impl OrogenicBeltGenerator {
    /// Create a new orogenic belt generator with configuration and planetary parameters
    pub fn new(config: OrogenicConfig, planetary_params: PlanetaryParams) -> Self {
        Self { config, planetary_params }
    }

    /// Generate orogenic belts from convergent plate boundaries
    ///
    /// Analyzes convergent boundaries and creates mountain-building zones based on:
    /// - Plate types (oceanic vs continental)
    /// - Convergence rate
    /// - Boundary geometry
    ///
    /// # Arguments
    /// * `boundaries` - Classified plate boundary segments
    /// * `plate_stats` - Statistics and metadata for each plate
    /// * `plate_map` - The tectonic plate map
    ///
    /// # Returns
    /// Vector of orogenic belt regions
    pub fn generate_orogens(
        &self,
        boundaries: &[BoundarySegment],
        plate_stats: &HashMap<u16, PlateStats>,
        plate_map: &TerrainMap<u16>,
    ) -> Vec<ProvinceRegion> {
        let mut regions = Vec::new();

        for (boundary_idx, boundary) in boundaries.iter().enumerate() {
            // Only process convergent boundaries
            if boundary.interaction_type != PlateInteraction::Convergent {
                continue;
            }

            // Check if convergence rate is sufficient
            if boundary.relative_velocity < self.config.min_convergence_rate {
                continue;
            }

            // Classify the type of orogen
            if let Some(orogen_type) = self.classify_orogen_type(boundary, plate_stats) {
                // Calculate dynamic width based on convergence rate
                let width_km = self.calculate_orogen_width(orogen_type, boundary.relative_velocity);

                // Create the orogenic belt region
                let region = self.create_belt_region(
                    boundary,
                    orogen_type,
                    width_km,
                    boundary.relative_velocity,
                    plate_map,
                    plate_stats,
                    boundary_idx,
                );

                regions.push(region);
            }
        }

        regions
    }

    /// Classify what type of orogen forms at a convergent boundary
    ///
    /// Based on the plate types of the two converging plates:
    /// - Continental + Continental → Collision Orogen (90% chance) or Accretionary (10%)
    /// - Oceanic + Continental → Subduction Orogen (70%) or Accretionary (30%)
    /// - Oceanic + Oceanic → None (island arcs, handled in Stage 2.3)
    ///
    /// # Arguments
    /// * `boundary` - The convergent boundary to classify
    /// * `plate_stats` - Plate metadata with type information
    ///
    /// # Returns
    /// The type of orogen, or None if plates not found
    fn classify_orogen_type(
        &self,
        boundary: &BoundarySegment,
        plate_stats: &HashMap<u16, PlateStats>,
    ) -> Option<GeologicProvince> {
        let plate_a_stats = plate_stats.get(&boundary.plate_a)?;
        let plate_b_stats = plate_stats.get(&boundary.plate_b)?;

        let type_a = plate_a_stats.plate_type;
        let type_b = plate_b_stats.plate_type;

        // Classify based on plate type combination
        match (type_a, type_b) {
            // Continental-Continental collision
            // Pure collision orogens (Himalayas, Alps)
            (PlateType::Continental, PlateType::Continental) => {
                Some(GeologicProvince::CollisionOrogen)
            }

            // Oceanic-Continental convergence
            // NOTE: Now handled in Arc Systems (Stage 2.3) - generates VolcanicArc + AccretionaryWedge
            // This orogenic belt generator only handles collision orogens now
            (PlateType::Oceanic, PlateType::Continental)
            | (PlateType::Continental, PlateType::Oceanic) => None,

            // Oceanic-Oceanic → Skip (island arcs)
            // These will be handled in Stage 2.3 (Arc and Basin Systems)
            (PlateType::Oceanic, PlateType::Oceanic) => None,
        }
    }

    /// Calculate the width of an orogenic belt
    ///
    /// Width is dynamically calculated from base width and convergence rate.
    ///
    /// # Arguments
    /// * `orogen_type` - Type of orogen being created
    /// * `convergence_rate` - Convergence rate in cm/year
    ///
    /// # Returns
    /// Width in kilometers
    fn calculate_orogen_width(&self, orogen_type: GeologicProvince, convergence_rate: f64) -> f64 {
        let base_width = match orogen_type {
            GeologicProvince::CollisionOrogen => self.config.collision_base_width_km,
            // NOTE: SubductionOrogen and AccretionaryOrogen removed
            // These are now handled by Arc Systems generator as VolcanicArc + AccretionaryWedge
            _ => self.config.collision_base_width_km, // Fallback (only collision orogens now)
        };

        self.config.calculate_width(base_width, convergence_rate)
    }

    /// Create an orogenic belt region from a boundary
    ///
    /// Expands the boundary pixels into a belt with the specified width.
    ///
    /// # Arguments
    /// * `boundary` - The boundary segment to expand
    /// * `orogen_type` - Type of orogen to create
    /// * `width_km` - Width of the belt in kilometers
    /// * `convergence_rate` - Convergence rate in cm/year
    /// * `plate_map` - The tectonic plate map (for distance calculations)
    /// * `plate_stats` - Plate statistics for filtering by type
    /// * `boundary_idx` - Index of the source boundary
    ///
    /// # Returns
    /// Province region representing the orogenic belt
    fn create_belt_region(
        &self,
        boundary: &BoundarySegment,
        orogen_type: GeologicProvince,
        width_km: f64,
        convergence_rate: f64,
        plate_map: &TerrainMap<u16>,
        plate_stats: &HashMap<u16, PlateStats>,
        boundary_idx: usize,
    ) -> ProvinceRegion {
        // Calculate km per pixel using the map's projection and planetary radius
        let km_per_pixel = plate_map.projection.km_per_pixel(self.planetary_params.radius_km);

        let width_pixels = (width_km / km_per_pixel).ceil() as usize;

        // Expand boundary pixels into a belt (CONTINENTAL PLATES ONLY)
        let belt_pixels = self.expand_boundary_pixels(&boundary.pixels, width_pixels, plate_map, plate_stats);

        // Create characteristics based on orogen type
        let characteristics = match orogen_type {
            GeologicProvince::CollisionOrogen => {
                ProvinceCharacteristics::collision_orogen(convergence_rate, width_km)
            }
            // NOTE: SubductionOrogen and AccretionaryOrogen removed
            // Arc systems (VolcanicArc, AccretionaryWedge) handled separately in Stage 2.3
            _ => ProvinceCharacteristics::collision_orogen(convergence_rate, width_km), // Fallback
        };

        ProvinceRegion::new(belt_pixels, characteristics, Some(boundary_idx))
    }

    /// Expand boundary pixels into a belt region (CONTINENTAL PLATES ONLY)
    ///
    /// Uses flood-fill expansion to create a belt perpendicular to the boundary.
    /// **IMPORTANT**: Only expands onto CONTINENTAL plates. Orogens cannot form on oceanic crust.
    ///
    /// # Arguments
    /// * `boundary_pixels` - Initial boundary pixels
    /// * `width_pixels` - How many pixels to expand on each side
    /// * `plate_map` - The map (for bounds checking and plate ID lookup)
    /// * `plate_stats` - Plate statistics for filtering by plate type
    ///
    /// # Returns
    /// Expanded set of pixels forming the belt (continental pixels only)
    fn expand_boundary_pixels(
        &self,
        boundary_pixels: &[(usize, usize)],
        width_pixels: usize,
        plate_map: &TerrainMap<u16>,
        plate_stats: &HashMap<u16, PlateStats>,
    ) -> Vec<(usize, usize)> {
        let mut result = HashSet::new();
        let mut to_expand: Vec<(usize, usize)> = boundary_pixels.to_vec();

        // Add initial boundary pixels
        for &pixel in boundary_pixels {
            result.insert(pixel);
        }

        // Expand iteratively, ONLY onto continental plates
        for _ in 0..width_pixels {
            let mut next_layer = Vec::new();

            for &(x, y) in &to_expand {
                // Get neighbors
                let neighbors = plate_map.get_neighbors(x, y);

                for (nx, ny) in neighbors {
                    if !result.contains(&(nx, ny)) {
                        // Check if this pixel is on a continental plate
                        let idx = ny * plate_map.width + nx;
                        if idx < plate_map.data.len() {
                            let plate_id = plate_map.data[idx];

                            // ONLY expand onto CONTINENTAL plates
                            if let Some(stats) = plate_stats.get(&plate_id) {
                                if stats.plate_type == PlateType::Continental {
                                    result.insert((nx, ny));
                                    next_layer.push((nx, ny));
                                }
                            }
                        }
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
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_config() {
        let config = OrogenicConfig::default();
        assert_eq!(config.collision_base_width_km, 1000.0);
        assert_eq!(config.subduction_base_width_km, 400.0);
        assert_eq!(config.min_convergence_rate, 2.0);
    }

    #[test]
    fn test_dynamic_width_calculation() {
        let config = OrogenicConfig::default();

        // At minimum convergence rate (2 cm/yr), width = base width
        let width_at_min = config.calculate_width(1000.0, 2.0);
        assert_eq!(width_at_min, 1000.0);

        // At 5 cm/yr, width should be larger
        let width_at_5 = config.calculate_width(1000.0, 5.0);
        assert!(width_at_5 > 1000.0);
        // Expected: 1000 * (1 + 0.15 * 3) = 1450
        assert!((width_at_5 - 1450.0).abs() < 0.1);

        // At very high rate, should clamp to max multiplier
        let width_at_high = config.calculate_width(1000.0, 50.0);
        assert!(width_at_high <= 1000.0 * config.max_width_multiplier);
    }

    #[test]
    fn test_orogen_classification_collision() {
        let generator = OrogenicBeltGenerator::new(
            OrogenicConfig::default(),
            PlanetaryParams::earth()
        );
        let mut plate_stats = HashMap::new();

        // Create two continental plates
        let stats_a = PlateStats {
            pixels: 10000,
            percentage: 10.0,
            area_km2: 1000000,
            seed: crate::tectonics::plates::PlateSeed::new(
                1, 100, 100, 0.0, 0.0, 0.0, 5.0
            ),
            plate_type: PlateType::Continental,
        };

        let stats_b = PlateStats {
            pixels: 9000,
            percentage: 9.0,
            area_km2: 900000,
            seed: crate::tectonics::plates::PlateSeed::new(
                2, 200, 200, 10.0, 10.0, 180.0, 5.0
            ),
            plate_type: PlateType::Continental,
        };

        plate_stats.insert(1, stats_a);
        plate_stats.insert(2, stats_b);

        let boundary = BoundarySegment {
            plate_a: 1,
            plate_b: 2,
            pixels: vec![(150, 150)],
            interaction_type: PlateInteraction::Convergent,
            relative_velocity: 5.0,
            length_km: 1000.0,
        };

        let orogen_type = generator.classify_orogen_type(&boundary, &plate_stats);
        // Continental-continental convergence produces collision orogens
        assert_eq!(
            orogen_type,
            Some(GeologicProvince::CollisionOrogen),
            "Continental-continental should produce collision orogen, got {:?}",
            orogen_type
        );
    }

    #[test]
    fn test_orogen_classification_oceanic_continental() {
        let generator = OrogenicBeltGenerator::new(
            OrogenicConfig::default(),
            PlanetaryParams::earth()
        );
        let mut plate_stats = HashMap::new();

        // Oceanic plate
        let stats_a = PlateStats {
            pixels: 5000,
            percentage: 5.0,
            area_km2: 500000,
            seed: crate::tectonics::plates::PlateSeed::new(
                1, 100, 100, 0.0, 0.0, 0.0, 8.0
            ),
            plate_type: PlateType::Oceanic,
        };

        // Continental plate
        let stats_b = PlateStats {
            pixels: 12000,
            percentage: 12.0,
            area_km2: 1200000,
            seed: crate::tectonics::plates::PlateSeed::new(
                2, 200, 200, 10.0, 10.0, 180.0, 3.0
            ),
            plate_type: PlateType::Continental,
        };

        plate_stats.insert(1, stats_a);
        plate_stats.insert(2, stats_b);

        let boundary = BoundarySegment {
            plate_a: 1,
            plate_b: 2,
            pixels: vec![(150, 150)],
            interaction_type: PlateInteraction::Convergent,
            relative_velocity: 8.0,
            length_km: 2000.0,
        };

        let orogen_type = generator.classify_orogen_type(&boundary, &plate_stats);
        // Oceanic-continental convergence now handled by Arc Systems (Stage 2.3)
        // Orogenic belt generator returns None for these
        assert_eq!(
            orogen_type,
            None,
            "Oceanic-continental should return None (handled by Arc Systems), got {:?}",
            orogen_type
        );
    }

    #[test]
    fn test_width_calculation_collision() {
        let generator = OrogenicBeltGenerator::new(
            OrogenicConfig::default(),
            PlanetaryParams::earth()
        );

        let collision_width = generator.calculate_orogen_width(
            GeologicProvince::CollisionOrogen,
            5.0,
        );

        // Collision orogens should be wide and dynamically scaled
        assert!(collision_width > 1000.0, "Collision width should be > 1000 km, got {}", collision_width);
    }
}
