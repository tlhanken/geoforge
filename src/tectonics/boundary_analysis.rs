//! Boundary analysis and classification for tectonic plates
//!
//! This module identifies plate boundaries and classifies them based on plate motion.
//! It provides the foundation for Stage 1.4 (Plate Motion & Boundary Classification)
//! and Stage 2 (Geologic Provinces).

use crate::map::terrain::TerrainMap;
use crate::map::spherical::SphericalPoint;
use crate::tectonics::plates::{PlateSeed, PlateInteraction};
use std::collections::HashMap;

/// Configuration for boundary analysis
#[derive(Debug, Clone)]
pub struct BoundaryAnalysisConfig {
    /// Minimum boundary length in pixels to keep (filters noise)
    pub min_boundary_length: usize,

    /// Angle threshold in degrees for convergent/divergent classification
    /// Plates with relative motion angle < this are convergent/divergent
    /// Plates with angle near 90° are transform
    pub angle_threshold_degrees: f64,
}

impl Default for BoundaryAnalysisConfig {
    fn default() -> Self {
        Self {
            min_boundary_length: 10,
            angle_threshold_degrees: 45.0,
        }
    }
}

impl BoundaryAnalysisConfig {
    /// Create a new configuration with custom settings
    pub fn new(min_boundary_length: usize, angle_threshold_degrees: f64) -> Self {
        Self {
            min_boundary_length,
            angle_threshold_degrees,
        }
    }
}

/// Represents a segment of boundary between two plates
#[derive(Debug, Clone)]
pub struct BoundarySegment {
    /// First plate ID
    pub plate_a: u16,

    /// Second plate ID
    pub plate_b: u16,

    /// Pixels that form this boundary segment
    pub pixels: Vec<(usize, usize)>,

    /// Type of plate interaction at this boundary
    pub interaction_type: PlateInteraction,

    /// Relative velocity magnitude (cm/year)
    pub relative_velocity: f64,

    /// Approximate length of boundary in kilometers
    pub length_km: f64,
}

impl BoundarySegment {
    /// Create a new boundary segment without classification
    pub fn new(plate_a: u16, plate_b: u16, pixels: Vec<(usize, usize)>) -> Self {
        Self {
            plate_a,
            plate_b,
            pixels,
            interaction_type: PlateInteraction::Transform, // Default, will be classified
            relative_velocity: 0.0,
            length_km: 0.0,
        }
    }

    /// Get the number of pixels in this boundary
    pub fn pixel_count(&self) -> usize {
        self.pixels.len()
    }

    /// Check if this boundary involves the given plate
    pub fn involves_plate(&self, plate_id: u16) -> bool {
        self.plate_a == plate_id || self.plate_b == plate_id
    }

    /// Get the other plate ID given one plate
    pub fn other_plate(&self, plate_id: u16) -> Option<u16> {
        if self.plate_a == plate_id {
            Some(self.plate_b)
        } else if self.plate_b == plate_id {
            Some(self.plate_a)
        } else {
            None
        }
    }
}

/// Analyzes plate boundaries to identify and classify boundary segments
pub struct BoundaryAnalyzer {
    config: BoundaryAnalysisConfig,
}

impl BoundaryAnalyzer {
    /// Create a new boundary analyzer with default configuration
    pub fn new() -> Self {
        Self {
            config: BoundaryAnalysisConfig::default(),
        }
    }

    /// Create a new boundary analyzer with custom configuration
    pub fn with_config(config: BoundaryAnalysisConfig) -> Self {
        Self { config }
    }

    /// Find all plate boundaries in the given tectonic map
    ///
    /// This scans the entire map to identify pixels where neighbors belong to
    /// different plates, then groups them into boundary segments.
    ///
    /// # Arguments
    /// * `plate_map` - The tectonic plate map to analyze
    ///
    /// # Returns
    /// Vector of boundary segments, each representing a contiguous boundary between two plates
    pub fn find_boundaries(&self, plate_map: &TerrainMap<u16>) -> Vec<BoundarySegment> {
        // Map from plate pair to boundary pixels
        let mut boundary_pixels: HashMap<(u16, u16), Vec<(usize, usize)>> = HashMap::new();

        // Scan entire map for boundary pixels
        for y in 0..plate_map.height {
            for x in 0..plate_map.width {
                let current_plate = plate_map.get(x, y).copied().unwrap_or(0);
                if current_plate == 0 {
                    continue;
                }

                // Check neighbors for different plate IDs
                let neighbors = plate_map.get_neighbors(x, y);
                for (nx, ny) in neighbors {
                    let neighbor_plate = plate_map.get(nx, ny).copied().unwrap_or(0);

                    if neighbor_plate != 0 && neighbor_plate != current_plate {
                        // Found a boundary pixel
                        // Normalize plate pair (smaller ID first)
                        let plate_pair = if current_plate < neighbor_plate {
                            (current_plate, neighbor_plate)
                        } else {
                            (neighbor_plate, current_plate)
                        };

                        boundary_pixels.entry(plate_pair)
                            .or_default()
                            .push((x, y));
                    }
                }
            }
        }

        // Convert to boundary segments, filtering by minimum length
        let mut segments = Vec::new();
        for ((plate_a, plate_b), mut pixels) in boundary_pixels {
            // Remove duplicates (pixels can be added multiple times from neighbor checks)
            pixels.sort_unstable();
            pixels.dedup();

            if pixels.len() >= self.config.min_boundary_length {
                segments.push(BoundarySegment::new(plate_a, plate_b, pixels));
            }
        }

        segments
    }

    /// Classify boundary types based on plate motion vectors
    ///
    /// This calculates relative motion between plates at each boundary and classifies
    /// them as convergent, divergent, or transform based on the angle of relative motion.
    ///
    /// # Arguments
    /// * `segments` - Mutable reference to boundary segments to classify
    /// * `plate_seeds` - Plate seed data with motion vectors
    /// * `plate_map` - The tectonic plate map (for calculating boundary lengths)
    pub fn classify_boundaries(
        &self,
        segments: &mut [BoundarySegment],
        plate_seeds: &[PlateSeed],
        plate_map: &TerrainMap<u16>,
    ) {
        use crate::tectonics::motion::PlateMotionAssigner;

        // Create lookup map for seeds
        let seed_map: HashMap<u16, &PlateSeed> = plate_seeds
            .iter()
            .map(|seed| (seed.id, seed))
            .collect();

        for segment in segments.iter_mut() {
            // Get seeds for both plates
            let seed_a = seed_map.get(&segment.plate_a);
            let seed_b = seed_map.get(&segment.plate_b);

            if let (Some(seed_a), Some(seed_b)) = (seed_a, seed_b) {
                // Use improved classification from motion module
                segment.interaction_type = PlateMotionAssigner::classify_boundary_interaction(
                    seed_a,
                    seed_b,
                    self.config.angle_threshold_degrees,
                );

                // Calculate relative velocity magnitude
                let boundary_center = self.get_boundary_center(&segment.pixels, plate_map);
                let (velocity, _angle) = PlateMotionAssigner::calculate_relative_motion(
                    seed_a,
                    seed_b,
                    &boundary_center,
                );
                segment.relative_velocity = velocity;

                // Calculate boundary length
                segment.length_km = self.calculate_boundary_length(&segment.pixels, plate_map);
            }
        }
    }

    /// Get the approximate center point of a boundary
    fn get_boundary_center(&self, pixels: &[(usize, usize)], plate_map: &TerrainMap<u16>) -> SphericalPoint {
        if pixels.is_empty() {
            return SphericalPoint::from_lat_lon(0.0, 0.0);
        }

        let sum_x: usize = pixels.iter().map(|(x, _)| x).sum();
        let sum_y: usize = pixels.iter().map(|(_, y)| y).sum();
        let center_x = sum_x / pixels.len();
        let center_y = sum_y / pixels.len();

        plate_map.get_spherical_point(center_x, center_y)
    }

    /// Calculate approximate boundary length in kilometers
    fn calculate_boundary_length(
        &self,
        boundary_pixels: &[(usize, usize)],
        plate_map: &TerrainMap<u16>,
    ) -> f64 {
        if boundary_pixels.is_empty() {
            return 0.0;
        }

        // Approximate: multiply pixel count by average pixel size
        // Better would be to trace the actual path
        let mut total_area = 0.0;
        for (x, y) in boundary_pixels {
            let (lat, _) = plate_map.projection.pixel_to_coords(*x, *y);
            total_area += plate_map.projection.pixel_area_km2(lat);
        }

        // Rough estimate: sqrt of total area
        // (assumes boundary is roughly 1 pixel wide)
        total_area.sqrt() * boundary_pixels.len() as f64 / (boundary_pixels.len() as f64).sqrt()
    }

    /// Get statistics about all boundaries
    pub fn boundary_statistics(&self, segments: &[BoundarySegment]) -> BoundaryStatistics {
        let total_boundaries = segments.len();
        let convergent = segments.iter().filter(|s| s.interaction_type == PlateInteraction::Convergent).count();
        let divergent = segments.iter().filter(|s| s.interaction_type == PlateInteraction::Divergent).count();
        let transform = segments.iter().filter(|s| s.interaction_type == PlateInteraction::Transform).count();

        let total_length_km: f64 = segments.iter().map(|s| s.length_km).sum();
        let avg_velocity: f64 = if !segments.is_empty() {
            segments.iter().map(|s| s.relative_velocity).sum::<f64>() / segments.len() as f64
        } else {
            0.0
        };

        BoundaryStatistics {
            total_boundaries,
            convergent_count: convergent,
            divergent_count: divergent,
            transform_count: transform,
            total_length_km,
            average_relative_velocity: avg_velocity,
        }
    }
}

impl Default for BoundaryAnalyzer {
    fn default() -> Self {
        Self::new()
    }
}

/// Statistics about plate boundaries
#[derive(Debug, Clone)]
pub struct BoundaryStatistics {
    pub total_boundaries: usize,
    pub convergent_count: usize,
    pub divergent_count: usize,
    pub transform_count: usize,
    pub total_length_km: f64,
    pub average_relative_velocity: f64,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::map::terrain::TerrainMap;

    #[test]
    fn test_boundary_segment_creation() {
        let pixels = vec![(0, 0), (1, 0), (2, 0)];
        let segment = BoundarySegment::new(1, 2, pixels.clone());

        assert_eq!(segment.plate_a, 1);
        assert_eq!(segment.plate_b, 2);
        assert_eq!(segment.pixel_count(), 3);
        assert!(segment.involves_plate(1));
        assert!(segment.involves_plate(2));
        assert!(!segment.involves_plate(3));
        assert_eq!(segment.other_plate(1), Some(2));
        assert_eq!(segment.other_plate(2), Some(1));
        assert_eq!(segment.other_plate(3), None);
    }

    #[test]
    fn test_boundary_analyzer_creation() {
        let analyzer = BoundaryAnalyzer::new();
        assert_eq!(analyzer.config.min_boundary_length, 10);
        assert_eq!(analyzer.config.angle_threshold_degrees, 45.0);

        let config = BoundaryAnalysisConfig::new(5, 30.0);
        let analyzer = BoundaryAnalyzer::with_config(config);
        assert_eq!(analyzer.config.min_boundary_length, 5);
        assert_eq!(analyzer.config.angle_threshold_degrees, 30.0);
    }

    #[test]
    fn test_find_boundaries_simple() {
        // Create a simple 4x4 map with two plates
        let mut map = TerrainMap::new(4, 4, 0u16);

        // Left half = plate 1, right half = plate 2
        for y in 0..4 {
            for x in 0..2 {
                map.set(x, y, 1);
            }
            for x in 2..4 {
                map.set(x, y, 2);
            }
        }

        let analyzer = BoundaryAnalyzer::new();
        let boundaries = analyzer.find_boundaries(&map);

        // Should find one boundary between plates 1 and 2
        assert_eq!(boundaries.len(), 1);
        assert_eq!(boundaries[0].plate_a, 1);
        assert_eq!(boundaries[0].plate_b, 2);

        // Boundary should have pixels along x=1 and x=2
        assert!(boundaries[0].pixel_count() > 0);
    }

    #[test]
    fn test_find_boundaries_filters_short() {
        // Create a tiny boundary that should be filtered
        let mut map = TerrainMap::new(10, 10, 1u16);

        // Single pixel of plate 2
        map.set(5, 5, 2);

        let config = BoundaryAnalysisConfig::new(10, 45.0);
        let analyzer = BoundaryAnalyzer::with_config(config);
        let boundaries = analyzer.find_boundaries(&map);

        // Should be filtered out (too short)
        assert_eq!(boundaries.len(), 0);
    }
}
