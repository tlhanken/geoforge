//! Boundary refinement for tectonic plates
//!
//! This module provides post-processing to add realistic irregularity to plate boundaries,
//! transforming perfect Voronoi boundaries into more natural, rough edges that resemble
//! real tectonic plates.

use crate::map::terrain::TerrainMap;
use crate::tectonics::plates::PlateInteraction;
use std::collections::HashSet;

/// Configuration for boundary refinement
#[derive(Debug, Clone)]
pub struct BoundaryRefinementConfig {
    /// Random seed for deterministic generation
    pub seed: u64,

    /// Noise scale - larger values create larger irregularities (default: 0.1)
    /// Range: 0.01 (fine detail) to 1.0 (large features)
    pub noise_scale: f64,

    /// Noise amplitude - how far boundaries can shift (default: 3.0)
    /// Measured in pixels. Range: 1.0 to 10.0
    pub noise_amplitude: f64,

    /// Number of octaves for multi-scale noise (default: 3)
    /// More octaves = more detail at different scales
    pub octaves: usize,

    /// Persistence - how much each octave contributes (default: 0.5)
    /// Range: 0.0 (first octave dominates) to 1.0 (all equal)
    pub persistence: f64,

    /// Post-refinement smoothing iterations (default: 1)
    /// Prevents over-jaggedness while maintaining character
    pub smoothing_iterations: usize,

    /// Boundary type variation - scale different boundary types differently
    pub convergent_scale: f64,  // default: 1.2 (more irregular)
    pub divergent_scale: f64,   // default: 1.0 (normal)
    pub transform_scale: f64,   // default: 0.8 (less irregular)
}

impl Default for BoundaryRefinementConfig {
    fn default() -> Self {
        Self {
            seed: 0,
            noise_scale: 0.015,       // Medium-scale features
            noise_amplitude: 20.0,    // Moderate warping
            octaves: 4,               // Multi-scale detail
            persistence: 0.5,
            smoothing_iterations: 1,
            convergent_scale: 1.2,
            divergent_scale: 1.0,
            transform_scale: 0.8,
        }
    }
}

impl BoundaryRefinementConfig {
    /// Create a new configuration with a specific seed
    pub fn with_seed(seed: u64) -> Self {
        Self {
            seed,
            ..Default::default()
        }
    }

    /// Set noise parameters
    pub fn with_noise(mut self, scale: f64, amplitude: f64, octaves: usize) -> Self {
        self.noise_scale = scale;
        self.noise_amplitude = amplitude;
        self.octaves = octaves;
        self
    }

    /// Set smoothing iterations
    pub fn with_smoothing(mut self, iterations: usize) -> Self {
        self.smoothing_iterations = iterations;
        self
    }

    /// Set boundary type variations
    pub fn with_boundary_scales(mut self, convergent: f64, divergent: f64, transform: f64) -> Self {
        self.convergent_scale = convergent;
        self.divergent_scale = divergent;
        self.transform_scale = transform;
        self
    }
}

/// Boundary refinement processor
pub struct BoundaryRefiner {
    config: BoundaryRefinementConfig,
}

impl BoundaryRefiner {
    /// Create a new boundary refiner with configuration
    pub fn new(config: BoundaryRefinementConfig) -> Self {
        Self { config }
    }

    /// Refine the boundaries of a plate map using domain warping
    pub fn refine_boundaries(&mut self, plate_map: &mut TerrainMap<u16>) {
        println!("Refining plate boundaries with domain warping...");
        println!("  Noise scale: {}", self.config.noise_scale);
        println!("  Amplitude: {} pixels", self.config.noise_amplitude);
        println!("  Octaves: {}", self.config.octaves);

        // Use domain warping: warp the entire Voronoi diagram
        // instead of just perturbing boundary pixels
        self.apply_domain_warping(plate_map);

        // Optional smoothing to prevent over-jaggedness
        if self.config.smoothing_iterations > 0 {
            println!("  Applying {} smoothing iterations", self.config.smoothing_iterations);
            self.smooth_boundaries(plate_map, self.config.smoothing_iterations);
        }

        println!("Boundary refinement complete!");
    }

    /// Identify all boundary pixels (pixels adjacent to different plates)
    fn identify_boundaries(&self, plate_map: &TerrainMap<u16>) -> HashSet<(usize, usize)> {
        let mut boundaries = HashSet::new();

        for y in 0..plate_map.height {
            for x in 0..plate_map.width {
                let current_plate = plate_map.data[plate_map.get_index(x, y)];

                // Check if any neighbor has a different plate ID
                let neighbors = plate_map.get_neighbors(x, y);
                for (nx, ny) in neighbors {
                    let neighbor_plate = plate_map.data[plate_map.get_index(nx, ny)];
                    if neighbor_plate != current_plate {
                        boundaries.insert((x, y));
                        break;
                    }
                }
            }
        }

        boundaries
    }

    /// Apply domain warping to warp the entire Voronoi diagram
    ///
    /// Instead of perturbing boundary pixels, this warps the lookup coordinates
    /// for the entire map, creating large-scale organic distortions
    fn apply_domain_warping(&mut self, plate_map: &mut TerrainMap<u16>) {
        // Create a copy of the original data
        let original_data = plate_map.data.clone();

        println!("  Applying domain warping to {} pixels...", plate_map.width * plate_map.height);

        // Process each pixel with warped lookup
        for y in 0..plate_map.height {
            for x in 0..plate_map.width {
                // Generate noise offsets for this pixel
                // Use two different noise samples for X and Y to avoid correlation
                let offset_x = self.multi_octave_noise(
                    x as f64 * self.config.noise_scale,
                    y as f64 * self.config.noise_scale,
                ) * self.config.noise_amplitude;

                let offset_y = self.multi_octave_noise(
                    x as f64 * self.config.noise_scale + 1000.0,  // Offset to decorrelate
                    y as f64 * self.config.noise_scale + 1000.0,
                ) * self.config.noise_amplitude;

                // Compute warped sample position
                let sample_x_f = x as f64 + offset_x;
                let sample_y_f = y as f64 + offset_y;

                // Wrap X coordinate (longitude wraps around)
                let sample_x = sample_x_f.rem_euclid(plate_map.width as f64) as usize;
                // Clamp Y coordinate (latitude doesn't wrap)
                let sample_y = sample_y_f.clamp(0.0, (plate_map.height - 1) as f64) as usize;

                // Sample from warped position
                let sampled_plate = original_data[plate_map.get_index(sample_x, sample_y)];

                // Assign to current position
                let idx = plate_map.get_index(x, y);
                plate_map.data[idx] = sampled_plate;
            }

            if y % 100 == 0 && y > 0 {
                println!("    Progress: {:.1}%", (y as f64 / plate_map.height as f64) * 100.0);
            }
        }
    }

    /// Multi-octave Perlin-like noise
    fn multi_octave_noise(&mut self, x: f64, y: f64) -> f64 {
        let mut total = 0.0;
        let mut amplitude = 1.0;
        let mut frequency = 1.0;
        let mut max_value = 0.0;

        for _ in 0..self.config.octaves {
            total += self.simple_noise(x * frequency, y * frequency) * amplitude;
            max_value += amplitude;

            amplitude *= self.config.persistence;
            frequency *= 2.0;
        }

        total / max_value
    }

    /// Simple noise function using value noise
    fn simple_noise(&self, x: f64, y: f64) -> f64 {
        // Grid cell coordinates
        let x0 = x.floor() as i64;
        let y0 = y.floor() as i64;
        let x1 = x0 + 1;
        let y1 = y0 + 1;

        // Local coordinates within cell
        let fx = x - x0 as f64;
        let fy = y - y0 as f64;

        // Smooth interpolation
        let sx = self.smoothstep(fx);
        let sy = self.smoothstep(fy);

        // Get random values at grid corners
        let n00 = self.hash_noise(x0, y0);
        let n10 = self.hash_noise(x1, y0);
        let n01 = self.hash_noise(x0, y1);
        let n11 = self.hash_noise(x1, y1);

        // Bilinear interpolation
        let nx0 = self.lerp(n00, n10, sx);
        let nx1 = self.lerp(n01, n11, sx);
        self.lerp(nx0, nx1, sy)
    }

    /// Hash function for noise generation
    fn hash_noise(&self, x: i64, y: i64) -> f64 {
        // Simple hash using prime numbers and the seed
        let mut hash = self.config.seed.wrapping_mul(2654435761);
        hash = hash.wrapping_add(x as u64).wrapping_mul(2654435761);
        hash = hash.wrapping_add(y as u64).wrapping_mul(2654435761);

        // Convert to -1.0 to 1.0 range
        (hash as f64 / u64::MAX as f64) * 2.0 - 1.0
    }

    /// Smooth interpolation (smoothstep)
    fn smoothstep(&self, t: f64) -> f64 {
        t * t * (3.0 - 2.0 * t)
    }

    /// Linear interpolation
    fn lerp(&self, a: f64, b: f64, t: f64) -> f64 {
        a + (b - a) * t
    }

    /// Smooth boundaries to prevent over-jaggedness
    fn smooth_boundaries(&mut self, plate_map: &mut TerrainMap<u16>, iterations: usize) {
        for iter in 0..iterations {
            let original_data = plate_map.data.clone();
            let mut changed = 0;

            for y in 1..plate_map.height - 1 {
                for x in 1..plate_map.width - 1 {
                    let idx = plate_map.get_index(x, y);
                    let current_plate = original_data[idx];

                    // Check if this is a boundary pixel
                    let neighbors = plate_map.get_neighbors(x, y);
                    let mut plate_counts = std::collections::HashMap::new();

                    for (nx, ny) in neighbors {
                        let neighbor_plate = original_data[plate_map.get_index(nx, ny)];
                        *plate_counts.entry(neighbor_plate).or_insert(0) += 1;
                    }

                    // If current plate is a minority among neighbors, consider changing it
                    if let Some((&most_common, &count)) = plate_counts.iter().max_by_key(|(_, &c)| c) {
                        // Only smooth very isolated pixels
                        if most_common != current_plate && count >= 6 {
                            plate_map.data[idx] = most_common;
                            changed += 1;
                        }
                    }
                }
            }

            if changed > 0 {
                println!("    Iteration {}/{}: smoothed {} pixels", iter + 1, iterations, changed);
            }
        }
    }
}

/// Convenience function to detect boundary type between two plates
pub fn detect_boundary_type(
    _plate_a_motion: (f64, f64),  // (direction_deg, speed_cm_yr)
    _plate_b_motion: (f64, f64),
) -> PlateInteraction {
    // TODO: Implement proper boundary type detection based on relative motion
    // For now, return a default
    PlateInteraction::Transform
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_config_creation() {
        let config = BoundaryRefinementConfig::with_seed(42)
            .with_noise(0.15, 4.0, 4)
            .with_smoothing(2);

        assert_eq!(config.seed, 42);
        assert_eq!(config.noise_scale, 0.15);
        assert_eq!(config.noise_amplitude, 4.0);
        assert_eq!(config.octaves, 4);
        assert_eq!(config.smoothing_iterations, 2);
    }

    #[test]
    fn test_boundary_identification() {
        // Create a simple 5x5 map with two plates
        let mut plate_map = TerrainMap::new(5, 5, 1u16);

        // Make right half plate 2
        for y in 0..5 {
            for x in 3..5 {
                let idx = plate_map.get_index(x, y);
                plate_map.data[idx] = 2;
            }
        }

        let config = BoundaryRefinementConfig::with_seed(42);
        let refiner = BoundaryRefiner::new(config);
        let boundaries = refiner.identify_boundaries(&plate_map);

        // Should have boundary pixels at x=2 and x=3
        assert!(boundaries.len() > 0, "Should have detected boundaries");
        assert!(boundaries.contains(&(2, 2)), "x=2 should be a boundary (next to plate 2)");
        assert!(boundaries.contains(&(3, 2)), "x=3 should be a boundary (next to plate 1)");

        // Note: Due to longitude wraparound in global projection, (0,0) might be detected as boundary
        // This is actually correct behavior for a global map where x=0 wraps to x=width-1
    }

    #[test]
    fn test_noise_determinism() {
        let config = BoundaryRefinementConfig::with_seed(123);
        let mut refiner1 = BoundaryRefiner::new(config.clone());
        let mut refiner2 = BoundaryRefiner::new(config);

        let noise1 = refiner1.multi_octave_noise(5.5, 10.2);
        let noise2 = refiner2.multi_octave_noise(5.5, 10.2);

        assert_eq!(noise1, noise2, "Noise should be deterministic with same seed");
    }
}
