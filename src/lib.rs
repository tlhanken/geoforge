//! Geoforge - Realistic geological and climate modeling for procedural world generation
//! 
//! This library provides tools for generating realistic geological features, climate patterns,
//! and biomes for procedural world generation. It starts with tectonic plate simulation
//! and builds up through geological domains, elevation, climate, and biomes.
//! 
//! # Quick Start
//! 
//! ```rust
//! use geoforge::TectonicPlateGenerator;
//! 
//! // Generate a world with 15 tectonic plates at 0.2° resolution
//! let mut generator = TectonicPlateGenerator::with_seed(1800, 900, 15, 42)?;
//! let plate_map = generator.generate("region_growing", true)?;
//! 
//! // Get statistics about the generated plates
//! let stats = generator.get_plate_stats();
//! for (plate_id, stat) in stats {
//!     println!("Plate {}: {:.1}% of surface", plate_id, stat.percentage);
//! }
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use rand::prelude::*;
use rand::rngs::StdRng;
use std::collections::{HashMap, VecDeque};
use std::fs;
use std::path::Path;

// Constants
const PLANET_RADIUS_KM: f64 = 6371.0;
const PLANET_SURFACE_AREA_KM2: f64 = 4.0 * std::f64::consts::PI * PLANET_RADIUS_KM * PLANET_RADIUS_KM;
const MIN_SEED_DISTANCE: f64 = 100.0;
const MAX_SEED_ATTEMPTS: usize = 10000;

/// Represents a tectonic plate seed point with motion information
#[derive(Debug, Clone)]
pub struct PlateSeed {
    pub id: u16,
    pub x: usize,
    pub y: usize,
    pub lon: f64,
    pub lat: f64,
    pub motion_direction: f64, // degrees
    pub motion_speed: f64,     // cm/year
}

/// Statistics about a generated tectonic plate
#[derive(Debug, Clone)]
pub struct PlateStats {
    pub pixels: usize,
    pub percentage: f64,
    pub area_km2: u64,
    pub seed: PlateSeed,
}

/// Errors that can occur during plate generation
#[derive(Debug)]
pub enum PlateError {
    SeedPlacementFailed,
    InvalidMethod(String),
    InvalidParameters(String),
}

impl std::fmt::Display for PlateError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PlateError::SeedPlacementFailed => write!(f, "Failed to place seeds with minimum distance"),
            PlateError::InvalidMethod(method) => write!(f, "Invalid generation method: {}", method),
            PlateError::InvalidParameters(msg) => write!(f, "Invalid parameters: {}", msg),
        }
    }
}

impl std::error::Error for PlateError {}

/// Generator for realistic tectonic plates using various algorithms
pub struct TectonicPlateGenerator {
    width: usize,
    height: usize,
    num_plates: usize,
    plate_map: Vec<u16>,
    plate_seeds: Vec<PlateSeed>,
    rng: StdRng,
    current_seed: u64,
}

impl TectonicPlateGenerator {
    /// Create a new tectonic plate generator with random seed
    pub fn new(width: usize, height: usize, num_plates: usize) -> Result<Self, PlateError> {
        Self::with_seed(width, height, num_plates, rand::random())
    }

    /// Create a new tectonic plate generator with specific seed
    pub fn with_seed(width: usize, height: usize, num_plates: usize, seed: u64) -> Result<Self, PlateError> {
        if width == 0 || height == 0 {
            return Err(PlateError::InvalidParameters("Width and height must be > 0".to_string()));
        }
        if num_plates == 0 || num_plates > u16::MAX as usize {
            return Err(PlateError::InvalidParameters("Invalid number of plates".to_string()));
        }

        Ok(Self {
            width,
            height,
            num_plates,
            plate_map: vec![0; width * height],
            plate_seeds: Vec::with_capacity(num_plates),
            rng: StdRng::seed_from_u64(seed),
            current_seed: seed,
        })
    }

    /// Get the current seed
    pub fn get_seed(&self) -> u64 {
        self.current_seed
    }

    /// Reset the generator with a new seed
    pub fn set_seed(&mut self, seed: u64) {
        self.current_seed = seed;
        self.rng = StdRng::seed_from_u64(seed);
        self.plate_map.fill(0);
        self.plate_seeds.clear();
    }

    /// Convert pixel coordinates to geographic coordinates
    fn x_to_lon(&self, x: usize) -> f64 {
        (x as f64 / self.width as f64) * 360.0 - 180.0
    }

    fn y_to_lat(&self, y: usize) -> f64 {
        90.0 - (y as f64 / self.height as f64) * 180.0
    }

    /// Fast approximation of spherical distance for small distances
    /// More accurate than euclidean, much faster than full haversine
    fn fast_spherical_distance(&self, lon1: f64, lat1: f64, lon2: f64, lat2: f64) -> f64 {
        // Handle longitude wraparound
        let mut d_lon = (lon2 - lon1).abs();
        if d_lon > 180.0 {
            d_lon = 360.0 - d_lon;
        }

        let d_lat = lat2 - lat1;
        let lat_avg = (lat1 + lat2) * 0.5;

        // Convert degree differences to radians
        let d_lon_rad = d_lon.to_radians();
        let d_lat_rad = d_lat.to_radians();
        let lat_avg_rad = lat_avg.to_radians();

        // Approximate distance using equirectangular projection
        let x = d_lon_rad * lat_avg_rad.cos();
        let y = d_lat_rad;
        (x * x + y * y).sqrt() * PLANET_RADIUS_KM
}

    /// Check if a point is too close to existing seeds
    fn too_close_to_existing_seed(&self, x: usize, y: usize) -> bool {
        self.plate_seeds.iter().any(|seed| {
            let dx = x as f64 - seed.x as f64;
            let dy = y as f64 - seed.y as f64;
            (dx * dx + dy * dy).sqrt() < MIN_SEED_DISTANCE
        })
    }

    /// Generate random seed points for plates with improved distribution
    pub fn generate_seeds(&mut self) -> Result<(), PlateError> {
        self.plate_seeds.clear();

        for i in 0..self.num_plates {
            let mut attempts = 0;
            let (x, y) = loop {
                attempts += 1;
                
                if attempts > MAX_SEED_ATTEMPTS {
                    // Fallback: reduce minimum distance requirement
                    let reduced_distance = MIN_SEED_DISTANCE * 0.5;
                    let mut fallback_found = false;
                    let mut fallback_coords = (0, 0);
                    
                    for _attempt in 0..1000 {
                        let x = self.rng.gen_range(0..self.width);
                        let y = self.rng.gen_range(0..self.height);
                        
                        let too_close = self.plate_seeds.iter().any(|seed| {
                            let dx = x as f64 - seed.x as f64;
                            let dy = y as f64 - seed.y as f64;
                            (dx * dx + dy * dy).sqrt() < reduced_distance
                        });
                        
                        if !too_close {
                            fallback_coords = (x, y);
                            fallback_found = true;
                            break;
                        }
                    }
                    
                    if fallback_found {
                        break fallback_coords;
                    } else {
                        return Err(PlateError::SeedPlacementFailed);
                    }
                }

                let x = self.rng.gen_range(0..self.width);
                let y = self.rng.gen_range(0..self.height);

                if !self.too_close_to_existing_seed(x, y) {
                    break (x, y);
                }
            };

            let lon = self.x_to_lon(x);
            let lat = self.y_to_lat(y);

            self.plate_seeds.push(PlateSeed {
                id: (i + 1) as u16,
                x,
                y,
                lon,
                lat,
                motion_direction: self.rng.gen_range(0.0..360.0),
                motion_speed: self.rng.gen_range(2.0..10.0),
            });
        }

        println!("Generated {} plate seeds", self.plate_seeds.len());
        Ok(())
    }

    /// Generate plates using optimized Voronoi approach
    pub fn generate_plates_voronoi(&mut self) {
        println!("Generating plates using Voronoi method...");

        // Pre-calculate seed coordinates for faster lookup
        let seed_coords: Vec<(f64, f64, u16)> = self.plate_seeds
            .iter()
            .map(|seed| (seed.lon, seed.lat, seed.id))
            .collect();

        for y in 0..self.height {
            let lat = self.y_to_lat(y);
            
            for x in 0..self.width {
                let lon = self.x_to_lon(x);
                let idx = y * self.width + x;

                let mut nearest_plate = 1;
                let mut min_distance = f64::INFINITY;

                // Find closest seed using fast distance approximation
                for &(seed_lon, seed_lat, seed_id) in &seed_coords {
                    let distance = self.fast_spherical_distance(lon, lat, seed_lon, seed_lat);
                    if distance < min_distance {
                        min_distance = distance;
                        nearest_plate = seed_id;
                    }
                }

                self.plate_map[idx] = nearest_plate;
            }

            if y % 100 == 0 {
                println!("Progress: {:.1}%", (y as f64 / self.height as f64) * 100.0);
            }
        }
    }

    /// Generate plates using region growing with deterministic randomness
    pub fn generate_plates_region_growing(&mut self) -> Result<(), PlateError> {
        println!("Generating plates using region growing method...");

        // Initialize with seeds
        self.plate_map.fill(0);
        let mut queue = VecDeque::with_capacity(self.width * self.height / 4);

        // Place seeds
        for seed in &self.plate_seeds {
            let idx = seed.y * self.width + seed.x;
            self.plate_map[idx] = seed.id;
            queue.push_back((seed.x, seed.y, seed.id));
        }

        // Directions for 8-connected neighbors
        const DIRECTIONS: [(isize, isize); 8] = [
            (-1, -1), (-1, 0), (-1, 1),
            (0, -1),           (0, 1),
            (1, -1),  (1, 0),  (1, 1)
        ];

        // Grow regions
        while let Some((x, y, plate_id)) = queue.pop_front() {
            for &(dx, dy) in &DIRECTIONS {
                let mut nx = x as isize + dx;
                let ny = y as isize + dy;

                // Handle longitude wraparound
                if nx < 0 {
                    nx = self.width as isize - 1;
                } else if nx >= self.width as isize {
                    nx = 0;
                }

                // Handle latitude bounds
                if ny < 0 || ny >= self.height as isize {
                    continue;
                }

                let nx = nx as usize;
                let ny = ny as usize;
                let idx = ny * self.width + nx;

                // If unassigned, assign to current plate
                if self.plate_map[idx] == 0 {
                    self.plate_map[idx] = plate_id;

                    // Add controlled randomness for natural boundaries
                    if self.rng.gen::<f64>() > 0.1 {
                        queue.push_back((nx, ny, plate_id));
                    }
                }
            }
        }

        Ok(())
    }

    /// Optimized boundary smoothing with wraparound handling
    pub fn smooth_boundaries(&mut self, iterations: usize) {
        println!("Smoothing boundaries ({} iterations)...", iterations);

        for iter in 0..iterations {
            let original_map = &self.plate_map.clone(); // Only clone once
            
            for y in 0..self.height {
                for x in 0..self.width {
                    let idx = y * self.width + x;
                    let mut neighbor_counts = HashMap::new();

                    // Count neighbors with proper wraparound
                    for dy in -1isize..=1 {
                        for dx in -1isize..=1 {
                            if dx == 0 && dy == 0 { continue; }

                            let mut nx = x as isize + dx;
                            let ny = y as isize + dy;

                            // Handle longitude wraparound
                            if nx < 0 {
                                nx = self.width as isize - 1;
                            } else if nx >= self.width as isize {
                                nx = 0;
                            }

                            // Handle latitude bounds
                            if ny < 0 || ny >= self.height as isize {
                                continue;
                            }

                            let neighbor_idx = ny as usize * self.width + nx as usize;
                            let neighbor_plate = original_map[neighbor_idx];
                            *neighbor_counts.entry(neighbor_plate).or_insert(0) += 1;
                        }
                    }

                    // Find most common neighbor
                    if let Some((&most_common, &max_count)) = neighbor_counts
                        .iter()
                        .max_by_key(|(_, &count)| count)
                    {
                        // Only change if there's strong consensus (>= 4 out of 8 neighbors)
                        if max_count >= 4 {
                            self.plate_map[idx] = most_common;
                        }
                    }
                }
            }

            println!("Completed smoothing iteration {}/{}", iter + 1, iterations);
        }
    }

    /// Main generation function with better error handling
    pub fn generate(&mut self, method: &str, smooth: bool) -> Result<&Vec<u16>, PlateError> {
        println!("Generating {} tectonic plates using {} method... (seed: {})", 
                 self.num_plates, method, self.current_seed);

        // Step 1: Generate seed points
        self.generate_seeds()?;

        // Step 2: Generate plate boundaries
        match method {
            "voronoi" => self.generate_plates_voronoi(),
            "region_growing" => self.generate_plates_region_growing()?,
            _ => return Err(PlateError::InvalidMethod(method.to_string())),
        }

        // Step 3: Smooth boundaries if requested
        if smooth {
            self.smooth_boundaries(3);
        }

        println!("Tectonic plate generation complete!");
        Ok(&self.plate_map)
    }

    /// Generate comprehensive statistics
    pub fn get_plate_stats(&self) -> HashMap<u16, PlateStats> {
        let mut plate_sizes = HashMap::new();
        
        // Count pixels for each plate
        for &pixel in &self.plate_map {
            if pixel > 0 {
                *plate_sizes.entry(pixel).or_insert(0) += 1;
            }
        }

        let total_pixels = self.width * self.height;
        let km_per_pixel = PLANET_SURFACE_AREA_KM2 as f64 / total_pixels as f64;

        // Generate statistics
        let mut stats = HashMap::new();
        for (&plate_id, &pixels) in &plate_sizes {
            if let Some(seed) = self.plate_seeds.iter().find(|s| s.id == plate_id) {
                stats.insert(
                    plate_id,
                    PlateStats {
                        pixels,
                        percentage: (pixels as f64 / total_pixels as f64) * 100.0,
                        area_km2: (pixels as f64 * km_per_pixel) as u64,
                        seed: seed.clone(),
                    },
                );
            }
        }

        stats
    }

    /// Get plate map dimensions and data
    pub fn get_plate_data(&self) -> (usize, usize, &Vec<u16>, &Vec<PlateSeed>) {
        (self.width, self.height, &self.plate_map, &self.plate_seeds)
    }

    /// Export as little-endian bytes for file writing
    pub fn export_raw_u16_le(&self) -> Vec<u8> {
        self.plate_map
            .iter()
            .flat_map(|&value| value.to_le_bytes())
            .collect()
    }

    /// Export raw binary data to file
    pub fn export_binary(&self, output_dir: &str, filename: &str) -> Result<(), Box<dyn std::error::Error>> {
        fs::create_dir_all(output_dir)?;
        let path = Path::new(output_dir).join(filename);
        let raw_bytes = self.export_raw_u16_le();
        fs::write(path, raw_bytes)?;
        Ok(())
    }

    /// Export all available formats
    pub fn export_all(&self, output_dir: &str, base_name: &str) -> Result<(), Box<dyn std::error::Error>> {
        fs::create_dir_all(output_dir)?;
        self.export_binary(output_dir, &format!("{}.bin", base_name))?;
        Ok(())
    }

    /// Get GeoTIFF geotransform parameters
    pub fn get_geotransform(&self) -> [f64; 6] {
        [
            -180.0,                              // Top-left X (longitude)
            360.0 / self.width as f64,          // Pixel width (degrees per pixel)
            0.0,                                // Row rotation (typically 0)
            90.0,                               // Top-left Y (latitude)
            0.0,                                // Column rotation (typically 0)
            -180.0 / self.height as f64,        // Pixel height (negative = north up)
        ]
    }

    /// Validate the generated plate map for consistency
    pub fn validate(&self) -> Result<(), PlateError> {
        // Check all pixels are assigned
        if self.plate_map.iter().any(|&p| p == 0) {
            return Err(PlateError::InvalidParameters("Unassigned pixels found".to_string()));
        }

        // Check all plate IDs are valid
        let max_plate_id = self.num_plates as u16;
        if self.plate_map.iter().any(|&p| p > max_plate_id) {
            return Err(PlateError::InvalidParameters("Invalid plate ID found".to_string()));
        }

        // Check all plates have at least some pixels
        let mut plate_counts = vec![0; self.num_plates + 1];
        for &plate_id in &self.plate_map {
            if plate_id > 0 && plate_id <= max_plate_id {
                plate_counts[plate_id as usize] += 1;
            }
        }

        for (i, &count) in plate_counts.iter().enumerate().skip(1) {
            if count == 0 {
                return Err(PlateError::InvalidParameters(format!("Plate {} has no pixels", i)));
            }
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_generator() {
        let generator = TectonicPlateGenerator::new(1800, 900, 15);
        assert!(generator.is_ok());

        let invalid = TectonicPlateGenerator::new(0, 900, 15);
        assert!(invalid.is_err());
    }

    #[test]
    fn test_coordinate_conversion() {
        let generator = TectonicPlateGenerator::new(1800, 900, 5).unwrap();

        // Test coordinate conversion
        assert_eq!(generator.x_to_lon(0), -180.0);
        assert_eq!(generator.x_to_lon(1800), 180.0);
        assert_eq!(generator.y_to_lat(0), 90.0);
        assert_eq!(generator.y_to_lat(900), -90.0);
    }

    #[test]
    fn test_distance_calculations() {
        let generator = TectonicPlateGenerator::new(3600, 1800, 5).unwrap();

        // 1 degree latitude at any longitude
        let expected_lat = std::f64::consts::PI * PLANET_RADIUS_KM / 180.0; // ≈ 111.19 km
        let dist_lat = generator.fast_spherical_distance(0.0, 0.0, 0.0, 1.0);
        assert!((dist_lat - expected_lat).abs() < 0.5);

        // 1 degree longitude at equator
        let expected_lon = std::f64::consts::PI * PLANET_RADIUS_KM / 180.0; // ≈ 111.19 km
        let dist_lon = generator.fast_spherical_distance(0.0, 0.0, 1.0, 0.0);
        assert!((dist_lon - expected_lon).abs() < 0.5);

        // Test wraparound (crossing date line)
        let dist_wrap = generator.fast_spherical_distance(-179.95, 0.0, 179.95, 0.0);
        let dist_no_wrap = generator.fast_spherical_distance(-0.05, 0.0, 0.05, 0.0);  // Same 0.1° distance
        let diff = (dist_wrap - dist_no_wrap).abs();
        println!("dist_wrap: {}, dist_no_wrap: {}, difference: {}", 
                 dist_wrap, dist_no_wrap, diff);
        assert!(diff < 0.5,
                "Distance difference too large: wrap={}, no_wrap={}, diff={}", 
                dist_wrap, dist_no_wrap, diff);
    }

    #[test]
    fn test_deterministic_generation() {
        let mut gen1 = TectonicPlateGenerator::with_seed(3600, 1800, 3, 123).unwrap();
        let mut gen2 = TectonicPlateGenerator::with_seed(3600, 1800, 3, 123).unwrap();

        let result1 = gen1.generate("region_growing", false);
        let result2 = gen2.generate("region_growing", false);

        assert!(result1.is_ok() && result2.is_ok());
        // With same seed, results should be identical
        assert_eq!(result1.unwrap(), result2.unwrap());
    }

    #[test]
    fn test_validation() {
        let mut generator = TectonicPlateGenerator::new(3600, 1800, 3).unwrap();
        let result = generator.generate("voronoi", false);
        assert!(result.is_ok());
        assert!(generator.validate().is_ok());
    }
}
