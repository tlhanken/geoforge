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
use std::collections::{HashMap, VecDeque, HashSet};
use std::fs;
use std::path::Path;

// Constants
const PLANET_RADIUS_KM: f64 = 6371.0;
const PLANET_SURFACE_AREA_KM2: f64 = 4.0 * std::f64::consts::PI * PLANET_RADIUS_KM * PLANET_RADIUS_KM;
const MIN_SEED_DISTANCE: f64 = 100.0;
const MAX_SEED_ATTEMPTS: usize = 10000;

// Define pole threshold - pixels within this many degrees of poles are treated specially
const POLE_THRESHOLD_DEGREES: f64 = 5.0;

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

    /// Check if a pixel is at or near a pole
    fn is_near_pole(&self, y: usize) -> bool {
        let lat = self.y_to_lat(y);
        lat.abs() > (90.0 - POLE_THRESHOLD_DEGREES)
    }

    /// Get the pixel area in km² accounting for latitude distortion
    fn get_pixel_area_km2(&self, y: usize) -> f64 {
        let lat = self.y_to_lat(y);
        let lat_rad = lat.to_radians();
        
        // Area of a spherical rectangle
        let delta_lat = 180.0 / self.height as f64;
        let delta_lon = 360.0 / self.width as f64;
        
        let lat1_rad = (lat + delta_lat / 2.0).to_radians();
        let lat2_rad = (lat - delta_lat / 2.0).to_radians();
        
        let area = PLANET_RADIUS_KM * PLANET_RADIUS_KM * delta_lon.to_radians() * 
                   (lat1_rad.sin() - lat2_rad.sin()).abs();
        
        area
    }

    /// Fast approximation of spherical distance with pole handling
    fn fast_spherical_distance(&self, lon1: f64, lat1: f64, lon2: f64, lat2: f64) -> f64 {
        // Special case: if both points are very close to the same pole
        if lat1.abs() > 89.9 && lat2.abs() > 89.9 && lat1.signum() == lat2.signum() {
            // At the poles, all longitudes converge to the same point
            return ((lat2 - lat1).abs() / 180.0) * std::f64::consts::PI * PLANET_RADIUS_KM;
        }
        
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
        let lon = self.x_to_lon(x);
        let lat = self.y_to_lat(y);
        
        self.plate_seeds.iter().any(|seed| {
            let dist = self.fast_spherical_distance(lon, lat, seed.lon, seed.lat);
            dist < MIN_SEED_DISTANCE
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
                        
                        let lon = self.x_to_lon(x);
                        let lat = self.y_to_lat(y);
                        
                        let too_close = self.plate_seeds.iter().any(|seed| {
                            let dist = self.fast_spherical_distance(lon, lat, seed.lon, seed.lat);
                            dist < reduced_distance
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

                // Generate with latitude bias to account for pole convergence
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

    /// Generate plates using optimized Voronoi approach with pole handling
    pub fn generate_plates_voronoi(&mut self) {
        println!("Generating plates using Voronoi method with pole handling...");

        // Pre-calculate seed coordinates for faster lookup
        let seed_coords: Vec<(f64, f64, u16)> = self.plate_seeds
            .iter()
            .map(|seed| (seed.lon, seed.lat, seed.id))
            .collect();

        // Handle poles specially - assign entire pole to nearest plate
        let mut north_pole_plate = 0u16;
        let mut south_pole_plate = 0u16;
        
        // Find closest plate to each pole
        let mut min_north_dist = f64::INFINITY;
        let mut min_south_dist = f64::INFINITY;
        
        for seed in &self.plate_seeds {
            let north_dist = self.fast_spherical_distance(seed.lon, seed.lat, 0.0, 90.0);
            let south_dist = self.fast_spherical_distance(seed.lon, seed.lat, 0.0, -90.0);
            
            if north_dist < min_north_dist {
                min_north_dist = north_dist;
                north_pole_plate = seed.id;
            }
            
            if south_dist < min_south_dist {
                min_south_dist = south_dist;
                south_pole_plate = seed.id;
            }
        }

        for y in 0..self.height {
            let lat = self.y_to_lat(y);
            
            // Check if we're at a pole
            if y == 0 {
                // North pole - assign entire row to the closest plate
                for x in 0..self.width {
                    let idx = y * self.width + x;
                    self.plate_map[idx] = north_pole_plate;
                }
                continue;
            } else if y == self.height - 1 {
                // South pole - assign entire row to the closest plate
                for x in 0..self.width {
                    let idx = y * self.width + x;
                    self.plate_map[idx] = south_pole_plate;
                }
                continue;
            }
            
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

    /// Get neighbors of a pixel with pole handling
    fn get_neighbors(&self, x: usize, y: usize) -> Vec<(usize, usize)> {
        let mut neighbors = Vec::new();
        
        // Special handling for poles
        if y == 0 {
            // North pole - all pixels in the top row are neighbors
            for nx in 0..self.width {
                neighbors.push((nx, 0));
            }
            // Also add the row below
            for nx in 0..self.width {
                neighbors.push((nx, 1));
            }
            return neighbors;
        } else if y == self.height - 1 {
            // South pole - all pixels in the bottom row are neighbors
            for nx in 0..self.width {
                neighbors.push((nx, self.height - 1));
            }
            // Also add the row above
            for nx in 0..self.width {
                neighbors.push((nx, self.height - 2));
            }
            return neighbors;
        }
        
        // Regular 8-connected neighbors for non-pole pixels
        const DIRECTIONS: [(isize, isize); 8] = [
            (-1, -1), (-1, 0), (-1, 1),
            (0, -1),           (0, 1),
            (1, -1),  (1, 0),  (1, 1)
        ];
        
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
            
            neighbors.push((nx, ny));
        }
        
        // Near poles, add extra longitude connections
        if self.is_near_pole(y) {
            let lat = self.y_to_lat(y);
            let convergence_factor = lat.to_radians().cos();
            
            // Add connections to more distant longitude neighbors based on convergence
            if convergence_factor < 0.5 {
                let skip = ((1.0 / convergence_factor).max(2.0) as usize).min(self.width / 8);
                for i in (skip..self.width).step_by(skip) {
                    let nx = (x + i) % self.width;
                    neighbors.push((nx, y));
                }
            }
        }
        
        neighbors
    }

    /// Generate plates using region growing with pole handling
    pub fn generate_plates_region_growing(&mut self) -> Result<(), PlateError> {
        println!("Generating plates using region growing method with pole handling...");

        // Initialize with seeds
        self.plate_map.fill(0);
        let mut queue = VecDeque::with_capacity(self.width * self.height / 4);

        // Place seeds
        for seed in &self.plate_seeds {
            let idx = seed.y * self.width + seed.x;
            self.plate_map[idx] = seed.id;
            queue.push_back((seed.x, seed.y, seed.id));
        }

        // Grow regions
        while let Some((x, y, plate_id)) = queue.pop_front() {
            let neighbors = self.get_neighbors(x, y);
            
            for (nx, ny) in neighbors {
                let idx = ny * self.width + nx;

                // If unassigned, assign to current plate
                if self.plate_map[idx] == 0 {
                    self.plate_map[idx] = plate_id;

                    // Add controlled randomness for natural boundaries
                    // Less randomness near poles for more stable boundaries
                    let randomness_threshold = if self.is_near_pole(ny) {
                        0.05
                    } else {
                        0.1
                    };
                    
                    if self.rng.gen::<f64>() > randomness_threshold {
                        queue.push_back((nx, ny, plate_id));
                    }
                }
            }
        }

        // Ensure poles are fully assigned to single plates
        self.assign_poles_to_single_plate();

        Ok(())
    }

    /// Ensure each pole is assigned to a single plate
    fn assign_poles_to_single_plate(&mut self) {
        // North pole
        let mut north_pole_counts: HashMap<u16, usize> = HashMap::new();
        for x in 0..self.width {
            let plate_id = self.plate_map[x];
            if plate_id > 0 {
                *north_pole_counts.entry(plate_id).or_insert(0) += 1;
            }
        }
        
        if let Some((&dominant_plate, _)) = north_pole_counts.iter().max_by_key(|(_, &count)| count) {
            for x in 0..self.width {
                self.plate_map[x] = dominant_plate;
            }
        }
        
        // South pole
        let mut south_pole_counts: HashMap<u16, usize> = HashMap::new();
        let south_row = self.height - 1;
        for x in 0..self.width {
            let idx = south_row * self.width + x;
            let plate_id = self.plate_map[idx];
            if plate_id > 0 {
                *south_pole_counts.entry(plate_id).or_insert(0) += 1;
            }
        }
        
        if let Some((&dominant_plate, _)) = south_pole_counts.iter().max_by_key(|(_, &count)| count) {
            for x in 0..self.width {
                let idx = south_row * self.width + x;
                self.plate_map[idx] = dominant_plate;
            }
        }
    }

    /// Optimized boundary smoothing with pole handling
    pub fn smooth_boundaries(&mut self, iterations: usize) {
        println!("Smoothing boundaries ({} iterations) with pole handling...", iterations);

        for iter in 0..iterations {
            let original_map = self.plate_map.clone();
            
            for y in 0..self.height {
                // Skip poles - they should remain uniform
                if y == 0 || y == self.height - 1 {
                    continue;
                }
                
                for x in 0..self.width {
                    let idx = y * self.width + x;
                    let mut neighbor_counts = HashMap::new();

                    // Get neighbors with pole-aware logic
                    let neighbors = self.get_neighbors(x, y);
                    
                    for (nx, ny) in neighbors {
                        let neighbor_idx = ny * self.width + nx;
                        let neighbor_plate = original_map[neighbor_idx];
                        *neighbor_counts.entry(neighbor_plate).or_insert(0) += 1;
                    }

                    // Find most common neighbor
                    if let Some((&most_common, &max_count)) = neighbor_counts
                        .iter()
                        .max_by_key(|(_, &count)| count)
                    {
                        // Adjust threshold based on latitude
                        let threshold = if self.is_near_pole(y) {
                            3  // Lower threshold near poles
                        } else {
                            4  // Normal threshold
                        };
                        
                        if max_count >= threshold {
                            self.plate_map[idx] = most_common;
                        }
                    }
                }
            }

            // Ensure poles remain assigned to single plates after smoothing
            self.assign_poles_to_single_plate();

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

    /// Generate comprehensive statistics with area-weighted calculations
    pub fn get_plate_stats(&self) -> HashMap<u16, PlateStats> {
        let mut plate_areas = HashMap::new();
        
        // Calculate area for each plate accounting for latitude
        for y in 0..self.height {
            let pixel_area = self.get_pixel_area_km2(y);
            
            for x in 0..self.width {
                let idx = y * self.width + x;
                let pixel = self.plate_map[idx];
                
                if pixel > 0 {
                    let entry = plate_areas.entry(pixel).or_insert((0, 0.0));
                    entry.0 += 1;  // pixel count
                    entry.1 += pixel_area;  // area
                }
            }
        }

        let total_pixels = self.width * self.height;

        // Generate statistics
        let mut stats = HashMap::new();
        for (&plate_id, &(pixels, area_km2)) in &plate_areas {
            if let Some(seed) = self.plate_seeds.iter().find(|s| s.id == plate_id) {
                stats.insert(
                    plate_id,
                    PlateStats {
                        pixels,
                        percentage: (area_km2 / PLANET_SURFACE_AREA_KM2) * 100.0,
                        area_km2: area_km2 as u64,
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
        
        // Always export binary
        self.export_binary(output_dir, &format!("{}.bin", base_name))?;
        
        #[cfg(feature = "export-png")]
        self.export_png(output_dir, &format!("{}.png", base_name))?;
        
        #[cfg(feature = "export-tiff")]
        self.export_geotiff(output_dir, &format!("{}.tiff", base_name))?;
        
        Ok(())
    }

    #[cfg(feature = "export-png")]
    pub fn export_png(&self, output_dir: &str, filename: &str) -> Result<(), Box<dyn std::error::Error>> {
        use image::{ImageBuffer, Rgb};
        
        let mut img = ImageBuffer::new(self.width as u32, self.height as u32);
        
        // Generate a color for each plate
        let mut rng = StdRng::seed_from_u64(42); // Consistent colors
        let colors: Vec<[u8; 3]> = (0..=self.num_plates).map(|_| {
            [rng.gen(), rng.gen(), rng.gen()]
        }).collect();
        
        // Draw the image
        for (y, row) in self.plate_map.chunks(self.width).enumerate() {
            for (x, &plate_id) in row.iter().enumerate() {
                let color = colors[plate_id as usize];
                img.put_pixel(x as u32, y as u32, Rgb(color));
            }
        }
        
        let path = Path::new(output_dir).join(filename);
        img.save(path)?;
        Ok(())
    }

    #[cfg(feature = "export-tiff")]
    pub fn export_geotiff(&self, output_dir: &str, filename: &str) -> Result<(), Box<dyn std::error::Error>> {
        use gdal::Dataset;
        use gdal::Driver;
        use gdal::spatial_ref::SpatialRef;
        
        let driver = Driver::get("GTiff")?;
        let path = Path::new(output_dir).join(filename);
        let mut ds = driver.create_with_band_type::<u16, _>(
            path,
            self.width as isize,
            self.height as isize,
            1
        )?;
        
        // Set geotransform and projection
        ds.set_geo_transform(&self.get_geotransform())?;
        let srs = SpatialRef::from_epsg(4326)?;  // WGS84
        ds.set_spatial_ref(&srs)?;
        
        // Write data
        let band = ds.rasterband(1)?;
        band.write_band_as_slice(0, 0, self.width as usize, self.height as usize, &self.plate_map)?;
        
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

        // Validate poles are uniform
        let north_pole_plate = self.plate_map[0];
        for x in 1..self.width {
            if self.plate_map[x] != north_pole_plate {
                return Err(PlateError::InvalidParameters("North pole not uniform".to_string()));
            }
        }
        
        let south_idx_start = (self.height - 1) * self.width;
        let south_pole_plate = self.plate_map[south_idx_start];
        for x in 1..self.width {
            if self.plate_map[south_idx_start + x] != south_pole_plate {
                return Err(PlateError::InvalidParameters("South pole not uniform".to_string()));
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
    fn test_pole_detection() {
        let generator = TectonicPlateGenerator::new(3600, 1800, 5).unwrap();
        
        // Test pole detection
        assert!(generator.is_near_pole(0));  // North pole
        assert!(generator.is_near_pole(1799));  // South pole
        assert!(!generator.is_near_pole(900));  // Equator
    }

    #[test]
    fn test_distance_calculations() {
        let generator = TectonicPlateGenerator::new(3600, 1800, 5).unwrap();

        // Test pole distance
        let pole_dist = generator.fast_spherical_distance(0.0, 89.95, 180.0, 89.95);
        assert!(pole_dist < 20.0);  // Should be very small at poles

        // 1 degree latitude at any longitude
        let expected_lat = std::f64::consts::PI * PLANET_RADIUS_KM / 180.0;
        let dist_lat = generator.fast_spherical_distance(0.0, 0.0, 0.0, 1.0);
        assert!((dist_lat - expected_lat).abs() < 0.5);
    }

    #[test]
    fn test_pole_uniformity() {
        let mut generator = TectonicPlateGenerator::with_seed(3600, 1800, 5, 123).unwrap();
        generator.generate("region_growing", false).unwrap();
        
        // Check north pole uniformity
        let north_pole_plate = generator.plate_map[0];
        for x in 0..generator.width {
            assert_eq!(generator.plate_map[x], north_pole_plate, 
                      "North pole not uniform at x={}", x);
        }
        
        // Check south pole uniformity
        let south_idx = (generator.height - 1) * generator.width;
        let south_pole_plate = generator.plate_map[south_idx];
        for x in 0..generator.width {
            assert_eq!(generator.plate_map[south_idx + x], south_pole_plate,
                      "South pole not uniform at x={}", x);
        }
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
