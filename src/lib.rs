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
use std::collections::{HashMap, VecDeque, HashSet, BinaryHeap};
use std::cmp::Ordering;
use std::fs;
use std::path::Path;

// Constants
const PLANET_RADIUS_KM: f64 = 6371.0;
const PLANET_SURFACE_AREA_KM2: f64 = 4.0 * std::f64::consts::PI * PLANET_RADIUS_KM * PLANET_RADIUS_KM;
const MIN_SEED_DISTANCE: f64 = 100.0;
const MAX_SEED_ATTEMPTS: usize = 10000;

// 3D point on unit sphere
#[derive(Debug, Clone, Copy)]
struct Point3D {
    x: f64,
    y: f64,
    z: f64,
}

impl Point3D {
    fn new(x: f64, y: f64, z: f64) -> Self {
        let len = (x * x + y * y + z * z).sqrt();
        Point3D {
            x: x / len,
            y: y / len,
            z: z / len,
        }
    }

    fn from_lat_lon(lat: f64, lon: f64) -> Self {
        let lat_rad = lat.to_radians();
        let lon_rad = lon.to_radians();
        
        Point3D::new(
            lat_rad.cos() * lon_rad.cos(),
            lat_rad.cos() * lon_rad.sin(),
            lat_rad.sin()
        )
    }

    fn to_lat_lon(&self) -> (f64, f64) {
        let lat = self.z.asin().to_degrees();
        let lon = self.y.atan2(self.x).to_degrees();
        (lat, lon)
    }

    fn geodesic_distance(&self, other: &Point3D) -> f64 {
        // Dot product gives cos of angle between unit vectors
        let dot = self.x * other.x + self.y * other.y + self.z * other.z;
        // Clamp to avoid numerical errors
        let dot = dot.max(-1.0).min(1.0);
        dot.acos()
    }

    fn geodesic_distance_km(&self, other: &Point3D) -> f64 {
        self.geodesic_distance(other) * PLANET_RADIUS_KM
    }
}

// Growth frontier item for priority queue
#[derive(Clone)]
struct GrowthPoint {
    point: Point3D,
    plate_id: u16,
    distance: f64,
}

impl PartialEq for GrowthPoint {
    fn eq(&self, other: &Self) -> bool {
        self.distance == other.distance
    }
}

impl Eq for GrowthPoint {}

impl PartialOrd for GrowthPoint {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        // Reverse ordering for min-heap behavior
        other.distance.partial_cmp(&self.distance)
    }
}

impl Ord for GrowthPoint {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap_or(Ordering::Equal)
    }
}

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
    point_3d: Point3D,
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
    // Cache 3D points for each pixel
    pixel_points: Vec<Point3D>,
    // Variable plate sizes
    plate_size_targets: Vec<f64>,
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

        let mut generator = Self {
            width,
            height,
            num_plates,
            plate_map: vec![0; width * height],
            plate_seeds: Vec::with_capacity(num_plates),
            rng: StdRng::seed_from_u64(seed),
            current_seed: seed,
            pixel_points: Vec::with_capacity(width * height),
            plate_size_targets: Vec::new(),
        };

        // Pre-calculate 3D points for all pixels
        generator.calculate_pixel_points();
        
        Ok(generator)
    }

    /// Pre-calculate 3D points for all pixels
    fn calculate_pixel_points(&mut self) {
        self.pixel_points.clear();
        
        for y in 0..self.height {
            let lat = self.y_to_lat(y);
            for x in 0..self.width {
                let lon = self.x_to_lon(x);
                self.pixel_points.push(Point3D::from_lat_lon(lat, lon));
            }
        }
    }

    /// Generate varied plate sizes using power law distribution
    fn generate_plate_sizes(&mut self) {
        self.plate_size_targets.clear();
        
        // Generate sizes with power law distribution (like real Earth)
        let mut sizes: Vec<f64> = Vec::new();
        let alpha = 1.5; // Power law exponent
        
        for i in 0..self.num_plates {
            // Power law: larger index = smaller plate
            let size = 1.0 / ((i + 1) as f64).powf(alpha);
            sizes.push(size);
        }
        
        // Shuffle to randomize which plates get which sizes
        sizes.shuffle(&mut self.rng);
        
        // Normalize so they sum to 1.0
        let sum: f64 = sizes.iter().sum();
        for size in sizes {
            self.plate_size_targets.push(size / sum);
        }
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
        self.plate_size_targets.clear();
    }

    /// Convert pixel coordinates to geographic coordinates
    fn x_to_lon(&self, x: usize) -> f64 {
        (x as f64 / self.width as f64) * 360.0 - 180.0
    }

    fn y_to_lat(&self, y: usize) -> f64 {
        90.0 - (y as f64 / self.height as f64) * 180.0
    }

    /// Get the pixel area in km² accounting for latitude distortion
    fn get_pixel_area_km2(&self, y: usize) -> f64 {
        let lat = self.y_to_lat(y);
        
        // Area of a spherical rectangle
        let delta_lat = 180.0 / self.height as f64;
        let delta_lon = 360.0 / self.width as f64;
        
        let lat1_rad = (lat + delta_lat / 2.0).to_radians();
        let lat2_rad = (lat - delta_lat / 2.0).to_radians();
        
        let area = PLANET_RADIUS_KM * PLANET_RADIUS_KM * delta_lon.to_radians() * 
                   (lat1_rad.sin() - lat2_rad.sin()).abs();
        
        area
    }

    /// Get pixel index from coordinates
    fn get_pixel_index(&self, x: usize, y: usize) -> usize {
        y * self.width + x
    }

    /// Get 3D point for pixel
    fn get_pixel_point(&self, x: usize, y: usize) -> &Point3D {
        &self.pixel_points[self.get_pixel_index(x, y)]
    }

    /// Generate random seed points for plates with improved distribution
    pub fn generate_seeds(&mut self) -> Result<(), PlateError> {
        self.plate_seeds.clear();
        
        // Generate plate sizes
        self.generate_plate_sizes();

        for i in 0..self.num_plates {
            let mut attempts = 0;
            let (x, y, point_3d) = loop {
                attempts += 1;
                
                if attempts > MAX_SEED_ATTEMPTS {
                    return Err(PlateError::SeedPlacementFailed);
                }

                // Generate random point on sphere with uniform distribution
                // Using inverse transform sampling for uniform sphere distribution
                let u = self.rng.gen::<f64>();
                let v = self.rng.gen::<f64>();
                
                let theta = 2.0 * std::f64::consts::PI * u;
                let phi = (2.0 * v - 1.0).acos();
                
                let x_3d = phi.sin() * theta.cos();
                let y_3d = phi.sin() * theta.sin();
                let z_3d = phi.cos();
                
                let point_3d = Point3D::new(x_3d, y_3d, z_3d);
                let (lat, lon) = point_3d.to_lat_lon();
                
                // Convert to pixel coordinates
                let x = ((lon + 180.0) / 360.0 * self.width as f64) as usize;
                let y = ((90.0 - lat) / 180.0 * self.height as f64) as usize;
                
                // Clamp to valid range
                let x = x.min(self.width - 1);
                let y = y.min(self.height - 1);
                
                // Check minimum distance to existing seeds
                let too_close = self.plate_seeds.iter().any(|seed| {
                    seed.point_3d.geodesic_distance_km(&point_3d) < MIN_SEED_DISTANCE
                });
                
                if !too_close {
                    break (x, y, point_3d);
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
                point_3d,
            });
        }

        println!("Generated {} plate seeds with varied sizes", self.plate_seeds.len());
        Ok(())
    }

    /// Generate plates using spherical Voronoi
    pub fn generate_plates_voronoi(&mut self) {
        println!("Generating plates using spherical Voronoi method...");

        for y in 0..self.height {
            for x in 0..self.width {
                let idx = self.get_pixel_index(x, y);
                let point = self.get_pixel_point(x, y);

                let mut nearest_plate = 1;
                let mut min_distance = f64::INFINITY;

                // Find closest seed using geodesic distance
                for seed in &self.plate_seeds {
                    let distance = point.geodesic_distance(&seed.point_3d);
                    if distance < min_distance {
                        min_distance = distance;
                        nearest_plate = seed.id;
                    }
                }

                self.plate_map[idx] = nearest_plate;
            }

            if y % 100 == 0 {
                println!("Progress: {:.1}%", (y as f64 / self.height as f64) * 100.0);
            }
        }
    }

    /// Generate plates using spherical region growing
    pub fn generate_plates_region_growing(&mut self) -> Result<(), PlateError> {
        println!("Generating plates using spherical region growing...");

        // Initialize
        self.plate_map.fill(0);
        let mut assigned_pixels = HashSet::new();
        let mut plate_pixel_counts = vec![0usize; self.num_plates + 1];
        let total_pixels = self.width * self.height;
        
        // Priority queue for growth frontier
        let mut frontier = BinaryHeap::new();
        
        // Place seeds and initialize frontier
        for (i, seed) in self.plate_seeds.iter().enumerate() {
            let idx = self.get_pixel_index(seed.x, seed.y);
            self.plate_map[idx] = seed.id;
            assigned_pixels.insert(idx);
            plate_pixel_counts[seed.id as usize] = 1;
            
            // Add neighbors to frontier
            self.add_neighbors_to_frontier(seed.x, seed.y, seed.id, &seed.point_3d, 
                                         &mut frontier, &assigned_pixels);
        }
        
        // Grow regions
        let mut pixels_assigned = self.num_plates;
        
        while let Some(growth_point) = frontier.pop() {
            // Find pixel closest to this 3D point
            let (lat, lon) = growth_point.point.to_lat_lon();
            let x = ((lon + 180.0) / 360.0 * self.width as f64) as usize;
            let y = ((90.0 - lat) / 180.0 * self.height as f64) as usize;
            
            // Clamp to valid range
            let x = x.min(self.width - 1);
            let y = y.min(self.height - 1);
            
            let idx = self.get_pixel_index(x, y);
            
            // Skip if already assigned
            if assigned_pixels.contains(&idx) {
                continue;
            }
            
            // Check if this plate has reached its target size
            let plate_idx = growth_point.plate_id as usize;
            let current_fraction = plate_pixel_counts[plate_idx] as f64 / total_pixels as f64;
            let target_fraction = self.plate_size_targets[plate_idx - 1];
            
            // Add randomness but bias based on target size
            let growth_probability = if current_fraction < target_fraction {
                0.95 // High probability if under target
            } else {
                // Decreasing probability as we exceed target
                0.95 * (target_fraction / current_fraction).min(1.0)
            };
            
            if self.rng.gen::<f64>() < growth_probability {
                // Assign pixel
                self.plate_map[idx] = growth_point.plate_id;
                assigned_pixels.insert(idx);
                plate_pixel_counts[plate_idx] += 1;
                pixels_assigned += 1;
                
                // Add neighbors to frontier
                self.add_neighbors_to_frontier(x, y, growth_point.plate_id, &growth_point.point,
                                             &mut frontier, &assigned_pixels);
            }
            
            // Progress update
            if pixels_assigned % 10000 == 0 {
                println!("Progress: {:.1}%", (pixels_assigned as f64 / total_pixels as f64) * 100.0);
            }
        }
        
        // Fill any remaining unassigned pixels with nearest plate
        self.fill_unassigned_pixels();
        
        Ok(())
    }

    /// Add neighbors to growth frontier
    fn add_neighbors_to_frontier(&self, x: usize, y: usize, plate_id: u16, 
                                center_point: &Point3D, frontier: &mut BinaryHeap<GrowthPoint>,
                                assigned: &HashSet<usize>) {
        // Get neighbors in pixel space
        let neighbors = self.get_pixel_neighbors(x, y);
        
        for (nx, ny) in neighbors {
            let idx = self.get_pixel_index(nx, ny);
            
            // Skip if already assigned
            if assigned.contains(&idx) {
                continue;
            }
            
            let neighbor_point = self.get_pixel_point(nx, ny);
            let distance = center_point.geodesic_distance(neighbor_point);
            
            frontier.push(GrowthPoint {
                point: *neighbor_point,
                plate_id,
                distance,
            });
        }
    }

    /// Get pixel neighbors with wraparound
    fn get_pixel_neighbors(&self, x: usize, y: usize) -> Vec<(usize, usize)> {
        let mut neighbors = Vec::new();
        
        for dy in -1i32..=1 {
            for dx in -1i32..=1 {
                if dx == 0 && dy == 0 {
                    continue;
                }
                
                let ny = y as i32 + dy;
                
                // Skip out of bounds latitude
                if ny < 0 || ny >= self.height as i32 {
                    continue;
                }
                
                // Handle longitude wraparound
                let mut nx = x as i32 + dx;
                if nx < 0 {
                    nx = self.width as i32 - 1;
                } else if nx >= self.width as i32 {
                    nx = 0;
                }
                
                neighbors.push((nx as usize, ny as usize));
            }
        }
        
        neighbors
    }

    /// Fill any remaining unassigned pixels
    fn fill_unassigned_pixels(&mut self) {
        for y in 0..self.height {
            for x in 0..self.width {
                let idx = self.get_pixel_index(x, y);
                
                if self.plate_map[idx] == 0 {
                    let point = self.get_pixel_point(x, y);
                    
                    // Find nearest plate
                    let mut nearest_plate = 1;
                    let mut min_distance = f64::INFINITY;
                    
                    for seed in &self.plate_seeds {
                        let distance = point.geodesic_distance(&seed.point_3d);
                        if distance < min_distance {
                            min_distance = distance;
                            nearest_plate = seed.id;
                        }
                    }
                    
                    self.plate_map[idx] = nearest_plate;
                }
            }
        }
    }

    /// Smooth boundaries with geodesic-aware smoothing
    pub fn smooth_boundaries(&mut self, iterations: usize) {
        println!("Smoothing boundaries ({} iterations)...", iterations);

        for iter in 0..iterations {
            let original_map = self.plate_map.clone();
            
            for y in 0..self.height {
                for x in 0..self.width {
                    let idx = self.get_pixel_index(x, y);
                    let point = self.get_pixel_point(x, y);
                    
                    // Get neighbors and their weights based on geodesic distance
                    let neighbors = self.get_pixel_neighbors(x, y);
                    let mut weighted_counts: HashMap<u16, f64> = HashMap::new();
                    let mut total_weight = 0.0;
                    
                    for (nx, ny) in neighbors {
                        let neighbor_idx = self.get_pixel_index(nx, ny);
                        let neighbor_point = self.get_pixel_point(nx, ny);
                        let neighbor_plate = original_map[neighbor_idx];
                        
                        // Weight by inverse geodesic distance
                        let distance = point.geodesic_distance(neighbor_point);
                        let weight = 1.0 / (distance + 0.01); // Add small epsilon to avoid division by zero
                        
                        *weighted_counts.entry(neighbor_plate).or_insert(0.0) += weight;
                        total_weight += weight;
                    }
                    
                    // Also consider current pixel with some weight
                    let current_plate = original_map[idx];
                    let self_weight = total_weight * 0.5; // Current pixel has weight equal to half of all neighbors
                    *weighted_counts.entry(current_plate).or_insert(0.0) += self_weight;
                    total_weight += self_weight;
                    
                    // Find plate with highest weighted count
                    if let Some((&most_common, &weight)) = weighted_counts
                        .iter()
                        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
                    {
                        // Only change if there's significant majority
                        if weight / total_weight > 0.6 {
                            self.plate_map[idx] = most_common;
                        }
                    }
                }
            }

            println!("Completed smoothing iteration {}/{}", iter + 1, iterations);
        }
    }

    /// Add realistic boundary noise using fractal patterns
    pub fn add_boundary_noise(&mut self, noise_scale: f64) {
        println!("Adding fractal boundary noise...");
        
        let mut changes = Vec::new();
        
        for y in 1..self.height-1 {
            for x in 0..self.width {
                let idx = self.get_pixel_index(x, y);
                let current_plate = self.plate_map[idx];
                
                // Check if this is a boundary pixel
                let neighbors = self.get_pixel_neighbors(x, y);
                let is_boundary = neighbors.iter().any(|&(nx, ny)| {
                    let n_idx = self.get_pixel_index(nx, ny);
                    self.plate_map[n_idx] != current_plate
                });
                
                if is_boundary {
                    // Use multiple octaves of noise for fractal pattern
                    let point = self.get_pixel_point(x, y);
                    let (lat, lon) = point.to_lat_lon();
                    
                    let mut noise_value = 0.0;
                    let mut amplitude = 1.0;
                    let mut frequency = 1.0;
                    
                    for _ in 0..4 { // 4 octaves
                        let nx = lon * frequency;
                        let ny = lat * frequency;
                        
                        // Simple pseudo-random based on position
                        let hash = ((nx * 12.9898 + ny * 78.233).sin() * 43758.5453).fract();
                        noise_value += hash * amplitude;
                        
                        amplitude *= 0.5;
                        frequency *= 2.0;
                    }
                    
                    // Randomly switch to a neighboring plate based on noise
                    if noise_value.abs() * noise_scale > self.rng.gen::<f64>() {
                        // Find a different neighboring plate
                        let different_plates: Vec<u16> = neighbors.iter()
                            .map(|&(nx, ny)| self.plate_map[self.get_pixel_index(nx, ny)])
                            .filter(|&p| p != current_plate)
                            .collect();
                        
                        if !different_plates.is_empty() {
                            let new_plate = different_plates[self.rng.gen_range(0..different_plates.len())];
                            changes.push((idx, new_plate));
                        }
                    }
                }
            }
        }
        
        // Apply changes
        for (idx, new_plate) in changes {
            self.plate_map[idx] = new_plate;
        }
        
        println!("Added noise to {} boundary pixels", changes.len());
    }

    /// Main generation function
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

        // Step 3: Add boundary noise for realism
        self.add_boundary_noise(0.3);

        // Step 4: Smooth boundaries if requested
        if smooth {
            self.smooth_boundaries(2);
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
                let idx = self.get_pixel_index(x, y);
                let pixel = self.plate_map[idx];
                
                if pixel > 0 {
                    let entry = plate_areas.entry(pixel).or_insert((0, 0.0));
                    entry.0 += 1;  // pixel count
                    entry.1 += pixel_area;  // area
                }
            }
        }

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
    fn test_3d_point_conversion() {
        // Test equator
        let p1 = Point3D::from_lat_lon(0.0, 0.0);
        assert!((p1.x - 1.0).abs() < 0.001);
        assert!(p1.y.abs() < 0.001);
        assert!(p1.z.abs() < 0.001);
        
        // Test north pole
        let p2 = Point3D::from_lat_lon(90.0, 0.0);
        assert!(p2.x.abs() < 0.001);
        assert!(p2.y.abs() < 0.001);
        assert!((p2.z - 1.0).abs() < 0.001);
        
        // Test round trip
        let (lat, lon) = p1.to_lat_lon();
        assert!(lat.abs() < 0.001);
        assert!(lon.abs() < 0.001);
    }

    #[test]
    fn test_geodesic_distance() {
        // Same point
        let p1 = Point3D::from_lat_lon(0.0, 0.0);
        let p2 = Point3D::from_lat_lon(0.0, 0.0);
        assert!(p1.geodesic_distance_km(&p2) < 0.001);
        
        // Opposite points on equator
        let p3 = Point3D::from_lat_lon(0.0, 0.0);
        let p4 = Point3D::from_lat_lon(0.0, 180.0);
        let dist = p3.geodesic_distance_km(&p4);
        let expected = std::f64::consts::PI * PLANET_RADIUS_KM; // Half circumference
        assert!((dist - expected).abs() < 1.0);
        
        // Pole to equator
        let p5 = Point3D::from_lat_lon(90.0, 0.0);
        let p6 = Point3D::from_lat_lon(0.0, 0.0);
        let dist2 = p5.geodesic_distance_km(&p6);
        let expected2 = std::f64::consts::PI * PLANET_RADIUS_KM / 2.0; // Quarter circumference
        assert!((dist2 - expected2).abs() < 1.0);
    }

    #[test]
    fn test_deterministic_generation() {
        let mut gen1 = TectonicPlateGenerator::with_seed(1800, 900, 5, 123).unwrap();
        let mut gen2 = TectonicPlateGenerator::with_seed(1800, 900, 5, 123).unwrap();

        let result1 = gen1.generate("region_growing", false);
        let result2 = gen2.generate("region_growing", false);

        assert!(result1.is_ok() && result2.is_ok());
        // With same seed, results should be identical
        assert_eq!(result1.unwrap(), result2.unwrap());
    }

    #[test]
    fn test_validation() {
        let mut generator = TectonicPlateGenerator::new(1800, 900, 5).unwrap();
        let result = generator.generate("voronoi", false);
        assert!(result.is_ok());
        assert!(generator.validate().is_ok());
    }
}
