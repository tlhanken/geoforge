//! WorldMap - Central data structure for holding and managing all world generation layers
//! 
//! This module provides the main WorldMap class that stores all generated layers
//! (tectonics, elevation, climate, biomes) and handles serialization/export.

use crate::map::TerrainMap;
use crate::tectonics::electrostatic::*;
use crate::tectonics::plates::{PlateSeed, PlateStats};
use rand::rngs::StdRng;
use rand::SeedableRng;
use std::collections::HashMap;
use std::fs;
use std::path::Path;

/// Main world map containing all generated layers
#[derive(Debug, Clone)]
pub struct WorldMap {
    pub width: usize,
    pub height: usize,
    pub seed: u64,
    
    // World generation layers
    pub tectonics: Option<TerrainMap<u16>>,
    pub elevation: Option<TerrainMap<f32>>,
    pub temperature: Option<TerrainMap<f32>>,
    pub precipitation: Option<TerrainMap<f32>>,
    pub biomes: Option<TerrainMap<u8>>,
    
    // Metadata for layers
    pub plate_seeds: Option<Vec<PlateSeed>>,
    pub plate_stats: Option<HashMap<u16, PlateStats>>,
    
    // Internal RNG for consistent generation
    rng: StdRng,
}

impl WorldMap {
    /// Create a new WorldMap with specified dimensions and seed
    pub fn new(width: usize, height: usize, seed: u64) -> Result<Self, Box<dyn std::error::Error>> {
        if width == 0 || height == 0 {
            return Err("Width and height must be > 0".into());
        }

        Ok(Self {
            width,
            height,
            seed,
            tectonics: None,
            elevation: None,
            temperature: None,
            precipitation: None,
            biomes: None,
            plate_seeds: None,
            plate_stats: None,
            rng: StdRng::seed_from_u64(seed),
        })
    }
    
    /// Create a new WorldMap with random seed
    pub fn new_random(width: usize, height: usize) -> Result<Self, Box<dyn std::error::Error>> {
        Self::new(width, height, rand::random())
    }
    
    /// Generate tectonic plates using electrostatic physics simulation
    pub fn generate_tectonics(&mut self, num_plates: usize, smooth: bool) -> Result<(), Box<dyn std::error::Error>> {
        println!("🗺️ Generating {} tectonic plates using electrostatic physics... (seed: {})", 
                 num_plates, self.seed);
        
        if num_plates == 0 || num_plates > u16::MAX as usize {
            return Err("Invalid number of plates".into());
        }
        
        // Generate electrostatic simulation
        let config = ElectrostaticConfig::default();
        let mut charges = generate_random_charges(num_plates, &mut self.rng, &config)?;
        simulate_equilibrium(&mut charges, &config)?;
        
        // Convert charges to seeds
        let plate_seeds = charges_to_seeds(&charges, self.width, self.height, &mut self.rng);
        
        // Generate tectonic plate map
        let mut tectonic_map = TerrainMap::new(self.width, self.height, 0u16);
        
        // Voronoi generation from equilibrium positions
        for y in 0..self.height {
            for x in 0..self.width {
                let point = tectonic_map.get_spherical_point(x, y);
                let mut nearest_plate = 1;
                let mut min_distance = f64::INFINITY;
                
                for seed in &plate_seeds {
                    let distance = point.distance_to(seed.spherical_point());
                    if distance < min_distance {
                        min_distance = distance;
                        nearest_plate = seed.id;
                    }
                }
                
                tectonic_map.set(x, y, nearest_plate);
            }
            
            if y % 100 == 0 {
                println!("Progress: {:.1}%", (y as f64 / self.height as f64) * 100.0);
            }
        }
        
        // Optional smoothing
        if smooth {
            self.smooth_tectonic_boundaries(&mut tectonic_map, 2);
        }
        
        // Calculate statistics
        let plate_stats = self.calculate_plate_stats(&tectonic_map, &plate_seeds);
        
        // Store results
        self.tectonics = Some(tectonic_map);
        self.plate_seeds = Some(plate_seeds);
        self.plate_stats = Some(plate_stats);
        
        println!("✅ Tectonic plate generation complete!");
        Ok(())
    }
    
    /// Smooth tectonic boundaries with geodesic-aware smoothing
    fn smooth_tectonic_boundaries(&self, tectonic_map: &mut TerrainMap<u16>, iterations: usize) {
        println!("Smoothing boundaries ({} iterations)...", iterations);

        for iter in 0..iterations {
            let original_data = tectonic_map.data.clone();
            
            for y in 0..self.height {
                for x in 0..self.width {
                    let idx = tectonic_map.get_index(x, y);
                    let point = tectonic_map.get_spherical_point(x, y);
                    
                    let neighbors = tectonic_map.get_neighbors(x, y);
                    let mut weighted_counts: HashMap<u16, f64> = HashMap::new();
                    let mut total_weight = 0.0;
                    
                    for (nx, ny) in neighbors {
                        let neighbor_idx = tectonic_map.get_index(nx, ny);
                        let neighbor_point = tectonic_map.get_spherical_point(nx, ny);
                        let neighbor_plate = original_data[neighbor_idx];
                        
                        let distance = point.distance_to(&neighbor_point);
                        let weight = 1.0 / (distance + 0.01);
                        
                        *weighted_counts.entry(neighbor_plate).or_insert(0.0) += weight;
                        total_weight += weight;
                    }
                    
                    let current_plate = original_data[idx];
                    let self_weight = total_weight * 0.5;
                    *weighted_counts.entry(current_plate).or_insert(0.0) += self_weight;
                    total_weight += self_weight;
                    
                    if let Some((&most_common, &weight)) = weighted_counts
                        .iter()
                        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
                    {
                        if weight / total_weight > 0.6 {
                            tectonic_map.data[idx] = most_common;
                        }
                    }
                }
            }

            println!("Completed smoothing iteration {}/{}", iter + 1, iterations);
        }
    }
    
    /// Calculate comprehensive plate statistics
    fn calculate_plate_stats(&self, tectonic_map: &TerrainMap<u16>, plate_seeds: &[PlateSeed]) -> HashMap<u16, PlateStats> {
        let mut plate_areas = HashMap::new();
        
        for y in 0..self.height {
            for x in 0..self.width {
                let (lat, _lon) = tectonic_map.projection.pixel_to_coords(x, y);
                let pixel_area = tectonic_map.projection.pixel_area_km2(lat);
                let pixel = tectonic_map.data[tectonic_map.get_index(x, y)];
                
                if pixel > 0 {
                    let entry = plate_areas.entry(pixel).or_insert((0, 0.0));
                    entry.0 += 1;
                    entry.1 += pixel_area;
                }
            }
        }
        
        let mut stats = HashMap::new();
        for (&plate_id, &(pixels, area_km2)) in &plate_areas {
            if let Some(seed) = plate_seeds.iter().find(|s| s.id == plate_id) {
                stats.insert(plate_id, PlateStats::new(pixels, area_km2, seed.clone()));
            }
        }
        
        stats
    }
    
    /// Get statistics for the tectonic layer
    pub fn get_tectonic_stats(&self) -> Option<&HashMap<u16, PlateStats>> {
        self.plate_stats.as_ref()
    }
    
    /// Export a single layer as PNG
    #[cfg(feature = "export-png")]
    pub fn export_layer_png(&self, layer: LayerType, output_dir: &str, filename: &str) -> Result<(), Box<dyn std::error::Error>> {
        use image::{ImageBuffer, Rgb, Luma};
        use rand::prelude::*;
        
        fs::create_dir_all(output_dir)?;
        let path = Path::new(output_dir).join(filename);
        
        match layer {
            LayerType::Tectonics => {
                if let Some(ref tectonic_map) = self.tectonics {
                    let mut img = ImageBuffer::new(self.width as u32, self.height as u32);
                    
                    // Generate consistent colors for each plate
                    let mut rng = StdRng::seed_from_u64(42);
                    let max_plate_id = tectonic_map.data.iter().max().unwrap_or(&0);
                    let colors: Vec<[u8; 3]> = (0..=*max_plate_id).map(|_| {
                        [rng.gen(), rng.gen(), rng.gen()]
                    }).collect();
                    
                    for (y, row) in tectonic_map.data.chunks(self.width).enumerate() {
                        for (x, &plate_id) in row.iter().enumerate() {
                            let color = colors[plate_id as usize];
                            img.put_pixel(x as u32, y as u32, Rgb(color));
                        }
                    }
                    
                    img.save(path)?;
                } else {
                    return Err("Tectonic layer not generated".into());
                }
            }
            LayerType::Elevation => {
                if let Some(ref elevation_map) = self.elevation {
                    let mut img = ImageBuffer::new(self.width as u32, self.height as u32);
                    
                    let min_val = elevation_map.data.iter().fold(f32::INFINITY, |a, &b| a.min(b));
                    let max_val = elevation_map.data.iter().fold(f32::NEG_INFINITY, |a, &b| a.max(b));
                    let range = max_val - min_val;
                    
                    for (y, row) in elevation_map.data.chunks(self.width).enumerate() {
                        for (x, &value) in row.iter().enumerate() {
                            let normalized = ((value - min_val) / range * 255.0) as u8;
                            img.put_pixel(x as u32, y as u32, Luma([normalized]));
                        }
                    }
                    
                    img.save(path)?;
                } else {
                    return Err("Elevation layer not generated".into());
                }
            }
            // Add other layer types as needed
            _ => return Err("Layer type not yet implemented".into()),
        }
        
        Ok(())
    }
    
    /// Serialize the entire WorldMap to a binary .map file
    pub fn save_to_file(&self, filepath: &str) -> Result<(), Box<dyn std::error::Error>> {
        use std::io::Write;
        
        let mut file = fs::File::create(filepath)?;
        
        // Write header
        file.write_all(b"GEOFORGE_MAP_V1\0")?;
        file.write_all(&(self.width as u32).to_le_bytes())?;
        file.write_all(&(self.height as u32).to_le_bytes())?;
        file.write_all(&self.seed.to_le_bytes())?;
        
        // Write layer flags
        let flags = self.get_layer_flags();
        file.write_all(&flags.to_le_bytes())?;
        
        // Write layer data
        if let Some(ref tectonic_map) = self.tectonics {
            for &value in &tectonic_map.data {
                file.write_all(&value.to_le_bytes())?;
            }
        }
        
        if let Some(ref elevation_map) = self.elevation {
            for &value in &elevation_map.data {
                file.write_all(&value.to_le_bytes())?;
            }
        }
        
        // Write plate seeds if available
        if let Some(ref seeds) = self.plate_seeds {
            file.write_all(&(seeds.len() as u32).to_le_bytes())?;
            for seed in seeds {
                file.write_all(&seed.id.to_le_bytes())?;
                file.write_all(&(seed.x as u32).to_le_bytes())?;
                file.write_all(&(seed.y as u32).to_le_bytes())?;
                file.write_all(&seed.lat.to_le_bytes())?;
                file.write_all(&seed.lon.to_le_bytes())?;
                file.write_all(&seed.motion_direction.to_le_bytes())?;
                file.write_all(&seed.motion_speed.to_le_bytes())?;
            }
        }
        
        println!("✅ WorldMap saved to {}", filepath);
        Ok(())
    }
    
    /// Load a WorldMap from a binary .map file
    pub fn load_from_file(filepath: &str) -> Result<Self, Box<dyn std::error::Error>> {
        use std::io::Read;
        
        let mut file = fs::File::open(filepath)?;
        let mut buffer = Vec::new();
        file.read_to_end(&mut buffer)?;
        
        let mut cursor = 0;
        
        // Read header
        let header = &buffer[cursor..cursor+16];
        if header != b"GEOFORGE_MAP_V1\0" {
            return Err("Invalid file format".into());
        }
        cursor += 16;
        
        let width = u32::from_le_bytes([buffer[cursor], buffer[cursor+1], buffer[cursor+2], buffer[cursor+3]]) as usize;
        cursor += 4;
        let height = u32::from_le_bytes([buffer[cursor], buffer[cursor+1], buffer[cursor+2], buffer[cursor+3]]) as usize;
        cursor += 4;
        let seed = u64::from_le_bytes([
            buffer[cursor], buffer[cursor+1], buffer[cursor+2], buffer[cursor+3],
            buffer[cursor+4], buffer[cursor+5], buffer[cursor+6], buffer[cursor+7]
        ]);
        cursor += 8;
        
        let flags = u32::from_le_bytes([buffer[cursor], buffer[cursor+1], buffer[cursor+2], buffer[cursor+3]]);
        cursor += 4;
        
        // Create the worldmap with proper validation
        if width == 0 || height == 0 {
            return Err("Invalid dimensions in file".into());
        }
        
        let mut worldmap = Self {
            width,
            height,
            seed,
            tectonics: None,
            elevation: None,
            temperature: None,
            precipitation: None,
            biomes: None,
            plate_seeds: None,
            plate_stats: None,
            rng: StdRng::seed_from_u64(seed),
        };
        
        // Read tectonic layer
        if flags & 0x01 != 0 {
            let mut tectonic_data = vec![0u16; width * height];
            for i in 0..width * height {
                tectonic_data[i] = u16::from_le_bytes([buffer[cursor], buffer[cursor+1]]);
                cursor += 2;
            }
            worldmap.tectonics = Some(TerrainMap::from_data(width, height, tectonic_data));
        }
        
        // Read elevation layer
        if flags & 0x02 != 0 {
            let mut elevation_data = vec![0f32; width * height];
            for i in 0..width * height {
                elevation_data[i] = f32::from_le_bytes([
                    buffer[cursor], buffer[cursor+1], buffer[cursor+2], buffer[cursor+3]
                ]);
                cursor += 4;
            }
            worldmap.elevation = Some(TerrainMap::from_data(width, height, elevation_data));
        }
        
        println!("✅ WorldMap loaded from {}", filepath);
        Ok(worldmap)
    }
    
    /// Get flags indicating which layers are present
    fn get_layer_flags(&self) -> u32 {
        let mut flags = 0u32;
        if self.tectonics.is_some() { flags |= 0x01; }
        if self.elevation.is_some() { flags |= 0x02; }
        if self.temperature.is_some() { flags |= 0x04; }
        if self.precipitation.is_some() { flags |= 0x08; }
        if self.biomes.is_some() { flags |= 0x10; }
        flags
    }
}

/// Available layers for export
#[derive(Debug, Clone, Copy)]
pub enum LayerType {
    Tectonics,
    Elevation,
    Temperature,
    Precipitation,
    Biomes,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_worldmap_creation() {
        let worldmap = WorldMap::new(1800, 900, 42).unwrap();
        assert_eq!(worldmap.width, 1800);
        assert_eq!(worldmap.height, 900);
        assert_eq!(worldmap.seed, 42);
    }

    #[test]
    fn test_tectonic_generation() {
        let mut worldmap = WorldMap::new(100, 50, 123).unwrap();
        let result = worldmap.generate_tectonics(5, true);
        assert!(result.is_ok());
        assert!(worldmap.tectonics.is_some());
        assert!(worldmap.plate_seeds.is_some());
        assert!(worldmap.plate_stats.is_some());
    }
    
    #[test]
    fn test_serialization_roundtrip() {
        let mut worldmap = WorldMap::new(100, 50, 42).unwrap();
        worldmap.generate_tectonics(5, false).unwrap();
        
        // Save to file
        let temp_file = "test_worldmap.map";
        assert!(worldmap.save_to_file(temp_file).is_ok());
        
        // Load from file
        let loaded_worldmap = WorldMap::load_from_file(temp_file).unwrap();
        
        // Verify dimensions and seed
        assert_eq!(loaded_worldmap.width, 100);
        assert_eq!(loaded_worldmap.height, 50);
        assert_eq!(loaded_worldmap.seed, 42);
        
        // Verify tectonic data is preserved
        assert!(loaded_worldmap.tectonics.is_some());
        let original_data = worldmap.tectonics.as_ref().unwrap();
        let loaded_data = loaded_worldmap.tectonics.as_ref().unwrap();
        assert_eq!(original_data.data, loaded_data.data);
        
        // Cleanup
        std::fs::remove_file(temp_file).ok();
    }
    
    #[test]
    fn test_serialization_invalid_file() {
        let result = WorldMap::load_from_file("nonexistent_file.map");
        assert!(result.is_err());
    }
    
    #[test]
    fn test_serialization_empty_worldmap() {
        let worldmap = WorldMap::new(50, 25, 999).unwrap();
        let temp_file = "test_empty_worldmap.map";
        
        // Save empty worldmap
        assert!(worldmap.save_to_file(temp_file).is_ok());
        
        // Load and verify
        let loaded = WorldMap::load_from_file(temp_file).unwrap();
        assert_eq!(loaded.width, 50);
        assert_eq!(loaded.height, 25);
        assert_eq!(loaded.seed, 999);
        assert!(loaded.tectonics.is_none());
        
        // Cleanup
        std::fs::remove_file(temp_file).ok();
    }
    
    #[test]
    fn test_seed_reproducibility() {
        // Test that same seed produces identical results
        let mut worldmap1 = WorldMap::new(50, 25, 12345).unwrap();
        let mut worldmap2 = WorldMap::new(50, 25, 12345).unwrap();
        
        worldmap1.generate_tectonics(5, false).unwrap();
        worldmap2.generate_tectonics(5, false).unwrap();
        
        // Compare tectonic data
        let data1 = &worldmap1.tectonics.as_ref().unwrap().data;
        let data2 = &worldmap2.tectonics.as_ref().unwrap().data;
        assert_eq!(data1, data2);
        
        // Compare plate seeds
        let seeds1 = worldmap1.plate_seeds.as_ref().unwrap();
        let seeds2 = worldmap2.plate_seeds.as_ref().unwrap();
        assert_eq!(seeds1.len(), seeds2.len());
        
        for (s1, s2) in seeds1.iter().zip(seeds2.iter()) {
            assert_eq!(s1.id, s2.id);
            assert!((s1.lat - s2.lat).abs() < 1e-10);
            assert!((s1.lon - s2.lon).abs() < 1e-10);
        }
    }
    
    #[test]
    fn test_failure_modes() {
        // Test invalid dimensions
        assert!(WorldMap::new(0, 100, 42).is_err());
        assert!(WorldMap::new(100, 0, 42).is_err());
        
        // Test invalid plate count
        let mut worldmap = WorldMap::new(100, 50, 42).unwrap();
        assert!(worldmap.generate_tectonics(0, false).is_err());
        assert!(worldmap.generate_tectonics(u16::MAX as usize + 1, false).is_err());
    }
    
    #[test]
    fn test_plate_size_ratios() {
        let mut worldmap = WorldMap::new(200, 100, 42).unwrap();
        worldmap.generate_tectonics(20, false).unwrap();
        
        let stats = worldmap.get_tectonic_stats().unwrap();
        let mut areas: Vec<u64> = stats.values().map(|s| s.area_km2).collect();
        areas.sort_by(|a, b| b.cmp(a));
        
        // Check that we have significant size variety
        assert!(areas.len() > 0);
        let largest = areas[0];
        let smallest = areas[areas.len() - 1];
        
        // Should have at least 10x size ratio
        assert!(largest as f64 / smallest as f64 > 10.0);
        
        // All areas should be positive
        for area in areas {
            assert!(area > 0);
        }
    }
}