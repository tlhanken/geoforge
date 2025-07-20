//! WorldMap - Central data structure for holding and managing all world generation layers
//! 
//! This module provides the main WorldMap class that stores all generated layers
//! (tectonics, elevation, climate, biomes) and handles serialization/export.

use crate::map::terrain::TerrainMap;
use crate::map::spherical::PlanetaryParams;
use crate::tectonics::TectonicPlateGenerator;
use crate::tectonics::plates::{PlateSeed, PlateStats};
use std::collections::HashMap;
use std::fs;

/// Main world map containing all generated layers
#[derive(Debug, Clone)]
pub struct WorldMap {
    pub width: usize,
    pub height: usize,
    pub seed: u64,
    pub planetary_params: PlanetaryParams,
    
    // World generation layers
    pub tectonics: Option<TerrainMap<u16>>,
    pub elevation: Option<TerrainMap<f32>>,
    pub temperature: Option<TerrainMap<f32>>,
    pub precipitation: Option<TerrainMap<f32>>,
    pub biomes: Option<TerrainMap<u8>>,
    
    // Metadata for layers
    pub plate_seeds: Option<Vec<PlateSeed>>,
    pub plate_stats: Option<HashMap<u16, PlateStats>>,
}

impl WorldMap {
    /// Create a new WorldMap with specified dimensions and seed (using Earth-like parameters)
    pub fn new(width: usize, height: usize, seed: u64) -> Result<Self, Box<dyn std::error::Error>> {
        Self::new_with_params(width, height, seed, PlanetaryParams::earth())
    }
    
    /// Create a new WorldMap with specified dimensions, seed, and planetary parameters
    pub fn new_with_params(
        width: usize, 
        height: usize, 
        seed: u64, 
        planetary_params: PlanetaryParams
    ) -> Result<Self, Box<dyn std::error::Error>> {
        if width == 0 || height == 0 {
            return Err("Width and height must be > 0".into());
        }

        Ok(Self {
            width,
            height,
            seed,
            planetary_params,
            tectonics: None,
            elevation: None,
            temperature: None,
            precipitation: None,
            biomes: None,
            plate_seeds: None,
            plate_stats: None,
        })
    }
    
    /// Create a Mars-like world
    pub fn new_mars(width: usize, height: usize, seed: u64) -> Result<Self, Box<dyn std::error::Error>> {
        Self::new_with_params(width, height, seed, PlanetaryParams::mars())
    }
    
    /// Create a Venus-like world
    pub fn new_venus(width: usize, height: usize, seed: u64) -> Result<Self, Box<dyn std::error::Error>> {
        Self::new_with_params(width, height, seed, PlanetaryParams::venus())
    }
    
    /// Create a new WorldMap with random seed
    pub fn new_random(width: usize, height: usize) -> Result<Self, Box<dyn std::error::Error>> {
        Self::new(width, height, rand::random())
    }
    
    /// Import tectonic plates from a PNG image where each color represents a different plate
    #[cfg(feature = "export-png")]
    pub fn import_tectonics_png(&mut self, png_path: &str) -> Result<(), Box<dyn std::error::Error>> {
        use image::io::Reader as ImageReader;
        use std::collections::HashMap;
        
        println!("🖼️ Importing tectonic plates from PNG: {}", png_path);
        
        // Load the PNG image
        let img = ImageReader::open(png_path)?.decode()?;
        let rgb_img = img.to_rgb8();
        let (img_width, img_height) = rgb_img.dimensions();
        
        // Ensure image dimensions match world map
        if img_width as usize != self.width || img_height as usize != self.height {
            return Err(format!(
                "Image dimensions ({}×{}) don't match world map ({}×{})",
                img_width, img_height, self.width, self.height
            ).into());
        }
        
        // Map colors to plate IDs
        let mut color_to_plate_id: HashMap<[u8; 3], u16> = HashMap::new();
        let mut next_plate_id = 1u16;
        let mut tectonic_data = vec![0u16; self.width * self.height];
        
        // Process each pixel
        for y in 0..self.height {
            for x in 0..self.width {
                let pixel = rgb_img.get_pixel(x as u32, y as u32);
                let color = [pixel[0], pixel[1], pixel[2]];
                
                // Get or assign plate ID for this color
                let plate_id = if let Some(&existing_id) = color_to_plate_id.get(&color) {
                    existing_id
                } else {
                    let new_id = next_plate_id;
                    color_to_plate_id.insert(color, new_id);
                    next_plate_id += 1;
                    
                    if next_plate_id > u16::MAX {
                        return Err("Too many unique colors in PNG (max 65535 plates)".into());
                    }
                    
                    new_id
                };
                
                let idx = y * self.width + x;
                tectonic_data[idx] = plate_id;
            }
        }
        
        let num_plates = color_to_plate_id.len();
        println!("📊 Detected {} unique colors/plates", num_plates);
        
        // Generate plate seeds from the imported data
        let plate_seeds = self.generate_seeds_from_plate_data(&tectonic_data)?;
        
        // Calculate statistics
        let tectonic_map = TerrainMap::from_data(self.width, self.height, tectonic_data);
        let plate_stats = self.calculate_plate_stats_from_map(&tectonic_map, &plate_seeds);
        
        // Store results
        self.tectonics = Some(tectonic_map);
        self.plate_seeds = Some(plate_seeds);
        self.plate_stats = Some(plate_stats);
        
        println!("✅ Tectonic plates imported successfully from PNG");
        Ok(())
    }

    /// Generate tectonic plates using electrostatic physics simulation
    pub fn generate_tectonics(&mut self, num_plates: usize, smooth: bool) -> Result<(), Box<dyn std::error::Error>> {
        // Use the dedicated TectonicPlateGenerator
        let mut generator = TectonicPlateGenerator::with_seed(
            self.width, 
            self.height, 
            num_plates, 
            self.seed
        ).map_err(|e| -> Box<dyn std::error::Error> { Box::new(e) })?;
        
        // Generate the plates
        generator.generate("electrostatic", smooth)
            .map_err(|e| -> Box<dyn std::error::Error> { Box::new(e) })?;
        
        // Extract data from generator
        let (width, height, plate_data, seeds) = generator.get_plate_data();
        let plate_stats = generator.get_plate_stats();
        
        // Store results in WorldMap
        self.tectonics = Some(TerrainMap::from_data(width, height, plate_data.clone()));
        self.plate_seeds = Some(seeds.clone());
        self.plate_stats = Some(plate_stats);
        
        Ok(())
    }
    
    
    /// Generate plate seeds from imported plate data by finding centroids
    fn generate_seeds_from_plate_data(&self, plate_data: &[u16]) -> Result<Vec<PlateSeed>, Box<dyn std::error::Error>> {
        use std::collections::HashMap;
        
        // Group pixels by plate ID to find centroids
        let mut plate_pixels: HashMap<u16, Vec<(usize, usize)>> = HashMap::new();
        
        for y in 0..self.height {
            for x in 0..self.width {
                let idx = y * self.width + x;
                let plate_id = plate_data[idx];
                
                if plate_id > 0 {
                    plate_pixels.entry(plate_id).or_insert_with(Vec::new).push((x, y));
                }
            }
        }
        
        let mut seeds = Vec::new();
        
        for (plate_id, pixels) in plate_pixels {
            if pixels.is_empty() {
                continue;
            }
            
            // Calculate centroid
            let sum_x: usize = pixels.iter().map(|(x, _)| x).sum();
            let sum_y: usize = pixels.iter().map(|(_, y)| y).sum();
            let count = pixels.len();
            
            let centroid_x = sum_x / count;
            let centroid_y = sum_y / count;
            
            // Convert to geographic coordinates
            let projection = crate::map::terrain::MapProjection::global_equirectangular(self.width, self.height);
            let (lat, lon) = projection.pixel_to_coords(centroid_x, centroid_y);
            
            // Create seed with default motion parameters
            let seed = PlateSeed::new(
                plate_id,
                centroid_x,
                centroid_y,
                lat,
                lon,
                0.0,  // Default motion direction
                2.5,  // Default motion speed (cm/year)
            );
            
            seeds.push(seed);
        }
        
        seeds.sort_by_key(|s| s.id);
        Ok(seeds)
    }
    
    /// Calculate plate statistics from a terrain map and seeds
    fn calculate_plate_stats_from_map(&self, tectonic_map: &TerrainMap<u16>, plate_seeds: &[PlateSeed]) -> HashMap<u16, PlateStats> {
        let mut plate_areas = HashMap::new();
        
        for y in 0..self.height {
            for x in 0..self.width {
                let (lat, _lon) = tectonic_map.projection.pixel_to_coords(x, y);
                let pixel_area = tectonic_map.projection.pixel_area_km2_with_radius(lat, self.planetary_params.radius_km);
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
    
    /// Export tectonic plates as PNG
    #[cfg(feature = "export-png")]
    pub fn export_tectonics_png(&self, output_dir: &str, filename: &str) -> Result<(), Box<dyn std::error::Error>> {
        use image::{ImageBuffer, Rgb};
        use rand::prelude::*;
        
        fs::create_dir_all(output_dir)?;
        let path = std::path::Path::new(output_dir).join(filename);
        
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
        
        Ok(())
    }

    /// Export all available layers as PNG files
    #[cfg(feature = "export-png")]
    pub fn export_all_png(&self, output_dir: &str, base_name: &str) -> Result<(), Box<dyn std::error::Error>> {
        fs::create_dir_all(output_dir)?;
        
        // Export tectonics if available
        if self.tectonics.is_some() {
            self.export_tectonics_png(output_dir, &format!("{}_tectonics.png", base_name))?;
            println!("✅ Exported tectonics: {}/{}_tectonics.png", output_dir, base_name);
        }
        
        // Export elevation if available (placeholder for future implementation)
        if self.elevation.is_some() {
            // TODO: Implement export_elevation_png
            println!("⚠️ Elevation export not yet implemented");
        }
        
        // Export temperature if available (placeholder for future implementation)
        if self.temperature.is_some() {
            // TODO: Implement export_temperature_png
            println!("⚠️ Temperature export not yet implemented");
        }
        
        // Export precipitation if available (placeholder for future implementation)
        if self.precipitation.is_some() {
            // TODO: Implement export_precipitation_png
            println!("⚠️ Precipitation export not yet implemented");
        }
        
        // Export biomes if available (placeholder for future implementation)
        if self.biomes.is_some() {
            // TODO: Implement export_biomes_png
            println!("⚠️ Biomes export not yet implemented");
        }
        
        if self.tectonics.is_none() && self.elevation.is_none() && 
           self.temperature.is_none() && self.precipitation.is_none() && 
           self.biomes.is_none() {
            println!("⚠️ No layers available to export");
        }
        
        Ok(())
    }

    /// Generate all available layers in sequence
    pub fn generate_all(&mut self, num_plates: usize, smooth: bool) -> Result<(), Box<dyn std::error::Error>> {
        println!("🌍 Generating all world layers...");
        
        // Generate tectonics (foundation layer)
        println!("⚡ Generating tectonic layer...");
        self.generate_tectonics(num_plates, smooth)?;
        
        // Generate elevation from tectonics (placeholder for future implementation)
        println!("🏔️ Elevation generation not yet implemented");
        // TODO: self.generate_elevation()?;
        
        // Generate climate from elevation (placeholder for future implementation)
        println!("🌡️ Temperature generation not yet implemented");
        // TODO: self.generate_temperature()?;
        
        println!("🌧️ Precipitation generation not yet implemented");
        // TODO: self.generate_precipitation()?;
        
        // Generate biomes from climate (placeholder for future implementation)
        println!("🌿 Biomes generation not yet implemented");
        // TODO: self.generate_biomes()?;
        
        println!("✅ All available layers generated!");
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
            planetary_params: PlanetaryParams::earth(), // Default to Earth parameters for loaded files
            tectonics: None,
            elevation: None,
            temperature: None,
            precipitation: None,
            biomes: None,
            plate_seeds: None,
            plate_stats: None,
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
    
    #[test]
    #[cfg(feature = "export-png")]
    fn test_png_export_import_roundtrip() {
        use std::fs;
        
        // Create test directory
        let test_dir = "outputs/examples/test_png_roundtrip";
        fs::create_dir_all(test_dir).unwrap();
        
        // Step 1: Generate world with seed 42 and export PNG
        let mut original_world = WorldMap::new(120, 60, 42).unwrap();
        original_world.generate_tectonics(6, false).unwrap();
        
        let png_path = format!("{}/test_plates.png", test_dir);
        original_world.export_tectonics_png(test_dir, "test_plates.png").unwrap();
        
        // Step 2: Import PNG into new world object
        let mut imported_world = WorldMap::new(120, 60, 999).unwrap(); // Different seed
        imported_world.import_tectonics_png(&png_path).unwrap();
        
        // Step 3: Compare the worlds
        let original_data = &original_world.tectonics.as_ref().unwrap().data;
        let imported_data = &imported_world.tectonics.as_ref().unwrap().data;
        
        // The plate boundaries should be preserved, but plate IDs might be remapped
        // Let's check that the number of plates and their relative sizes are preserved
        let original_stats = original_world.get_tectonic_stats().unwrap();
        let imported_stats = imported_world.get_tectonic_stats().unwrap();
        
        // Should have same number of plates
        assert_eq!(original_stats.len(), imported_stats.len());
        
        // Should have same total area coverage
        let original_total_area: u64 = original_stats.values().map(|s| s.area_km2).sum();
        let imported_total_area: u64 = imported_stats.values().map(|s| s.area_km2).sum();
        assert_eq!(original_total_area, imported_total_area);
        
        // Should have same area distribution (when sorted by size)
        let mut original_areas: Vec<u64> = original_stats.values().map(|s| s.area_km2).collect();
        let mut imported_areas: Vec<u64> = imported_stats.values().map(|s| s.area_km2).collect();
        original_areas.sort_by(|a, b| b.cmp(a));
        imported_areas.sort_by(|a, b| b.cmp(a));
        assert_eq!(original_areas, imported_areas);
        
        // Step 4: Verify boundary preservation by checking if plate groupings are preserved
        // Create a mapping from original plate IDs to imported plate IDs based on overlapping pixels
        use std::collections::HashMap;
        let mut plate_mapping: HashMap<u16, u16> = HashMap::new();
        
        for (orig_id, orig_stat) in original_stats {
            // Find the imported plate that has the same area as this original plate
            for (imp_id, imp_stat) in imported_stats {
                if orig_stat.area_km2 == imp_stat.area_km2 && !plate_mapping.values().any(|&v| v == *imp_id) {
                    plate_mapping.insert(*orig_id, *imp_id);
                    break;
                }
            }
        }
        
        // Verify that every original plate was mapped to an imported plate
        assert_eq!(plate_mapping.len(), original_stats.len());
        
        // Step 5: Verify exact boundary preservation using the mapping
        let mut boundaries_preserved = true;
        for (_i, (&orig_plate, &imp_plate)) in original_data.iter().zip(imported_data.iter()).enumerate() {
            let expected_imp_plate = plate_mapping[&orig_plate];
            if imp_plate != expected_imp_plate {
                boundaries_preserved = false;
                break;
            }
        }
        
        assert!(boundaries_preserved, "Plate boundaries were not exactly preserved through PNG roundtrip");
        
        // Check that every pixel has a valid plate assignment
        for &plate_id in imported_data {
            assert!(plate_id > 0);
            assert!(imported_stats.contains_key(&plate_id));
        }
        
        // Cleanup
        fs::remove_dir_all(test_dir).ok();
    }
    
    #[test]
    fn test_generate_all() {
        let mut world = WorldMap::new(60, 30, 123).unwrap();
        
        // Test generate_all function
        let result = world.generate_all(4, false);
        assert!(result.is_ok());
        
        // Should have generated tectonics
        assert!(world.tectonics.is_some());
        assert!(world.plate_seeds.is_some());
        assert!(world.plate_stats.is_some());
        
        // Other layers should still be None (not yet implemented)
        assert!(world.elevation.is_none());
        assert!(world.temperature.is_none());
        assert!(world.precipitation.is_none());
        assert!(world.biomes.is_none());
        
        // Verify tectonic data is valid
        if let Some(stats) = world.get_tectonic_stats() {
            assert_eq!(stats.len(), 4);
            let total_area: u64 = stats.values().map(|s| s.area_km2).sum();
            assert!(total_area > 0);
        }
    }
}