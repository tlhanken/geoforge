//! Main tectonic plate generator

use crate::map::terrain::{TerrainMap, PlateMap};
use crate::tectonics::plates::{PlateSeed, PlateStats};
use crate::tectonics::electrostatic::*;
use crate::tectonics::boundary_refinement::{BoundaryRefiner, BoundaryRefinementConfig};
use crate::tectonics::PlateError;
use rand::rngs::StdRng;
use rand::SeedableRng;
use std::collections::HashMap;
use std::fs;
use std::path::Path;


/// Generator for realistic tectonic plates using electrostatic physics simulation
pub struct TectonicPlateGenerator {
    pub width: usize,
    pub height: usize,
    pub num_plates: usize,
    plate_map: PlateMap,
    plate_seeds: Vec<PlateSeed>,
    rng: StdRng,
    current_seed: u64,
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

        let generator = Self {
            width,
            height,
            num_plates,
            plate_map: TerrainMap::new(width, height, 0u16),
            plate_seeds: Vec::with_capacity(num_plates),
            rng: StdRng::seed_from_u64(seed),
            current_seed: seed,
            plate_size_targets: Vec::new(),
        };

        Ok(generator)
    }

    /// Get the current seed
    pub fn get_seed(&self) -> u64 {
        self.current_seed
    }

    /// Reset the generator with a new seed
    pub fn set_seed(&mut self, seed: u64) {
        self.current_seed = seed;
        self.rng = StdRng::seed_from_u64(seed);
        self.plate_map.data.fill(0);
        self.plate_seeds.clear();
        self.plate_size_targets.clear();
    }

    /// Generate plates using electrostatic simulation
    pub fn generate_plates_electrostatic(&mut self) -> Result<(), PlateError> {
        println!("Generating plates using electrostatic simulation...");
        
        // Generate random charges
        let config = ElectrostaticConfig::default();
        let mut charges = generate_random_charges(self.num_plates, &mut self.rng, &config)?;
        
        // Simulate to equilibrium
        simulate_equilibrium(&mut charges, &config)?;
        
        // Convert charges to seeds
        self.plate_seeds = charges_to_seeds(&charges, self.width, self.height, &mut self.rng);
        
        // Generate plates using Voronoi from equilibrium positions
        for y in 0..self.height {
            for x in 0..self.width {
                let point = self.plate_map.get_spherical_point(x, y);

                let mut nearest_plate = 1;
                let mut min_distance = f64::INFINITY;

                // Find closest seed using geodesic distance
                for seed in &self.plate_seeds {
                    let distance = point.distance_to(seed.spherical_point());
                    if distance < min_distance {
                        min_distance = distance;
                        nearest_plate = seed.id;
                    }
                }

                self.plate_map.set(x, y, nearest_plate);
            }

            if y % 100 == 0 {
                println!("Progress: {:.1}%", (y as f64 / self.height as f64) * 100.0);
            }
        }
        
        Ok(())
    }


    /// Apply boundary refinement to add realistic irregularity to plate edges (Stage 1.2)
    ///
    /// This post-processes the Voronoi boundaries from Stage 1.1 to create more natural-looking
    /// plate boundaries with roughness and variation, while keeping the core generation isolated.
    ///
    /// # Arguments
    /// * `config` - Configuration for boundary refinement parameters
    ///
    /// # Example
    /// ```no_run
    /// use geoforge::{TectonicPlateGenerator, BoundaryRefinementConfig};
    ///
    /// let mut generator = TectonicPlateGenerator::with_seed(1800, 900, 15, 42)?;
    /// generator.generate("electrostatic")?;
    ///
    /// // Apply boundary refinement (Stage 1.2)
    /// let config = BoundaryRefinementConfig::with_seed(42)
    ///     .with_noise(0.1, 3.0, 3)
    ///     .with_smoothing(1);
    /// generator.refine_boundaries(config);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn refine_boundaries(&mut self, config: BoundaryRefinementConfig) {
        let mut refiner = BoundaryRefiner::new(config);
        refiner.refine_boundaries(&mut self.plate_map);
    }

    /// Main generation function using electrostatic physics simulation
    pub fn generate(&mut self, _method: &str) -> Result<&Vec<u16>, PlateError> {
        println!("Generating {} tectonic plates using electrostatic physics... (seed: {})",
                 self.num_plates, self.current_seed);

        // Generate plate boundaries using electrostatic simulation
        self.generate_plates_electrostatic()?;

        println!("Tectonic plate generation complete!");
        Ok(&self.plate_map.data)
    }

    /// Generate comprehensive statistics with area-weighted calculations
    pub fn get_plate_stats(&self) -> HashMap<u16, PlateStats> {
        let mut plate_areas = HashMap::new();
        
        // Calculate area for each plate accounting for latitude
        for y in 0..self.height {
            for x in 0..self.width {
                let (lat, _lon) = self.plate_map.projection.pixel_to_coords(x, y);
                let pixel_area = self.plate_map.projection.pixel_area_km2(lat);
                let pixel = self.plate_map.data[self.plate_map.get_index(x, y)];
                
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
                    PlateStats::new(pixels, area_km2, seed.clone()),
                );
            }
        }

        stats
    }

    /// Get plate map dimensions and data
    pub fn get_plate_data(&self) -> (usize, usize, &Vec<u16>, &Vec<PlateSeed>) {
        (self.width, self.height, &self.plate_map.data, &self.plate_seeds)
    }

    /// Export as little-endian bytes for file writing
    pub fn export_raw_u16_le(&self) -> Vec<u8> {
        self.plate_map.data
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
        use rand::prelude::*;
        
        let mut img = ImageBuffer::new(self.width as u32, self.height as u32);
        
        // Generate a color for each plate
        let mut rng = StdRng::seed_from_u64(42); // Consistent colors
        let colors: Vec<[u8; 3]> = (0..=self.num_plates).map(|_| {
            [rng.gen(), rng.gen(), rng.gen()]
        }).collect();
        
        // Draw the image
        for (y, row) in self.plate_map.data.chunks(self.width).enumerate() {
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
        band.write_band_as_slice(0, 0, self.width as usize, self.height as usize, &self.plate_map.data)?;
        
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
        if self.plate_map.data.contains(&0) {
            return Err(PlateError::InvalidParameters("Unassigned pixels found".to_string()));
        }

        // Check all plate IDs are valid
        let max_plate_id = self.num_plates as u16;
        if self.plate_map.data.iter().any(|&p| p > max_plate_id) {
            return Err(PlateError::InvalidParameters("Invalid plate ID found".to_string()));
        }

        // Check all plates have at least some pixels
        let mut plate_counts = vec![0; self.num_plates + 1];
        for &plate_id in &self.plate_map.data {
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
    fn test_deterministic_generation() {
        let mut gen1 = TectonicPlateGenerator::with_seed(1800, 900, 5, 123).unwrap();
        let mut gen2 = TectonicPlateGenerator::with_seed(1800, 900, 5, 123).unwrap();

        let result1 = gen1.generate("electrostatic");
        let result2 = gen2.generate("electrostatic");

        assert!(result1.is_ok() && result2.is_ok());
        // With same seed, results should be identical
        assert_eq!(result1.unwrap(), result2.unwrap());
    }

    #[test]
    fn test_different_seeds_produce_different_maps() {
        let mut gen1 = TectonicPlateGenerator::with_seed(180, 90, 5, 123).unwrap();
        let mut gen2 = TectonicPlateGenerator::with_seed(180, 90, 5, 456).unwrap();

        let result1 = gen1.generate("electrostatic");
        let result2 = gen2.generate("electrostatic");

        assert!(result1.is_ok() && result2.is_ok());
        // With different seeds, results should be different
        assert_ne!(result1.unwrap(), result2.unwrap());
    }

    #[test]
    fn test_validation() {
        let mut generator = TectonicPlateGenerator::new(1800, 900, 5).unwrap();
        let result = generator.generate("electrostatic");
        assert!(result.is_ok());
        assert!(generator.validate().is_ok());
    }
}