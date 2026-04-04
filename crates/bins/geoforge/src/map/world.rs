//! WorldMap - Central data structure for holding and managing all world generation layers
//! 
//! This module provides the main WorldMap class that stores all generated layers
//! (tectonics, elevation, climate, biomes) and handles serialization/export.

use crate::map::terrain::TerrainMap;
use crate::map::spherical::PlanetaryParams;
use crate::tectonics::{TectonicPlateGenerator, GenerationMethod};
use crate::tectonics::plates::{PlateSeed, PlateStats, PlateType};
use crate::tectonics::boundary_refinement::{BoundaryRefiner, BoundaryRefinementConfig};
use crate::tectonics::island_removal::{IslandRemover, IslandRemovalConfig, IslandRemovalStats};
use crate::tectonics::boundary_analysis::{BoundaryAnalyzer, BoundaryAnalysisConfig, BoundarySegment, BoundaryStatistics};
use crate::tectonics::motion::{PlateMotionAssigner, PlateMotionConfig};
use crate::geology::provinces::ProvinceRegion;
use std::collections::HashMap;
use std::fs;
use serde::{Serialize, Deserialize};
use bincode;

/// Metadata for Stage 1: Tectonic Foundation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TectonicMetadata {
    /// Plate seed points with motion information
    pub plate_seeds: Vec<PlateSeed>,

    /// Statistics for each plate (area, percentage, type)
    pub plate_stats: HashMap<u16, PlateStats>,

    /// Identified boundary segments between plates
    pub plate_boundaries: Vec<BoundarySegment>,

    /// Statistics about all boundaries
    pub boundary_stats: Option<BoundaryStatistics>,
}

impl TectonicMetadata {
    /// Create new tectonic metadata with seeds and stats
    pub fn new(plate_seeds: Vec<PlateSeed>, plate_stats: HashMap<u16, PlateStats>) -> Self {
        Self {
            plate_seeds,
            plate_stats,
            plate_boundaries: Vec::new(),
            boundary_stats: None,
        }
    }

    /// Assign plate types based on size distribution
    pub fn assign_plate_types(&mut self) {
        // Calculate size percentiles first
        let total_plates = self.plate_stats.len();
        let mut percentiles = HashMap::new();

        for (plate_id, stats) in &self.plate_stats {
            let smaller_count = self.plate_stats.values()
                .filter(|s| s.area_km2 < stats.area_km2)
                .count();
            let percentile = if total_plates > 0 {
                smaller_count as f64 / total_plates as f64
            } else {
                0.5
            };
            percentiles.insert(*plate_id, percentile);
        }

        // Now assign types based on percentiles
        for (plate_id, stats) in self.plate_stats.iter_mut() {
            if let Some(&percentile) = percentiles.get(plate_id) {
                stats.plate_type = PlateType::from_size_percentile(percentile);
            }
        }
    }
    /// Write complete TectonicMetadata to file
    pub fn write_to_file(file: &mut fs::File, metadata: &TectonicMetadata) -> Result<(), Box<dyn std::error::Error>> {
        bincode::serialize_into(file, metadata)?;
        Ok(())
    }

    /// Read complete TectonicMetadata from buffer
    pub fn read_from_buffer(buffer: &[u8], cursor: &mut usize) -> Result<TectonicMetadata, Box<dyn std::error::Error>> {
        let mut reader = std::io::Cursor::new(&buffer[*cursor..]);
        let metadata: TectonicMetadata = bincode::deserialize_from(&mut reader)?;
        *cursor += reader.position() as usize;
        Ok(metadata)
    }
}

/// Main world map containing all generated layers
#[derive(Debug, Clone)]
pub struct WorldMap {
    pub width: usize,
    pub height: usize,
    pub seed: u64,
    pub planetary_params: PlanetaryParams,

    // World generation layers
    pub tectonics: Option<TerrainMap<u16>>,
    pub geology: Option<Vec<ProvinceRegion>>,  // Stage 2: Geological provinces
    pub elevation: Option<TerrainMap<f32>>,
    pub temperature: Option<TerrainMap<f32>>,
    pub precipitation: Option<TerrainMap<f32>>,
    pub biomes: Option<TerrainMap<u8>>,

    // Metadata for layers
    pub tectonic_metadata: Option<TectonicMetadata>,
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
            geology: None,
            elevation: None,
            temperature: None,
            precipitation: None,
            biomes: None,
            tectonic_metadata: None,
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

    /// Access tectonic operations module
    ///
    /// Provides organized access to all tectonic plate operations including:
    /// - `generate_plates()` - Create plates using electrostatic physics
    /// - `roughen_boundaries()` - Add realistic irregularity
    /// - `deisland()` - Remove plate fragments
    /// - `analyze()` - Classify boundaries
    /// - `generate()` - Run complete pipeline
    /// - `export()` / `import()` - Save/load data
    ///
    /// # Example
    /// ```no_run
    /// use geoforge::WorldMap;
    ///
    /// let mut world = WorldMap::new(1800, 900, 42)?;
    ///
    /// // Complete pipeline:
    /// world.tectonics().generate(20)?;
    /// world.tectonics().export("outputs")?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn tectonics(&mut self) -> crate::map::tectonics_module::TectonicsModule<'_> {
        crate::map::tectonics_module::TectonicsModule::new(self)
    }

    /// Import tectonic plates with motion vectors from PNG (plate_motion.png)
    ///
    /// This imports both plate boundaries AND motion vectors from a single PNG file.
    /// Each unique color represents a plate with specific motion (direction and speed).
    ///
    /// The PNG must use the same HSV encoding as export_plate_motion_png:
    /// - Hue (0-360°) = Motion direction
    /// - Saturation (50-100%) = Motion speed (1-10 cm/year)
    ///
    /// This is the recommended way to import tectonic data from PNG files, as it
    /// preserves complete plate information including motion vectors.
    ///
    /// # Arguments
    /// * `png_path` - Path to plate_motion.png file
    ///
    /// # Example
    /// ```no_run
    /// # use geoforge::WorldMap;
    /// let mut world = WorldMap::new(1800, 900, 0)?;
    /// world.import_png("outputs/plate_motion.png")?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    #[cfg(feature = "export-png")]
    pub fn import_png(&mut self, png_path: &str) -> Result<(), Box<dyn std::error::Error>> {
        use image::io::Reader as ImageReader;
        use std::collections::HashMap;

        println!("🖼️ Importing tectonic plates with motion from PNG: {}", png_path);

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

        // Map colors to (plate_id, direction, speed)
        let mut color_to_plate: HashMap<[u8; 3], (u16, f64, f64)> = HashMap::new();
        let mut next_plate_id = 1u16;
        let mut tectonic_data = vec![0u16; self.width * self.height];

        // Process each pixel
        for y in 0..self.height {
            for x in 0..self.width {
                let pixel = rgb_img.get_pixel(x as u32, y as u32);
                let color = [pixel[0], pixel[1], pixel[2]];

                // Get or assign plate data for this color
                let (plate_id, _direction, _speed) = if let Some(&existing) = color_to_plate.get(&color) {
                    existing
                } else {
                    // New color - extract motion from it
                    let (direction, speed) = rgb_to_motion(color);
                    let new_id = next_plate_id;
                    color_to_plate.insert(color, (new_id, direction, speed));
                    next_plate_id += 1;

                    if next_plate_id > u16::MAX {
                        return Err("Too many unique colors in PNG (max 65535 plates)".into());
                    }

                    (new_id, direction, speed)
                };

                let idx = y * self.width + x;
                tectonic_data[idx] = plate_id;
            }
        }

        let num_plates = color_to_plate.len();
        println!("📊 Detected {} unique plates with motion vectors", num_plates);

        // Generate plate seeds from the imported data
        let mut plate_pixels: HashMap<u16, Vec<(usize, usize)>> = HashMap::new();

        for y in 0..self.height {
            for x in 0..self.width {
                let idx = y * self.width + x;
                let plate_id = tectonic_data[idx];

                if plate_id > 0 {
                    plate_pixels.entry(plate_id).or_insert_with(Vec::new).push((x, y));
                }
            }
        }

        let mut plate_seeds = Vec::new();

        // Extract motion info for each plate from the color mapping
        let mut plate_motion: HashMap<u16, (f64, f64)> = HashMap::new();
        for (_color, (plate_id, direction, speed)) in color_to_plate {
            plate_motion.insert(plate_id, (direction, speed));
        }

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

            // Get motion from color mapping
            let (motion_direction, motion_speed) = plate_motion.get(&plate_id)
                .copied()
                .unwrap_or((0.0, 2.5)); // Default fallback

            // Create seed with motion from PNG
            let seed = PlateSeed::new(
                plate_id,
                centroid_x,
                centroid_y,
                lat,
                lon,
                motion_direction,
                motion_speed,
            );

            plate_seeds.push(seed);
        }

        plate_seeds.sort_by_key(|s| s.id);

        println!("  Motion vectors imported:");
        for seed in plate_seeds.iter().take(5) {
            println!("    Plate {}: {:.0}° at {:.1} cm/year",
                     seed.id, seed.motion_direction, seed.motion_speed);
        }
        if plate_seeds.len() > 5 {
            println!("    ... and {} more plates", plate_seeds.len() - 5);
        }

        // Calculate statistics
        let tectonic_map = TerrainMap::from_data(self.width, self.height, tectonic_data);
        let plate_stats = self.calculate_plate_stats_from_map(&tectonic_map, &plate_seeds);

        // Create metadata and assign plate types
        let mut metadata = TectonicMetadata::new(plate_seeds, plate_stats);
        metadata.assign_plate_types();

        // Store results (motion already assigned from PNG)
        self.tectonics = Some(tectonic_map);
        self.tectonic_metadata = Some(metadata);

        println!("✅ Tectonic plates with motion vectors imported successfully");
        Ok(())
    }

    /// Generate tectonic plates using electrostatic physics simulation
    pub fn generate_tectonics(&mut self, num_plates: usize) -> Result<(), Box<dyn std::error::Error>> {
        // Use the dedicated TectonicPlateGenerator
        let mut generator = TectonicPlateGenerator::with_seed(
            self.width,
            self.height,
            num_plates,
            self.seed
        ).map_err(|e| -> Box<dyn std::error::Error> { Box::new(e) })?;

        // Generate the plates
        generator.generate(GenerationMethod::Electrostatic)
            .map_err(|e| -> Box<dyn std::error::Error> { Box::new(e) })?;

        // Extract data from generator
        let (width, height, plate_data, seeds) = generator.get_plate_data();
        let plate_stats = generator.get_plate_stats();

        // Create metadata and assign plate types
        let mut metadata = TectonicMetadata::new(seeds.clone(), plate_stats);
        metadata.assign_plate_types();

        // Assign realistic plate motion vectors
        let motion_config = PlateMotionConfig::with_seed(self.seed);
        let mut motion_assigner = PlateMotionAssigner::with_config(motion_config);
        motion_assigner.assign_motion(&mut metadata.plate_seeds);

        // Store results in WorldMap
        self.tectonics = Some(TerrainMap::from_data(width, height, plate_data.clone()));
        self.tectonic_metadata = Some(metadata);

        Ok(())
    }

    /// Apply boundary refinement to add realistic irregularity to plate edges (Stage 1.2)
    ///
    /// This must be called after `generate_tectonics()` to post-process the boundaries.
    /// If no configuration is provided, uses default parameters with the world seed.
    ///
    /// # Arguments
    /// * `config` - Optional configuration for boundary refinement. If None, uses defaults.
    ///
    /// # Example
    /// ```no_run
    /// use geoforge::{WorldMap, BoundaryRefinementConfig};
    ///
    /// let mut world = WorldMap::new(1800, 900, 42)?;
    /// world.generate_tectonics(15)?;
    ///
    /// // Apply boundary refinement with custom settings
    /// let config = BoundaryRefinementConfig::with_seed(42)
    ///     .with_noise(0.1, 3.0, 3);
    /// world.refine_boundaries(Some(config))?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn refine_boundaries(&mut self, config: Option<BoundaryRefinementConfig>) -> Result<(), Box<dyn std::error::Error>> {
        // Ensure tectonics have been generated
        if self.tectonics.is_none() {
            return Err("Tectonics must be generated before refining boundaries".into());
        }

        // Use provided config or create default with world seed
        let config = config.unwrap_or_else(|| BoundaryRefinementConfig::with_seed(self.seed));

        // Apply refinement
        let mut refiner = BoundaryRefiner::new(config);
        let tectonics = self.tectonics.as_mut().unwrap(); // Safe because we checked above
        refiner.refine_boundaries(tectonics);

        // Recalculate plate stats after refinement
        if self.tectonic_metadata.is_some() && self.tectonics.is_some() {
            let seeds = self.tectonic_metadata.as_ref().unwrap().plate_seeds.clone();
            let tectonics = self.tectonics.as_ref().unwrap();
            let new_stats = self.calculate_plate_stats_from_map(tectonics, &seeds);

            // Update metadata
            let metadata = self.tectonic_metadata.as_mut().unwrap();
            metadata.plate_stats = new_stats;
            metadata.assign_plate_types();
        }

        Ok(())
    }

    /// Remove plate islands to ensure all plates are contiguous (Stage 1.3)
    ///
    /// This should be called after `refine_boundaries()` to clean up any isolated plate fragments
    /// that may have been created by the boundary refinement process. It ensures that each plate
    /// consists of a single contiguous region by reassigning isolated fragments to surrounding plates.
    ///
    /// # Arguments
    /// * `config` - Optional configuration for island removal. If None, uses defaults.
    ///
    /// # Returns
    /// Statistics about the island removal process
    ///
    /// # Example
    /// ```no_run
    /// use geoforge::{WorldMap, IslandRemovalConfig};
    ///
    /// let mut world = WorldMap::new(1800, 900, 42)?;
    /// world.generate_tectonics(15)?;
    /// world.refine_boundaries(None)?;
    ///
    /// // Remove any plate islands created by boundary refinement
    /// let stats = world.remove_islands(None)?;
    /// println!("Removed {} islands from {} plates", stats.islands_removed, stats.plates_with_islands);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn remove_islands(&mut self, config: Option<IslandRemovalConfig>) -> Result<IslandRemovalStats, Box<dyn std::error::Error>> {
        // Ensure tectonics have been generated
        if self.tectonics.is_none() {
            return Err("Tectonics must be generated before removing islands".into());
        }

        // Use provided config or create default
        let config = config.unwrap_or_default();

        // Apply island removal
        let mut remover = IslandRemover::new(config);
        let tectonics = self.tectonics.as_mut().unwrap(); // Safe because we checked above
        let stats = remover.remove_islands(tectonics);

        // Recalculate plate stats after island removal
        if self.tectonic_metadata.is_some() && self.tectonics.is_some() {
            let seeds = self.tectonic_metadata.as_ref().unwrap().plate_seeds.clone();
            let tectonics = self.tectonics.as_ref().unwrap();
            let new_stats = self.calculate_plate_stats_from_map(tectonics, &seeds);

            // Update metadata
            let metadata = self.tectonic_metadata.as_mut().unwrap();
            metadata.plate_stats = new_stats;
            metadata.assign_plate_types();
        }

        Ok(stats)
    }

    /// Analyze plate boundaries and classify their types (Stage 1.4)
    ///
    /// This identifies all boundaries between plates and classifies them as convergent,
    /// divergent, or transform based on relative plate motion. Results are stored in
    /// the tectonic metadata.
    ///
    /// # Arguments
    /// * `config` - Optional configuration for boundary analysis. If None, uses defaults.
    ///
    /// # Returns
    /// Statistics about the analyzed boundaries
    ///
    /// # Example
    /// ```no_run
    /// use geoforge::WorldMap;
    ///
    /// let mut world = WorldMap::new(1800, 900, 42)?;
    /// world.generate_tectonics(15)?;
    /// world.refine_boundaries(None)?;
    /// world.remove_islands(None)?;
    ///
    /// // Analyze boundaries
    /// let boundary_stats = world.analyze_boundaries(None)?;
    /// println!("Found {} boundaries", boundary_stats.total_boundaries);
    /// println!("  Convergent: {}", boundary_stats.convergent_count);
    /// println!("  Divergent: {}", boundary_stats.divergent_count);
    /// println!("  Transform: {}", boundary_stats.transform_count);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn analyze_boundaries(
        &mut self,
        config: Option<BoundaryAnalysisConfig>
    ) -> Result<BoundaryStatistics, Box<dyn std::error::Error>> {
        // Ensure tectonics have been generated
        if self.tectonics.is_none() {
            return Err("Tectonics must be generated before analyzing boundaries".into());
        }

        if self.tectonic_metadata.is_none() {
            return Err("Tectonic metadata missing. Generate tectonics first.".into());
        }

        // Use provided config or create default
        let config = config.unwrap_or_default();

        // Create analyzer and find boundaries
        let analyzer = BoundaryAnalyzer::with_config(config);
        let tectonics = self.tectonics.as_ref().unwrap();
        let mut boundaries = analyzer.find_boundaries(tectonics);

        // Classify boundaries based on plate motion
        let metadata = self.tectonic_metadata.as_ref().unwrap();
        analyzer.classify_boundaries(&mut boundaries, &metadata.plate_seeds, tectonics);

        // Get statistics
        let stats = analyzer.boundary_statistics(&boundaries);

        // Store boundaries in metadata
        if let Some(ref mut metadata) = self.tectonic_metadata {
            metadata.plate_boundaries = boundaries;
            metadata.boundary_stats = Some(stats.clone());
        }

        Ok(stats)
    }

    /// Generate geological provinces from tectonic data (Stage 2.1: Orogenic Belts)
    ///
    /// This analyzes convergent plate boundaries and creates orogenic belts (mountain-building zones).
    /// Requires that tectonics have been generated and boundaries have been analyzed.
    ///
    /// # Arguments
    /// * `config` - Optional configuration for orogenic belt generation
    ///
    /// # Returns
    /// Vector of province regions representing orogenic belts
    ///
    /// # Example
    /// ```no_run
    /// use geoforge::{WorldMap, OrogenicConfig};
    ///
    /// let mut world = WorldMap::new(1800, 900, 42)?;
    /// world.tectonics().generate_plates(15)?;
    /// world.analyze_boundaries(None)?;
    ///
    /// // Generate orogenic belts with default config
    /// let orogens = world.generate_geology(None)?;
    /// println!("Generated {} orogenic belts", orogens.len());
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    /// Generate comprehensive geological provinces (Complete Stage 2)
    ///
    /// This generates ALL geological provinces from tectonic data:
    /// - Stage 2.1: Orogenic belts (collision, subduction, accretionary, extensional)
    /// - Stage 2.2: Large Igneous Provinces (flood basalts, plateaus, hotspots)
    /// - Stage 2.3: Arc and Basin Systems (volcanic arcs, forearc/backarc basins)
    /// - Stage 2.4: Stable Continental Regions (cratons, platforms, intracratonic basins)
    /// - Stage 2.5: Extensional Zones (continental rifts, extended crust)
    /// - Stage 2.6: Oceanic Domains (mid-ocean ridges, abyssal plains, trenches)
    ///
    /// # Arguments
    /// * `config` - Optional configuration (uses default if None)
    ///
    /// # Returns
    /// Vector of all geological provinces
    pub fn generate_geology(
        &mut self,
        config: Option<crate::geology::GeologyConfig>
    ) -> Result<Vec<ProvinceRegion>, Box<dyn std::error::Error>> {
        // Ensure tectonics exist
        let plate_map = self.tectonics.as_ref()
            .ok_or("Tectonics not generated. Call tectonics().generate_plates() first.")?;

        // Ensure metadata exists with boundaries
        let metadata = self.tectonic_metadata.as_ref()
            .ok_or("Tectonic metadata not found. Call tectonics().generate_plates() first.")?;

        if metadata.plate_boundaries.is_empty() {
            return Err("Boundaries not analyzed. Call analyze_boundaries() first.".into());
        }

        // Generate all geological provinces
        let config = config.unwrap_or_default();
        let generator = crate::geology::GeologyGenerator::new(
            config,
            self.seed,
            self.planetary_params.clone()
        );

        let provinces = generator.generate_all_provinces(
            &metadata.plate_boundaries,
            &metadata.plate_stats,
            &metadata.plate_seeds,
            plate_map,
        );

        // Store in geology field
        self.geology = Some(provinces.clone());

        Ok(provinces)
    }

    /// Get the generated geological provinces
    pub fn get_geology(&self) -> Option<&Vec<ProvinceRegion>> {
        self.geology.as_ref()
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
        self.tectonic_metadata.as_ref().map(|m| &m.plate_stats)
    }

    /// Get plate seeds for the tectonic layer
    pub fn get_tectonic_seeds(&self) -> Option<&Vec<PlateSeed>> {
        self.tectonic_metadata.as_ref().map(|m| &m.plate_seeds)
    }

    /// Get full tectonic metadata
    pub fn get_tectonic_metadata(&self) -> Option<&TectonicMetadata> {
        self.tectonic_metadata.as_ref()
    }



    /// Generate all available layers in sequence
    pub fn generate_all(&mut self, num_plates: usize) -> Result<(), Box<dyn std::error::Error>> {
        println!("🌍 Generating all world layers...");

        // Generate tectonics (foundation layer)
        println!("⚡ Generating tectonic layer...");
        self.generate_tectonics(num_plates)?;
        
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
            geology: None,
            elevation: None,
            temperature: None,
            precipitation: None,
            biomes: None,
            tectonic_metadata: None,
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

        // Read tectonic metadata if tectonic layer exists and there's more data
        if flags & 0x01 != 0 && cursor < buffer.len() {
            worldmap.tectonic_metadata = Some(TectonicMetadata::read_from_buffer(&buffer, &mut cursor)?);
        }

        println!("✅ WorldMap loaded from {}", filepath);
        Ok(worldmap)
    }


    
    /// Get flags indicating which layers are present
    pub fn get_layer_flags(&self) -> u32 {
        let mut flags = 0u32;
        if self.tectonics.is_some() { flags |= 0x01; }
        if self.elevation.is_some() { flags |= 0x02; }
        if self.temperature.is_some() { flags |= 0x04; }
        if self.precipitation.is_some() { flags |= 0x08; }
        if self.biomes.is_some() { flags |= 0x10; }
        flags
    }
}




/// Convert RGB color back to plate motion (direction and speed)
///
/// This is the inverse of motion_to_rgb, used for importing from plate_motion.png
#[cfg(feature = "export-png")]
fn rgb_to_motion(rgb: [u8; 3]) -> (f64, f64) {
    use crate::utils::color::rgb_to_hsv;

    let (hue, saturation, _value) = rgb_to_hsv(rgb[0], rgb[1], rgb[2]);

    // Direction: hue maps directly to degrees (0-360)
    let direction_deg = hue;

    // Speed: reverse the saturation mapping
    // saturation = 0.5 + (speed - 1.0) / 9.0 * 0.5
    // saturation - 0.5 = (speed - 1.0) / 9.0 * 0.5
    // (saturation - 0.5) * 2.0 = (speed - 1.0) / 9.0
    // (saturation - 0.5) * 2.0 * 9.0 = speed - 1.0
    // speed = (saturation - 0.5) * 18.0 + 1.0
    let speed_cm_year = (saturation - 0.5) * 18.0 + 1.0;
    let speed_cm_year = speed_cm_year.clamp(1.0, 10.0);

    (direction_deg, speed_cm_year)
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::export::MapExporter;

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
        let result = worldmap.generate_tectonics(5);
        assert!(result.is_ok());
        assert!(worldmap.tectonics.is_some());
        assert!(worldmap.tectonic_metadata.is_some());

        let metadata = worldmap.tectonic_metadata.as_ref().unwrap();
        assert!(!metadata.plate_seeds.is_empty());
        assert!(!metadata.plate_stats.is_empty());
    }
    
    #[test]
    fn test_serialization_roundtrip() {
        let mut worldmap = WorldMap::new(100, 50, 42).unwrap();
        worldmap.generate_tectonics(5).unwrap();
        
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
        
        worldmap1.generate_tectonics(5).unwrap();
        worldmap2.generate_tectonics(5).unwrap();
        
        // Compare tectonic data
        let data1 = &worldmap1.tectonics.as_ref().unwrap().data;
        let data2 = &worldmap2.tectonics.as_ref().unwrap().data;
        assert_eq!(data1, data2);
        
        // Compare plate seeds
        let seeds1 = &worldmap1.tectonic_metadata.as_ref().unwrap().plate_seeds;
        let seeds2 = &worldmap2.tectonic_metadata.as_ref().unwrap().plate_seeds;
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
        assert!(worldmap.generate_tectonics(0).is_err());
        assert!(worldmap.generate_tectonics(u16::MAX as usize + 1).is_err());
    }
    
    #[test]
    fn test_plate_size_ratios() {
        let mut worldmap = WorldMap::new(200, 100, 42).unwrap();
        worldmap.generate_tectonics(20).unwrap();
        
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
        original_world.generate_tectonics(6).unwrap();

        let png_path = format!("{}/test_motion.png", test_dir);
        original_world.export_plate_motion_png(test_dir, "test_motion.png").unwrap();

        // Step 2: Import PNG into new world object
        let mut imported_world = WorldMap::new(120, 60, 999).unwrap(); // Different seed
        imported_world.import_png(&png_path).unwrap();
        
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
        let result = world.generate_all(4);
        assert!(result.is_ok());
        
        // Should have generated tectonics
        assert!(world.tectonics.is_some());
        assert!(world.tectonic_metadata.is_some());
        
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