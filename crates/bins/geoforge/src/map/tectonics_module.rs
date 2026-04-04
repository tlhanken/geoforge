//! Tectonic operations module for WorldMap
//!
//! This module provides a clean API for all tectonic plate operations.

use crate::map::world::WorldMap;
use crate::tectonics::{
    BoundaryAnalysisConfig, BoundaryStatistics, IslandRemovalConfig, IslandRemovalStats,
    BoundaryRefinementConfig,
};
use crate::io::export::MapExporter;
use std::error::Error;

/// Tectonic operations interface for WorldMap
///
/// Provides organized access to all tectonic plate generation, refinement,
/// and analysis operations.
///
/// # Example
/// ```no_run
/// use geoforge::WorldMap;
///
/// let mut world = WorldMap::new(1800, 900, 42)?;
///
/// // Complete pipeline:
/// world.tectonics().generate(20)?;
///
/// // Or individual steps:
/// world.tectonics().generate_plates(20)?;
/// world.tectonics().roughen_boundaries(None)?;
/// world.tectonics().deisland(None)?;
/// world.tectonics().analyze(None)?;
///
/// // Export all visualizations:
/// world.tectonics().export("outputs")?;
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub struct TectonicsModule<'a> {
    world: &'a mut WorldMap,
}

impl<'a> TectonicsModule<'a> {
    /// Create a new TectonicsModule (internal use only)
    pub(crate) fn new(world: &'a mut WorldMap) -> Self {
        Self { world }
    }

    /// Generate tectonic plates using electrostatic physics simulation
    ///
    /// This is Stage 1.1 of the pipeline. Creates realistic tectonic plates with
    /// Earth-like size distribution using Coulomb's law physics.
    ///
    /// # Arguments
    /// * `plate_count` - Number of plates to generate (typically 10-30 for Earth-like worlds)
    ///
    /// # Returns
    /// * `Ok(())` - Plates generated successfully
    /// * `Err(_)` - Generation failed
    ///
    /// # Example
    /// ```no_run
    /// # use geoforge::WorldMap;
    /// # let mut world = WorldMap::new(360, 180, 42)?;
    /// world.tectonics().generate_plates(15)?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn generate_plates(&mut self, plate_count: usize) -> Result<(), Box<dyn Error>> {
        self.world.generate_tectonics(plate_count)
    }

    /// Roughen plate boundaries to add realistic irregularity
    ///
    /// This is Stage 1.2 of the pipeline. Uses 3D spherical domain warping with
    /// multi-octave Perlin noise to create natural-looking jagged boundaries.
    ///
    /// # Arguments
    /// * `config` - Optional configuration (None uses defaults)
    ///
    /// # Returns
    /// * `Ok(())` - Boundaries refined successfully
    /// * `Err(_)` - Refinement failed or no plates exist
    ///
    /// # Example
    /// ```no_run
    /// # use geoforge::WorldMap;
    /// # let mut world = WorldMap::new(360, 180, 42)?;
    /// # world.tectonics().generate_plates(10)?;
    /// world.tectonics().roughen_boundaries(None)?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn roughen_boundaries(
        &mut self,
        config: Option<BoundaryRefinementConfig>,
    ) -> Result<(), Box<dyn Error>> {
        self.world.refine_boundaries(config)
    }

    /// Remove plate islands to ensure contiguous plates
    ///
    /// This is Stage 1.3 of the pipeline. Uses flood-fill to identify disconnected
    /// plate fragments and reassigns them to neighboring plates.
    ///
    /// # Arguments
    /// * `config` - Optional configuration (None uses defaults)
    ///
    /// # Returns
    /// * `Ok(IslandRemovalStatistics)` - Statistics about islands removed
    /// * `Err(_)` - Island removal failed or no plates exist
    ///
    /// # Example
    /// ```no_run
    /// # use geoforge::WorldMap;
    /// # let mut world = WorldMap::new(360, 180, 42)?;
    /// # world.tectonics().generate_plates(10)?;
    /// # world.tectonics().roughen_boundaries(None)?;
    /// let stats = world.tectonics().deisland(None)?;
    /// println!("Removed {} islands", stats.islands_removed);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn deisland(
        &mut self,
        config: Option<IslandRemovalConfig>,
    ) -> Result<IslandRemovalStats, Box<dyn Error>> {
        self.world.remove_islands(config)
    }

    /// Analyze plate boundaries and classify by interaction type
    ///
    /// This is Stage 1.4 of the pipeline. Identifies boundaries, calculates relative
    /// motion, and classifies as convergent, divergent, or transform.
    ///
    /// # Arguments
    /// * `config` - Optional configuration (None uses defaults)
    ///
    /// # Returns
    /// * `Ok(BoundaryStatistics)` - Statistics about boundaries found
    /// * `Err(_)` - Analysis failed or no plates exist
    ///
    /// # Example
    /// ```no_run
    /// # use geoforge::WorldMap;
    /// # let mut world = WorldMap::new(360, 180, 42)?;
    /// # world.tectonics().generate_plates(10)?;
    /// # world.tectonics().roughen_boundaries(None)?;
    /// # world.tectonics().deisland(None)?;
    /// let stats = world.tectonics().analyze(None)?;
    /// println!("Found {} convergent boundaries", stats.convergent_count);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn analyze(
        &mut self,
        config: Option<BoundaryAnalysisConfig>,
    ) -> Result<BoundaryStatistics, Box<dyn Error>> {
        self.world.analyze_boundaries(config)
    }

    /// Run complete tectonic generation pipeline (all 4 stages)
    ///
    /// Executes the full Stage 1 pipeline:
    /// 1. Generate plates (electrostatic physics)
    /// 2. Roughen boundaries (domain warping)
    /// 3. Remove islands (contiguity)
    /// 4. Analyze boundaries (classification)
    ///
    /// # Arguments
    /// * `plate_count` - Number of plates to generate
    ///
    /// # Returns
    /// * `Ok(BoundaryStatistics)` - Final boundary statistics
    /// * `Err(_)` - Pipeline failed at some stage
    ///
    /// # Example
    /// ```no_run
    /// # use geoforge::WorldMap;
    /// let mut world = WorldMap::new(1800, 900, 42)?;
    /// let stats = world.tectonics().generate(20)?;
    /// println!("Generated {} boundaries", stats.total_boundaries);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn generate(&mut self, plate_count: usize) -> Result<BoundaryStatistics, Box<dyn Error>> {
        // Stage 1.1: Generate plates
        self.generate_plates(plate_count)?;

        // Stage 1.2: Roughen boundaries
        self.roughen_boundaries(None)?;

        // Stage 1.3: Remove islands
        self.deisland(None)?;

        // Stage 1.4: Analyze boundaries
        let stats = self.analyze(None)?;

        Ok(stats)
    }

    /// Export all tectonic visualizations and data
    ///
    /// Exports:
    /// - `tectonics.png` - Color-coded plates
    /// - `tectonics_boundaries.png` - Boundary types (red/blue/green)
    /// - `tectonics_motion.png` - Motion vectors (hue=direction, saturation=speed)
    /// - `world.map` - Complete binary world data
    ///
    /// # Arguments
    /// * `output_dir` - Directory to save files (will be created if needed)
    ///
    /// # Returns
    /// * `Ok(())` - All files exported successfully
    /// * `Err(_)` - Export failed
    ///
    /// # Example
    /// ```no_run
    /// # use geoforge::WorldMap;
    /// # let mut world = WorldMap::new(360, 180, 42)?;
    /// # world.tectonics().generate(10)?;
    /// world.tectonics().export("outputs")?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    #[cfg(feature = "export-png")]
    pub fn export(&self, output_dir: &str) -> Result<(), Box<dyn Error>> {
        use std::fs;

        // Create output directory
        fs::create_dir_all(output_dir)?;

        // Export PNGs
        self.world.export_tectonics_png(output_dir, "tectonics.png")?;
        self.world.export_boundaries_png(output_dir, "tectonics_boundaries.png")?;
        self.world.export_plate_motion_png(output_dir, "tectonics_motion.png")?;

        // Export binary
        let map_path = format!("{}/world.map", output_dir);
        self.world.save_to_file(&map_path)?;

        Ok(())
    }

    /// Export all tectonic visualizations and data (no PNG support)
    #[cfg(not(feature = "export-png"))]
    pub fn export(&self, output_dir: &str) -> Result<(), Box<dyn Error>> {
        use std::fs;

        // Create output directory
        fs::create_dir_all(output_dir)?;

        // Export only binary
        let map_path = format!("{}/world.map", output_dir);
        self.world.save_to_file(&map_path)?;

        Ok(())
    }

    /// Import tectonic data from binary file
    ///
    /// Loads a previously saved world.map file, restoring all tectonic data
    /// including plates, boundaries, motion vectors, and statistics.
    ///
    /// # Arguments
    /// * `path` - Path to world.map file
    ///
    /// # Returns
    /// * `Ok(())` - Data loaded successfully
    /// * `Err(_)` - Import failed
    ///
    /// # Example
    /// ```no_run
    /// # use geoforge::WorldMap;
    /// let mut world = WorldMap::new(1800, 900, 0)?;
    /// world.tectonics().import("outputs/world.map")?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn import(&mut self, path: &str) -> Result<(), Box<dyn Error>> {
        // Load the world from file - this replaces the entire world state
        let loaded_world = WorldMap::load_from_file(path)?;

        // Replace the current world's data
        *self.world = loaded_world;

        Ok(())
    }

    /// Import tectonic plates with motion vectors from PNG
    ///
    /// Loads both plate boundaries AND motion vectors from a single PNG file (tectonics_motion.png).
    /// Each unique color encodes both the plate identity and its motion:
    /// - Hue (color) = Motion direction (0-360°)
    /// - Saturation (brightness) = Motion speed (1-10 cm/year)
    ///
    /// This is the standard way to import tectonic data from PNG files.
    ///
    /// # Arguments
    /// * `path` - Path to tectonics_motion.png file
    ///
    /// # Returns
    /// * `Ok(())` - Plates and motion imported successfully
    /// * `Err(_)` - Import failed
    ///
    /// # Example
    /// ```no_run
    /// # use geoforge::WorldMap;
    /// # let mut world = WorldMap::new(1800, 900, 0)?;
    /// // Import complete tectonic data including motion vectors
    /// world.tectonics().import_png("outputs/tectonics_motion.png")?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    #[cfg(feature = "export-png")]
    pub fn import_png(&mut self, path: &str) -> Result<(), Box<dyn Error>> {
        self.world.import_png(path)
    }
}
