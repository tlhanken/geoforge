use crate::map::world::WorldMap;
#[cfg(feature = "export-png")]
use crate::geology::provinces::GeologicProvince;
#[cfg(feature = "export-png")]
use image::{ImageBuffer, Rgb};
use std::error::Error;
use std::fs;
use std::path::Path;


/// Trait for exporting map data to various formats
pub trait MapExporter {
    /// Save the entire world state to a binary .map file
    fn save_to_file(&self, filepath: &str) -> Result<(), Box<dyn Error>>;

    /// Export tectonic plates visualization
    #[cfg(feature = "export-png")]
    fn export_tectonics_png(&self, output_dir: &str, filename: &str) -> Result<(), Box<dyn Error>>;

    /// Export boundary visualization
    #[cfg(feature = "export-png")]
    fn export_boundaries_png(&self, output_dir: &str, filename: &str) -> Result<(), Box<dyn Error>>;

    /// Export plate motion visualization
    #[cfg(feature = "export-png")]
    fn export_plate_motion_png(&self, output_dir: &str, filename: &str) -> Result<(), Box<dyn Error>>;

    /// Export motion reference visualization
    #[cfg(feature = "export-png")]
    fn export_motion_reference_png(&self, output_dir: &str, filename: &str) -> Result<(), Box<dyn Error>>;

    /// Export plate type visualization (Oceanic vs Continental)
    #[cfg(feature = "export-png")]
    fn export_plate_types_png(&self, output_dir: &str, filename: &str) -> Result<(), Box<dyn Error>>;

    /// Export geological provinces visualization (Stage 2)
    #[cfg(feature = "export-png")]
    fn export_geology_png(&self, output_dir: &str, filename: &str) -> Result<(), Box<dyn Error>>;

    /// Export geological provinces with boundary overlays
    #[cfg(feature = "export-png")]
    fn export_geology_with_boundaries_png(&self, output_dir: &str, filename: &str) -> Result<(), Box<dyn Error>>;

    /// Export all available visualizations
    #[cfg(feature = "export-png")]
    fn export_all_png(&self, output_dir: &str, base_name: &str) -> Result<(), Box<dyn Error>>;
}

impl MapExporter for WorldMap {
    fn save_to_file(&self, filepath: &str) -> Result<(), Box<dyn Error>> {
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
        
        // Write complete tectonic metadata if available
        if let Some(ref metadata) = self.tectonic_metadata {
            crate::map::world::TectonicMetadata::write_to_file(&mut file, metadata)?;
        }

        println!("✅ WorldMap saved to {}", filepath);
        Ok(())
    }

    #[cfg(feature = "export-png")]
    fn export_tectonics_png(&self, output_dir: &str, filename: &str) -> Result<(), Box<dyn Error>> {
        use rand::prelude::*;
        
        fs::create_dir_all(output_dir)?;
        let path = Path::new(output_dir).join(filename);
        
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

    #[cfg(feature = "export-png")]
    fn export_boundaries_png(&self, output_dir: &str, filename: &str) -> Result<(), Box<dyn Error>> {
        fs::create_dir_all(output_dir)?;
        let path = Path::new(output_dir).join(filename);

        if let (Some(ref tectonic_map), Some(ref metadata)) = (&self.tectonics, &self.tectonic_metadata) {
            if metadata.plate_boundaries.is_empty() {
                return Err("Boundaries not analyzed. Call analyze_boundaries() first.".into());
            }

            let mut img = ImageBuffer::new(self.width as u32, self.height as u32);

            // Background: plates in grayscale
            for (y, row) in tectonic_map.data.chunks(self.width).enumerate() {
                for (x, &plate_id) in row.iter().enumerate() {
                    // Vary grayscale based on plate ID for subtle differentiation
                    let gray_value = (((plate_id as f64) * 37.5).rem_euclid(128.0) + 64.0) as u8;
                    img.put_pixel(x as u32, y as u32, Rgb([gray_value, gray_value, gray_value]));
                }
            }

            // Overlay: boundaries colored by type
            for boundary in &metadata.plate_boundaries {
                use crate::tectonics::plates::PlateInteraction;

                let color = match boundary.interaction_type {
                    PlateInteraction::Convergent => [255, 0, 0],     // Red
                    PlateInteraction::Divergent => [0, 128, 255],    // Blue
                    PlateInteraction::Transform => [0, 255, 0],      // Green
                };

                for (x, y) in &boundary.pixels {
                    img.put_pixel(*x as u32, *y as u32, Rgb(color));
                }
            }

            img.save(path)?;
        } else {
            return Err("Tectonic layer not generated or metadata missing".into());
        }

        Ok(())
    }

    #[cfg(feature = "export-png")]
    fn export_plate_motion_png(&self, output_dir: &str, filename: &str) -> Result<(), Box<dyn Error>> {
        fs::create_dir_all(output_dir)?;
        let path = Path::new(output_dir).join(filename);

        if let (Some(ref tectonic_map), Some(ref metadata)) = (&self.tectonics, &self.tectonic_metadata) {
            let mut img = ImageBuffer::new(self.width as u32, self.height as u32);

            // Create motion color map for each plate
            let mut plate_colors: std::collections::HashMap<u16, [u8; 3]> = std::collections::HashMap::new();

            for seed in &metadata.plate_seeds {
                let color = motion_to_rgb(seed.motion_direction, seed.motion_speed);
                plate_colors.insert(seed.id, color);
            }

            // Fill image with plate colors
            for (y, row) in tectonic_map.data.chunks(self.width).enumerate() {
                for (x, &plate_id) in row.iter().enumerate() {
                    if plate_id > 0 {
                        if let Some(&color) = plate_colors.get(&plate_id) {
                            img.put_pixel(x as u32, y as u32, Rgb(color));
                        }
                    } else {
                        img.put_pixel(x as u32, y as u32, Rgb([0, 0, 0]));
                    }
                }
            }

            img.save(path)?;
        } else {
            return Err("Tectonic layer not generated or metadata missing".into());
        }

        Ok(())
    }

    #[cfg(feature = "export-png")]
    fn export_motion_reference_png(&self, output_dir: &str, filename: &str) -> Result<(), Box<dyn Error>> {
        fs::create_dir_all(output_dir)?;
        let path = Path::new(output_dir).join(filename);

        // Create 400x400 reference image
        let size = 400;
        let center = size as f64 / 2.0;
        let radius = 150.0;

        let mut img = ImageBuffer::new(size, size);

        // Fill with white background
        for y in 0..size {
            for x in 0..size {
                img.put_pixel(x, y, Rgb([255, 255, 255]));
            }
        }

        // Draw color wheel
        for y in 0..size {
            for x in 0..size {
                let dx = x as f64 - center;
                let dy = y as f64 - center;
                let dist = (dx * dx + dy * dy).sqrt();

                // Only draw within circle
                if dist < radius && dist > 50.0 {
                    // Calculate angle (0° = North, clockwise)
                    let angle = dy.atan2(dx).to_degrees() + 90.0;
                    let direction = if angle < 0.0 { angle + 360.0 } else { angle };

                    // Use medium speed for uniform saturation
                    let color = motion_to_rgb(direction, 5.0);
                    img.put_pixel(x, y, Rgb(color));
                }
            }
        }

        img.save(&path)?;
        Ok(())
    }

    #[cfg(feature = "export-png")]
    fn export_plate_types_png(&self, output_dir: &str, filename: &str) -> Result<(), Box<dyn Error>> {
        use rand::prelude::*;
        use crate::tectonics::plates::PlateType;

        fs::create_dir_all(output_dir)?;
        let path = Path::new(output_dir).join(filename);

        if let (Some(ref tectonic_map), Some(ref metadata)) = (&self.tectonics, &self.tectonic_metadata) {
            let mut img = ImageBuffer::new(self.width as u32, self.height as u32);

            // Generate colors for each plate based on its type
            let mut plate_colors: std::collections::HashMap<u16, [u8; 3]> = std::collections::HashMap::new();
            let mut rng = StdRng::seed_from_u64(self.seed);

            for (plate_id, stats) in &metadata.plate_stats {
                let color = match stats.plate_type {
                    // Oceanic: Cool blue/cyan
                    PlateType::Oceanic => {
                        let hue = rng.gen_range(190..230) as f64;
                        let sat = rng.gen_range(70..100) as f64 / 100.0;
                        let val = rng.gen_range(60..95) as f64 / 100.0;
                        crate::utils::color::hsv_to_rgb(hue, sat, val)
                    }
                    // Continental: Warm red/orange
                    PlateType::Continental => {
                        let hue = rng.gen_range(0..50) as f64; // Red to orange (avoiding yellow-green)
                        let sat = rng.gen_range(70..100) as f64 / 100.0; // High saturation
                        let val = rng.gen_range(60..95) as f64 / 100.0;  // Bright
                        crate::utils::color::hsv_to_rgb(hue, sat, val)
                    }
                };
                plate_colors.insert(*plate_id, color);
            }

            // Paint the map
            for (y, row) in tectonic_map.data.chunks(self.width).enumerate() {
                for (x, &plate_id) in row.iter().enumerate() {
                    if plate_id > 0 {
                        if let Some(&color) = plate_colors.get(&plate_id) {
                            img.put_pixel(x as u32, y as u32, Rgb(color));
                        }
                    } else {
                        img.put_pixel(x as u32, y as u32, Rgb([0, 0, 0]));
                    }
                }
            }

            img.save(path)?;
        } else {
            return Err("Tectonic layer not generated or metadata missing".into());
        }

        Ok(())
    }

    #[cfg(feature = "export-png")]
    fn export_geology_png(&self, output_dir: &str, filename: &str) -> Result<(), Box<dyn Error>> {
        fs::create_dir_all(output_dir)?;
        let path = Path::new(output_dir).join(filename);

        if let Some(ref geology) = self.geology {
            let mut img = ImageBuffer::from_pixel(
                self.width as u32,
                self.height as u32,
                Rgb([255, 255, 255])
            );

            for region in geology {
                let color = get_province_color(region.characteristics.province_type);

                for &(x, y) in &region.pixels {
                    if x < self.width && y < self.height {
                        img.put_pixel(x as u32, y as u32, Rgb(color));
                    }
                }
            }

            img.save(path)?;
        } else {
             return Err("Geological provinces not generated. Call generate_geology() first.".into());
        }

        Ok(())
    }

    #[cfg(feature = "export-png")]
    fn export_geology_with_boundaries_png(&self, output_dir: &str, filename: &str) -> Result<(), Box<dyn Error>> {
        fs::create_dir_all(output_dir)?;
        let path = Path::new(output_dir).join(filename);

        if let (Some(ref geology), Some(ref metadata)) = (&self.geology, &self.tectonic_metadata) {
            let mut img = ImageBuffer::from_pixel(
                self.width as u32,
                self.height as u32,
                Rgb([255, 255, 255])
            );

            // First layer: Geology
            for region in geology {
                let color = get_province_color(region.characteristics.province_type);

                for &(x, y) in &region.pixels {
                    if x < self.width && y < self.height {
                        img.put_pixel(x as u32, y as u32, Rgb(color));
                    }
                }
            }

            // Second layer: Boundaries
            if !metadata.plate_boundaries.is_empty() {
                for boundary in &metadata.plate_boundaries {
                    use crate::tectonics::plates::PlateInteraction;

                    let color = match boundary.interaction_type {
                        PlateInteraction::Convergent => [255, 0, 0],     // Red
                        PlateInteraction::Divergent => [0, 128, 255],    // Blue
                        PlateInteraction::Transform => [0, 255, 0],      // Green
                    };

                    for (x, y) in &boundary.pixels {
                        if *x < self.width && *y < self.height {
                            img.put_pixel(*x as u32, *y as u32, Rgb(color));
                        }
                    }
                }
            }

            img.save(path)?;
        } else {
             return Err("Missing geology or tectonic metadata".into());
        }

        Ok(())
    }

    #[cfg(feature = "export-png")]
    fn export_all_png(&self, output_dir: &str, base_name: &str) -> Result<(), Box<dyn Error>> {
        fs::create_dir_all(output_dir)?;
        
        if self.tectonics.is_some() {
            self.export_tectonics_png(output_dir, &format!("{}_tectonics.png", base_name))?;
            self.export_boundaries_png(output_dir, &format!("{}_boundaries.png", base_name))?;
            self.export_plate_motion_png(output_dir, &format!("{}_plate_motion.png", base_name))?;
            self.export_plate_types_png(output_dir, &format!("{}_plate_types.png", base_name))?;
        }

        if self.geology.is_some() {
            self.export_geology_png(output_dir, &format!("{}_geology.png", base_name))?;
            self.export_geology_with_boundaries_png(output_dir, &format!("{}_geology_boundaries.png", base_name))?;
        }

        Ok(())
    }
}

// Helpers
#[cfg(feature = "export-png")]
fn motion_to_rgb(direction_deg: f64, speed_cm_year: f64) -> [u8; 3] {
    use crate::utils::color::hsv_to_rgb;

    let hue = direction_deg;
    let saturation = 0.5 + (speed_cm_year - 1.0) / 9.0 * 0.5;
    let saturation = saturation.clamp(0.5, 1.0);
    let value = 0.9;

    hsv_to_rgb(hue, saturation, value)
}

#[cfg(feature = "export-png")]
fn get_province_color(province_type: GeologicProvince) -> [u8; 3] {
    match province_type {
        GeologicProvince::CollisionOrogen => [50, 150, 80],
        GeologicProvince::AccretionaryWedge => [255, 200, 100],
        
        GeologicProvince::ContinentalFloodBasalt => [150, 50, 150],
        GeologicProvince::OceanicPlateau => [180, 100, 180],
        GeologicProvince::HotspotTrack => [200, 120, 200],
        
        GeologicProvince::VolcanicArc => [255, 50, 50],
        GeologicProvince::ForearcBasin => [160, 160, 160],
        GeologicProvince::BackarcBasin => [180, 180, 180],
        
        GeologicProvince::Craton => [255, 140, 60],
        GeologicProvince::Platform => [255, 150, 200],
        GeologicProvince::IntracratonicBasin => [160, 120, 180],
        
        GeologicProvince::ContinentalRift => [255, 220, 60],
        GeologicProvince::ExtendedCrust => [255, 240, 120],
        
        GeologicProvince::MidOceanRidge => [100, 200, 200],
        GeologicProvince::AbyssalPlain => [70, 130, 180],
        GeologicProvince::OceanTrench => [0, 50, 100],
        GeologicProvince::OceanicFractureZone => [120, 180, 170],
        GeologicProvince::OceanicHotspotTrack => [100, 80, 180],
        GeologicProvince::ContinentalHotspotTrack => [120, 60, 200],
    }
}
