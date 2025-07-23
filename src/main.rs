use geoforge::WorldMap;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🌍 Geoforge - Realistic World Generation");
    println!("=====================================");

    // Check for PNG import argument
    let args: Vec<String> = std::env::args().collect();
    if args.len() == 3 && args[1] == "--import-png" {
        return import_png_mode(&args[2]);
    }

    // Create a new world map
    let seed = 097243067;
    println!("\n🗺️ Creating new world map (1900x900, seed: {})...", seed);
    let mut world = WorldMap::new(1800, 900, seed)?;
    
    // Generate tectonic plates using electrostatic physics
    println!("\n⚡ Generating tectonic layer...");
    world.generate_tectonics(20, false)?;
    
    // Generate geology domains based on tectonics
    println!("\n🗻 Generating geology layer...");
    world.generate_geology()?;
    
    // Show statistics
    if let Some(stats) = world.get_tectonic_stats() {
        println!("\n📊 Generated {} plates with realistic size distribution:", stats.len());
        
        // Sort by size for better display
        let mut sorted_stats: Vec<_> = stats.iter().collect();
        sorted_stats.sort_by(|a, b| b.1.area_km2.partial_cmp(&a.1.area_km2).unwrap());
        
        for (i, (plate_id, stat)) in sorted_stats.iter().enumerate().take(5) {
            let category = match i {
                0 => "Superplate",
                1..=2 => "Major plate",
                3..=4 => "Medium plate",
                _ => "Small plate",
            };
            println!("  {} {}: {:.1}% of surface ({:.0} km²)", 
                     category, plate_id, stat.percentage, stat.area_km2);
        }
        
        if sorted_stats.len() > 5 {
            println!("  ... and {} smaller plates", sorted_stats.len() - 5);
        }
    }
    
    // Export individual layers
    println!("\n💾 Exporting layers...");
    std::fs::create_dir_all("outputs")?;
    
    #[cfg(feature = "export-png")]
    {
        world.export_all_png("outputs", "world")?;
        println!("✅ All layers exported as PNG files");
    }
    
    #[cfg(not(feature = "export-png"))]
    {
        println!("⚠️ PNG export disabled. Use --features export-png for visualization files");
    }
    
    // Save complete world map to binary file
    world.save_to_file("outputs/world.map")?;
    println!("✅ Complete world map saved to binary file");
    
    println!("\n🎉 WORLD GENERATION COMPLETE!");
    println!("Files created in outputs/ directory:");
    println!("  • world.map - Complete world data (binary)");
    
    #[cfg(feature = "export-png")]
    {
        println!("  • world_tectonics.png - Tectonic plates visualization");
        println!("  • world_geology.png - Geological domains visualization");
    }
    
    Ok(())
}

#[cfg(feature = "export-png")]
fn import_png_mode(png_path: &str) -> Result<(), Box<dyn std::error::Error>> {
    println!("🖼️ PNG Import Mode");
    println!("=================");
    
    // Try to determine dimensions from the PNG file
    use image::io::Reader as ImageReader;
    use image::GenericImageView;
    let img = ImageReader::open(png_path)?.decode()?;
    let (width, height) = img.dimensions();
    
    println!("📐 Detected PNG dimensions: {}×{}", width, height);
    
    // Create world map with matching dimensions
    let mut world = WorldMap::new(width as usize, height as usize, 0)?;
    
    // Import the PNG
    world.import_tectonics_png(png_path)?;
    
    // Show statistics
    if let Some(stats) = world.get_tectonic_stats() {
        println!("\n📊 Imported {} plates:", stats.len());
        
        let mut sorted_stats: Vec<_> = stats.iter().collect();
        sorted_stats.sort_by(|a, b| b.1.area_km2.cmp(&a.1.area_km2));
        
        for (i, (plate_id, stat)) in sorted_stats.iter().enumerate().take(10) {
            let category = match i {
                0 => "🌍 Largest",
                1 => "🏔️  2nd largest",
                2 => "⛰️  3rd largest",
                _ => "🗻 Plate",
            };
            println!("  {} {}: {:.1}% ({} km²)", 
                     category, plate_id, stat.percentage, stat.area_km2);
        }
        
        if sorted_stats.len() > 10 {
            println!("  ... and {} more plates", sorted_stats.len() - 10);
        }
    }
    
    // Export results
    std::fs::create_dir_all("outputs/imported")?;
    world.save_to_file("outputs/imported/world.map")?;
    println!("\n💾 Saved imported world to: outputs/imported/world.map");
    
    world.export_tectonics_png("outputs/imported", "plates.png")?;
    println!("🎨 Exported visualization: outputs/imported/plates.png");
    
    println!("\n🎉 PNG import completed successfully!");
    println!("\n💡 Usage: geoforge --import-png <path-to-png>");
    
    Ok(())
}

#[cfg(not(feature = "export-png"))]
fn import_png_mode(_png_path: &str) -> Result<(), Box<dyn std::error::Error>> {
    println!("❌ PNG import requires --features export-png");
    println!("   Run: cargo run --features export-png -- --import-png <path>");
    Ok(())
}