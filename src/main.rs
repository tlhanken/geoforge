use geoforge::{WorldMap, LayerType};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🌍 Geoforge - Realistic World Generation");
    println!("=====================================");

    // Create a new world map
    let seed = 097243067;
    println!("\n🗺️ Creating new world map (1800x900, seed: {})...", seed);
    let mut world = WorldMap::new(1800, 900, seed)?;
    
    // Generate tectonic plates using electrostatic physics
    println!("\n⚡ Generating tectonic layer...");
    world.generate_tectonics(20, true)?;
    
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
        world.export_layer_png(LayerType::Tectonics, "outputs", "tectonics.png")?;
        println!("✅ Tectonics layer exported as PNG");
    }
    
    // Save complete world map to binary file
    world.save_to_file("outputs/world.map")?;
    println!("✅ Complete world map saved to binary file");
    
    // Demonstrate loading
    println!("\n🔄 Testing world map loading...");
    let loaded_world = WorldMap::load_from_file("outputs/world.map")?;
    println!("✅ World map loaded successfully ({}x{}, seed: {})", 
             loaded_world.width, loaded_world.height, loaded_world.seed);
    
    println!("\n🎉 WORLD GENERATION COMPLETE!");
    println!("Files created in outputs/ directory:");
    println!("  • world.map - Complete world data (binary)");
    
    #[cfg(feature = "export-png")]
    println!("  • tectonics.png - Tectonic plates visualization");
    
    Ok(())
}