use geoforge::TectonicPlateGenerator;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🌍 Geoforge - Tectonic Plate Generation");
    println!("======================================");

    // Create a test world
    let mut generator = TectonicPlateGenerator::with_seed(3600, 1800, 15, 42)?;
    println!("✅ Generator created with seed: {}", generator.get_seed());
    
    // Generate plates
    generator.generate("region_growing", true)?;
    println!("✅ Plates generated successfully");
    
    // Test export methods
    generator.export_binary("outputs", "test.bin")?;
    println!("✅ Binary export successful");
    
    generator.export_all("outputs", "test_world")?;
    println!("✅ Export all successful");
    
    // Show statistics
    let stats = generator.get_plate_stats();
    println!("✅ Generated {} plates with statistics", stats.len());
    
    for (plate_id, stat) in stats.iter().take(3) {
        println!("  Plate {}: {:.1}% of surface ({} km²)", 
                 plate_id, stat.percentage, stat.area_km2);
    }
    
    println!("\n🎉 ALL TESTS PASSED!");
    println!("Files created in outputs/ directory:");
    println!("  • test.bin");
    println!("  • test_world.bin");
    
    Ok(())
}