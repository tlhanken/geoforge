use geoforge::TectonicPlateGenerator;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🌍 Geoforge - Realistic Tectonic Plate Generation");
    println!("===============================================");

    // Generate a realistic world with 20 tectonic plates
    let seed = 097243067;
    println!("\n⚡ Generating realistic tectonic plates using electrostatic physics...");
    let mut generator = TectonicPlateGenerator::with_seed(1800, 900, 20, seed)?;
    
    // Generate plates with extreme Earth-like size variety
    generator.generate("electrostatic", true)?;
    println!("✅ Plates generated successfully");
    
    // Export in all available formats
    generator.export_all("outputs", "tectonic_plates")?;
    println!("✅ Export complete");
    
    // Show statistics
    let stats = generator.get_plate_stats();
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
    
    println!("\n🎉 GENERATION COMPLETE!");
    println!("Files created in outputs/ directory:");
    println!("  • tectonic_plates.bin - Raw binary data");
    
    #[cfg(feature = "export-png")]
    println!("  • tectonic_plates.png - Visualization");
    
    #[cfg(feature = "export-tiff")]
    println!("  • tectonic_plates.tiff - GeoTIFF for GIS");
    
    Ok(())
}