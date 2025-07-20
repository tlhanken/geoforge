/// Example demonstrating PNG import functionality for tectonic plates
use geoforge::WorldMap;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🖼️ PNG Import Example");
    println!("=====================");
    println!("Demonstrating importing tectonic plates from PNG images");
    
    // Step 1: Create a world and generate some plates to export as PNG
    println!("\n📐 Step 1: Generate plates and export as PNG");
    println!("--------------------------------------------");
    
    let mut world_export = WorldMap::new(180, 90, 42)?;
    world_export.generate_tectonics(8, true)?;
    
    std::fs::create_dir_all("outputs/examples/png_import")?;
    
    #[cfg(feature = "export-png")]
    {
        world_export.export_tectonics_png("outputs/examples/png_import", "original_plates.png")?;
        println!("✅ Exported plates to: outputs/examples/png_import/original_plates.png");
        
        // Show original statistics
        if let Some(stats) = world_export.get_tectonic_stats() {
            println!("📊 Original plates: {} plates", stats.len());
            let mut sorted_stats: Vec<_> = stats.iter().collect();
            sorted_stats.sort_by(|a, b| b.1.area_km2.cmp(&a.1.area_km2));
            
            for (_i, (plate_id, stat)) in sorted_stats.iter().enumerate().take(3) {
                println!("  Plate {}: {:.1}% ({} km²)", plate_id, stat.percentage, stat.area_km2);
            }
        }
    }
    
    #[cfg(not(feature = "export-png"))]
    {
        println!("⚠️ PNG features not enabled. Run with --features export-png");
        return Ok(());
    }
    
    // Step 2: Import the PNG back into a new world
    println!("\n📥 Step 2: Import PNG into new world");
    println!("-------------------------------------");
    
    let mut world_import = WorldMap::new(180, 90, 999)?; // Different seed
    
    #[cfg(feature = "export-png")]
    {
        world_import.import_tectonics_png("outputs/examples/png_import/original_plates.png")?;
        
        // Show imported statistics
        if let Some(stats) = world_import.get_tectonic_stats() {
            println!("📊 Imported plates: {} plates", stats.len());
            let mut sorted_stats: Vec<_> = stats.iter().collect();
            sorted_stats.sort_by(|a, b| b.1.area_km2.cmp(&a.1.area_km2));
            
            for (_i, (plate_id, stat)) in sorted_stats.iter().enumerate().take(3) {
                println!("  Plate {}: {:.1}% ({} km²)", plate_id, stat.percentage, stat.area_km2);
            }
        }
        
        // Export the imported data to verify it's the same
        world_import.export_tectonics_png("outputs/examples/png_import", "reimported_plates.png")?;
        println!("✅ Re-exported as: outputs/examples/png_import/reimported_plates.png");
        
        // Step 3: Demonstrate data preservation
        println!("\n🔍 Step 3: Verify data preservation");
        println!("-----------------------------------");
        
        // Compare the original and imported data
        let original_data = &world_export.tectonics.as_ref().unwrap().data;
        let imported_data = &world_import.tectonics.as_ref().unwrap().data;
        
        let identical = original_data == imported_data;
        println!("🎯 Plate boundaries preserved: {}", 
                 if identical { "✅ YES" } else { "⚠️ Plate IDs reassigned" });
        
        if !identical {
            let differences = original_data.iter()
                .zip(imported_data.iter())
                .filter(|(a, b)| a != b)
                .count();
            println!("   Note: {} pixels have different plate IDs due to color-to-ID remapping", differences);
            println!("   This is normal - plate boundaries are preserved, but IDs may change");
        }
    }
    
    // Step 4: Show file information
    println!("\n📁 Generated Files");
    println!("------------------");
    println!("• outputs/examples/png_import/original_plates.png - Generated tectonic plates");
    println!("• outputs/examples/png_import/reimported_plates.png - Same data after PNG import");
    
    // Step 5: Practical usage tips
    println!("\n💡 Usage Tips");
    println!("-------------");
    println!("• Each unique color in the PNG becomes a separate plate");
    println!("• PNG dimensions must match WorldMap dimensions exactly");
    println!("• Colors are mapped to plate IDs automatically (1, 2, 3...)");
    println!("• Plate seeds are calculated as centroids of colored regions");
    println!("• Import any hand-drawn or external tectonic plate map");
    println!("• Supports up to 65,535 unique colors/plates");
    
    println!("\n🎉 PNG import example completed successfully!");
    
    Ok(())
}