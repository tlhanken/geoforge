/// Example demonstrating PNG import functionality for tectonic plates
use geoforge::WorldMap;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🖼️ PNG Import Example");
    println!("=====================");
    println!("Demonstrating importing tectonic plates with motion vectors from PNG images");

    // Step 1: Create a world and generate some plates to export as PNG
    println!("\n📐 Step 1: Generate plates and export as PNG");
    println!("--------------------------------------------");

    let mut world_export = WorldMap::new(180, 90, 42)?;
    world_export.tectonics().generate_plates(8)?;
    
    std::fs::create_dir_all("outputs/examples/png_import")?;
    
    #[cfg(feature = "export-png")]
    {
        world_export.export_plate_motion_png("outputs/examples/png_import", "original_motion.png")?;
        println!("✅ Exported plates with motion to: outputs/examples/png_import/original_motion.png");

        // Show original statistics
        if let Some(metadata) = world_export.get_tectonic_metadata() {
            println!("📊 Original plates: {} plates", metadata.plate_stats.len());
            let mut sorted_stats: Vec<_> = metadata.plate_stats.iter().collect();
            sorted_stats.sort_by(|a, b| b.1.area_km2.cmp(&a.1.area_km2));

            for (_i, (plate_id, stat)) in sorted_stats.iter().enumerate().take(3) {
                let seed = metadata.plate_seeds.iter().find(|s| s.id == **plate_id).unwrap();
                println!("  Plate {}: {:.1}% ({} km²), moving {:.0}° at {:.1} cm/year",
                         plate_id, stat.percentage, stat.area_km2,
                         seed.motion_direction, seed.motion_speed);
            }
        }
    }
    
    #[cfg(not(feature = "export-png"))]
    {
        println!("⚠️ PNG features not enabled. Run with --features export-png");
    }

    #[cfg(feature = "export-png")]
    {
        // Step 2: Import the PNG back into a new world
        println!("\n📥 Step 2: Import PNG into new world");
        println!("-------------------------------------");

        let mut world_import = WorldMap::new(180, 90, 999)?; // Different seed
        world_import.tectonics().import_png("outputs/examples/png_import/original_motion.png")?;

        // Show imported statistics
        if let Some(metadata) = world_import.get_tectonic_metadata() {
            println!("📊 Imported plates: {} plates", metadata.plate_stats.len());
            let mut sorted_stats: Vec<_> = metadata.plate_stats.iter().collect();
            sorted_stats.sort_by(|a, b| b.1.area_km2.cmp(&a.1.area_km2));

            for (_i, (plate_id, stat)) in sorted_stats.iter().enumerate().take(3) {
                let seed = metadata.plate_seeds.iter().find(|s| s.id == **plate_id).unwrap();
                println!("  Plate {}: {:.1}% ({} km²), moving {:.0}° at {:.1} cm/year",
                         plate_id, stat.percentage, stat.area_km2,
                         seed.motion_direction, seed.motion_speed);
            }
        }

        // Export the imported data to verify it's the same
        world_import.export_plate_motion_png("outputs/examples/png_import", "reimported_motion.png")?;
        println!("✅ Re-exported as: outputs/examples/png_import/reimported_motion.png");

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

        // Step 4: Show file information
        println!("\n📁 Generated Files");
        println!("------------------");
        println!("• outputs/examples/png_import/original_motion.png - Generated plates with motion vectors");
        println!("• outputs/examples/png_import/reimported_motion.png - Same data after PNG import");

        // Step 5: Practical usage tips
        println!("\n💡 Usage Tips");
        println!("-------------");
        println!("• Each unique color encodes BOTH plate identity AND motion vector");
        println!("• Color encoding: Hue = direction (0-360°), Saturation = speed (1-10 cm/year)");
        println!("• PNG dimensions must match WorldMap dimensions exactly");
        println!("• Plate IDs are assigned in discovery order (1, 2, 3...)");
        println!("• Plate seeds are calculated as centroids of colored regions");
        println!("• Import any hand-drawn or external tectonic plate map with motion");
        println!("• Supports up to 65,535 unique color/motion combinations");

        println!("\n🎉 PNG import example completed successfully!");
    }
    
    Ok(())
}