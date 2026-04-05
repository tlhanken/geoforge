/// Example demonstrating the generate_all function
use geoforge::{MapExporter, WorldMap};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🌍 Generate All Layers Example");
    println!("==============================");

    // Create a world
    let mut world = WorldMap::new(180, 90, 42)?;

    // Generate all available layers at once
    world.generate_all(8)?;

    // Show what was generated
    println!("\n📊 Generated Layers:");
    if world.tectonics.is_some() {
        println!("✅ Tectonics layer generated");
        if let Some(stats) = world.get_tectonic_stats() {
            println!("   {} plates created", stats.len());
        }
    }

    if world.elevation.is_some() {
        println!("✅ Elevation layer generated");
    } else {
        println!("⏳ Elevation layer - not yet implemented");
    }

    if world.temperature.is_some() {
        println!("✅ Temperature layer generated");
    } else {
        println!("⏳ Temperature layer - not yet implemented");
    }

    if world.precipitation.is_some() {
        println!("✅ Precipitation layer generated");
    } else {
        println!("⏳ Precipitation layer - not yet implemented");
    }

    if world.biomes.is_some() {
        println!("✅ Biomes layer generated");
    } else {
        println!("⏳ Biomes layer - not yet implemented");
    }

    // Export all generated layers
    #[cfg(feature = "export-png")]
    {
        println!("\n📤 Exporting all layers...");
        world.export_all_png("outputs/examples/generate_all", "complete_world")?;
    }

    #[cfg(not(feature = "export-png"))]
    {
        println!("\n⚠️ PNG export features not enabled. Run with --features export-png");
    }

    // Save the complete world
    world.save_to_file("outputs/examples/generate_all/complete_world.map")?;
    println!("💾 Saved complete world to: outputs/examples/generate_all/complete_world.map");

    println!("\n🎉 Generate all example completed!");
    println!("📁 Check outputs/examples/generate_all/ directory for output files");

    Ok(())
}
