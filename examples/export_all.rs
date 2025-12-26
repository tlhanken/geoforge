/// Example demonstrating the export_all_png function
use geoforge::{WorldMap, MapExporter};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("📤 Export All PNG Example");
    println!("=========================");
    
    // Create and generate a world
    let mut world = WorldMap::new(180, 90, 42)?;
    world.generate_tectonics(8)?;
    
    println!("✅ Generated world with tectonics");
    
    #[cfg(feature = "export-png")]
    {
        // Export all available layers
        world.export_all_png("outputs/examples/export_all", "world")?;
        
        println!("\n📁 Check outputs/examples/export_all/ directory for:");
        println!("  • world_tectonics.png - Tectonic plates");
        println!("  • Future layers will appear here as they're implemented");
    }
    
    #[cfg(not(feature = "export-png"))]
    {
        println!("⚠️ PNG features not enabled. Run with --features export-png");
    }
    
    println!("\n🎉 Export all example completed!");
    
    Ok(())
}