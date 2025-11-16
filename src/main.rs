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
    println!("\n🗺️ Creating new world map (1800x900, seed: {})...", seed);
    let mut world = WorldMap::new(1800, 900, seed)?;

    // Run complete Stage 1 pipeline using new API
    println!("\n🌍 Running complete Stage 1: Tectonic Foundation pipeline...");
    println!("\n⚡ Stage 1.1: Generating tectonic plates...");
    world.tectonics().generate_plates(20)?;

    println!("🎨 Stage 1.2: Refining plate boundaries...");
    world.tectonics().roughen_boundaries(None)?;

    println!("🏝️  Stage 1.3: Removing plate islands...");
    let island_stats = world.tectonics().deisland(None)?;
    if island_stats.islands_removed > 0 {
        println!("   Removed {} islands ({} pixels reassigned)",
                 island_stats.islands_removed, island_stats.pixels_reassigned);
    }

    println!("🔍 Stage 1.4: Analyzing plate boundaries...");
    let boundary_stats = world.tectonics().analyze(None)?;
    println!("   Found {} plate boundaries:", boundary_stats.total_boundaries);
    println!("   • Convergent (colliding):    {}", boundary_stats.convergent_count);
    println!("   • Divergent (spreading):     {}", boundary_stats.divergent_count);
    println!("   • Transform (sliding):       {}", boundary_stats.transform_count);
    println!("   • Total length: {:.0} km", boundary_stats.total_length_km);
    println!("   • Avg relative velocity: {:.2} cm/year", boundary_stats.average_relative_velocity);

    // Show statistics
    if let Some(metadata) = world.get_tectonic_metadata() {
        println!("\n📊 Plate Statistics ({} plates total):", metadata.plate_stats.len());

        // Sort by size for better display
        let mut sorted_stats: Vec<_> = metadata.plate_stats.iter().collect();
        sorted_stats.sort_by(|a, b| b.1.area_km2.partial_cmp(&a.1.area_km2).unwrap());

        for (i, (plate_id, stat)) in sorted_stats.iter().enumerate().take(5) {
            let category = match i {
                0 => "Superplate",
                1..=2 => "Major",
                3..=4 => "Medium",
                _ => "Small",
            };
            let plate_type_icon = match stat.plate_type {
                geoforge::PlateType::Continental => "🏔️",
                geoforge::PlateType::Oceanic => "🌊",
                geoforge::PlateType::Mixed => "🏝️",
            };

            // Find motion info
            if let Some(seed) = metadata.plate_seeds.iter().find(|s| s.id == **plate_id) {
                println!("  {} {} {}: {:.1}% ({:.0} km²) - moving {:.0}° at {:.1} cm/yr",
                         plate_type_icon, category, plate_id, stat.percentage, stat.area_km2,
                         seed.motion_direction, seed.motion_speed);
            }
        }

        if sorted_stats.len() > 5 {
            println!("  ... and {} smaller plates", sorted_stats.len() - 5);
        }

        // Show plate type distribution
        let mut oceanic = 0;
        let mut continental = 0;
        let mut mixed = 0;
        for stat in metadata.plate_stats.values() {
            match stat.plate_type {
                geoforge::PlateType::Oceanic => oceanic += 1,
                geoforge::PlateType::Continental => continental += 1,
                geoforge::PlateType::Mixed => mixed += 1,
            }
        }
        println!("\n  Plate Types: {} continental, {} oceanic, {} mixed",
                 continental, oceanic, mixed);
    }

    // Stage 2: Geological Provinces (Orogenic Belts)
    println!("\n🏔️  Stage 2.1: Generating orogenic belts...");
    let orogens = world.generate_geology(None)?;

    // Show orogen statistics
    let mut collision = 0;
    let mut subduction = 0;
    let mut accretionary = 0;
    for orogen in &orogens {
        match orogen.characteristics.province_type {
            geoforge::GeologicProvince::CollisionOrogen => collision += 1,
            geoforge::GeologicProvince::SubductionOrogen => subduction += 1,
            geoforge::GeologicProvince::AccretionaryOrogen => accretionary += 1,
        }
    }
    println!("   Generated {} orogenic belts:", orogens.len());
    println!("   • Collision (continental-continental):  {}", collision);
    println!("   • Subduction (oceanic-continental):     {}", subduction);
    println!("   • Accretionary (mixed/terrane):         {}", accretionary);

    // Export all tectonic data using new API
    println!("\n💾 Exporting visualizations...");
    world.tectonics().export("outputs")?;

    #[cfg(feature = "export-png")]
    {
        println!("✅ Plate boundaries exported: outputs/tectonics.png");
        println!("✅ Boundary types exported: outputs/boundaries.png");
        println!("   (Red=convergent, Blue=divergent, Green=transform)");
        println!("✅ Plate motion exported: outputs/plate_motion.png");
        println!("   (Color=direction, Brightness=speed)");

        // Export geology
        world.export_geology_png("outputs", "geology.png")?;
        println!("✅ Geological provinces exported: outputs/geology.png");
        println!("   (Red=collision, Orange=subduction, Yellow=accretionary)");
    }

    println!("✅ Complete world data saved: outputs/world.map");

    println!("\n🎉 STAGES 1-2: TECTONIC & GEOLOGICAL FOUNDATION COMPLETE!");
    println!("\nPipeline executed:");
    println!("  ✅ Stage 1.1: Core Plate Generation (electrostatic physics)");
    println!("  ✅ Stage 1.2: Boundary Refinement (realistic irregularity)");
    println!("  ✅ Stage 1.3: Island Removal (contiguous plates)");
    println!("  ✅ Stage 1.4: Boundary Analysis (motion & classification)");
    println!("  ✅ Stage 2.1: Orogenic Belts (mountain-building zones)");
    println!("\nFiles created in outputs/ directory:");
    println!("  • world.map - Complete world data (binary)");

    #[cfg(feature = "export-png")]
    {
        println!("  • tectonics.png - Tectonic plates (color-coded)");
        println!("  • boundaries.png - Boundary types (red/blue/green)");
        println!("  • plate_motion.png - Motion vectors (hue=direction, sat=speed)");
        println!("  • geology.png - Orogenic belts (red/orange/yellow)");
        println!("\n📖 Motion Visualization Color Key:");
        println!("  • Red → Eastward    • Yellow → Northward");
        println!("  • Cyan → Westward   • Blue → Southward");
        println!("  • Brighter = faster, Grayer = slower");
        println!("\n🏔️  Geology Visualization Color Key:");
        println!("  • Red = Collision orogens (continental-continental)");
        println!("  • Orange = Subduction orogens (oceanic-continental)");
        println!("  • Yellow = Accretionary orogens (mixed/terrane)");
    }

    #[cfg(not(feature = "export-png"))]
    println!("\n💡 Run with --features export-png for visualization output");

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

    // Import the PNG using new API
    world.tectonics().import_png(png_path)?;
    
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