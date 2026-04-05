/// Enhanced example showcasing the geoforge library capabilities
use geoforge::{GenerationMethod, MapExporter, TectonicPlateGenerator, WorldMap};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🌍 Geoforge Library Example");
    println!("============================");
    println!("Demonstrating realistic geological world generation");

    // Example 1: Direct tectonic plate generation (low-level API)
    println!("\n📐 Example 1: Direct Tectonic Generation");
    println!("----------------------------------------");

    let mut generator = TectonicPlateGenerator::with_seed(360, 180, 8, 42)?;
    println!("• Creating 8 plates on 360×180 grid (1° resolution)");
    println!("• Using electrostatic physics simulation");
    println!("• Seed: {}", generator.get_seed());

    generator.generate(GenerationMethod::Electrostatic)?;
    generator.validate()?;

    // Show plate statistics
    let stats = generator.get_plate_stats();
    println!(
        "\n📊 Generated {} plates with Earth-like size distribution:",
        stats.len()
    );

    let mut sorted_stats: Vec<_> = stats.into_iter().collect();
    sorted_stats.sort_by(|a, b| b.1.area_km2.cmp(&a.1.area_km2));

    for (i, (plate_id, stat)) in sorted_stats.iter().enumerate().take(5) {
        let category = match i {
            0 => "🌍 Superplate",
            1 => "🏔️  Major plate",
            2 => "⛰️  Large plate",
            _ => "🗻 Medium plate",
        };
        println!(
            "  {} {}: {:.1}% ({:>8} km²)",
            category, plate_id, stat.percentage, stat.area_km2
        );
    }

    if sorted_stats.len() > 5 {
        println!("  ... and {} smaller plates", sorted_stats.len() - 5);
    }

    // Example 2: High-level WorldMap API (recommended)
    println!("\n🗺️ Example 2: WorldMap API (Recommended)");
    println!("------------------------------------------");

    let mut world = WorldMap::new(180, 90, 123456)?;
    println!("• Creating world map (180×90, 2° resolution)");
    println!("• Seed: {}", world.seed);

    // Generate tectonic layer
    world.generate_tectonics(12)?;
    println!("• Generated 12 tectonic plates");

    // Show world statistics
    if let Some(tectonic_stats) = world.get_tectonic_stats() {
        let total_plates = tectonic_stats.len();
        let areas: Vec<u64> = tectonic_stats.values().map(|s| s.area_km2).collect();
        let largest = areas.iter().max().unwrap_or(&0);
        let smallest = areas.iter().min().unwrap_or(&0);
        let ratio = if *smallest > 0 {
            *largest as f64 / *smallest as f64
        } else {
            0.0
        };

        println!(
            "  📈 Plate size variety: {:.1}x ratio (largest/smallest)",
            ratio
        );
        println!("  🌐 Total plates: {}", total_plates);
        println!(
            "  📏 World dimensions: {}×{} pixels",
            world.width, world.height
        );
    }

    // Example 3: Demonstrate deterministic generation
    println!("\n🔄 Example 3: Deterministic Generation");
    println!("--------------------------------------");

    let seed = 999;
    let mut world1 = WorldMap::new(90, 45, seed)?;
    let mut world2 = WorldMap::new(90, 45, seed)?;

    world1.generate_tectonics(6)?;
    world2.generate_tectonics(6)?;

    // Compare the results
    let data1 = &world1.tectonics.as_ref().unwrap().data;
    let data2 = &world2.tectonics.as_ref().unwrap().data;

    let identical = data1 == data2;
    println!(
        "• Same seed ({}) produces identical results: {}",
        seed,
        if identical { "✅ YES" } else { "❌ NO" }
    );

    // Example 4: Show memory efficiency
    println!("\n💾 Example 4: Memory Efficiency");
    println!("-------------------------------");

    let (_width, _height, plate_data, seeds) = generator.get_plate_data();
    let data_size_kb = plate_data.len() * std::mem::size_of::<u16>() / 1024;
    let seeds_size_bytes = seeds.len() * std::mem::size_of::<geoforge::PlateSeed>();

    println!(
        "• Plate data: {:.1} KB ({} pixels)",
        data_size_kb,
        plate_data.len()
    );
    println!(
        "• Seed data: {} bytes ({} seeds)",
        seeds_size_bytes,
        seeds.len()
    );
    println!(
        "• Total memory: {:.1} KB",
        data_size_kb as f64 + seeds_size_bytes as f64 / 1024.0
    );

    // Example 5: Export capabilities
    println!("\n📤 Example 5: Export Capabilities");
    println!("---------------------------------");

    std::fs::create_dir_all("outputs/examples/simple")?;

    // Export binary data
    generator.export_binary("outputs/examples/simple", "plates.bin")?;
    println!("• Binary export: outputs/examples/simple/plates.bin");

    // Save complete world
    world.save_to_file("outputs/examples/simple/world.map")?;
    println!("• World save: outputs/examples/simple/world.map");

    #[cfg(feature = "export-png")]
    {
        world.export_tectonics_png("outputs/examples/simple", "tectonics.png")?;
        println!("• PNG visualization: outputs/examples/simple/tectonics.png");
    }

    #[cfg(not(feature = "export-png"))]
    println!("• PNG export available with --features export-png");

    println!("\n🎉 All examples completed successfully!");
    println!("\n💡 Next steps:");
    println!("  • Run with --features export-png for visualizations");
    println!("  • Check outputs/examples/simple/ directory for generated files");
    println!("  • Try different seeds and parameters");
    println!("  • Use WorldMap API for multi-layer world generation");

    Ok(())
}
