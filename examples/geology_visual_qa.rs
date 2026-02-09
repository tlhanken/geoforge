/// Visual QA for Stage 2 Geologic Provinces
///
/// This example generates test worlds to visually validate all 18 province types.
/// Exports PNGs for manual inspection of colors, widths, and geological accuracy.

use geoforge::{WorldMap, MapExporter};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🌋 Geologic Provinces Visual QA");
    println!("================================\n");

    std::fs::create_dir_all("outputs/geology_qa")?;

    // Test Case 1: General World - All Province Types
    println!("📍 Test Case 1: General World (seed 12345)");
    println!("   Purpose: Showcase variety of all province types");
    generate_test_world(
        12345,
        1800,
        900,
        20,
        "01_general_world",
        "General world with mix of continental and oceanic features"
    )?;

    // Test Case 2: Collision-Heavy World
    println!("\n📍 Test Case 2: Collision-Heavy World (seed 777)");
    println!("   Purpose: Validate collision orogens (mountain belts)");
    generate_test_world(
        777,
        1800,
        900,
        15,
        "02_collision_heavy",
        "World optimized for collision orogens (continent-continent boundaries)"
    )?;

    // Test Case 3: Arc Systems World
    println!("\n📍 Test Case 3: Subduction Arc World (seed 333)");
    println!("   Purpose: Validate trenches, arcs, wedges, and basins");
    generate_test_world(
        333,
        1800,
        900,
        18,
        "03_arc_systems",
        "World showcasing subduction zone features (trenches, arcs, basins)"
    )?;

    // Test Case 4: Large Continental Plates
    println!("\n📍 Test Case 4: Large Continents (seed 555)");
    println!("   Purpose: Validate shields, platforms, intracratonic basins");
    generate_test_world(
        555,
        1800,
        900,
        8,
        "04_large_continents",
        "Fewer, larger plates to showcase stable continental features"
    )?;

    // Test Case 5: Oceanic World
    println!("\n📍 Test Case 5: Oceanic-Dominated (seed 999)");
    println!("   Purpose: Validate mid-ocean ridges, fracture zones, hotspots");
    generate_test_world(
        999,
        1800,
        900,
        25,
        "05_oceanic_world",
        "Many small plates to maximize oceanic features"
    )?;

    // Test Case 6: High Latitude World
    println!("\n📍 Test Case 6: Polar Regions (seed 2468)");
    println!("   Purpose: Validate spherical projection at high latitudes");
    generate_test_world(
        2468,
        1800,
        900,
        15,
        "06_polar_regions",
        "Test spherical-aware width scaling near poles"
    )?;

    println!("\n✅ All test worlds generated successfully!");
    println!("\n📂 Output directory: outputs/geology_qa/");
    println!("\n🔍 Visual Inspection Checklist:");
    println!("   □ Province Colors:");
    println!("     • Collision orogens: Browns/whites (high elevation)");
    println!("     • Volcanic arcs: Reds/oranges (active volcanism)");
    println!("     • Ocean trenches: Deep blues (deepest ocean)");
    println!("     • Cratons: Oranges (ancient continental cores)");
    println!("     • Platforms: Pinks (stable continental)");
    println!("     • Abyssal plains: Blues (oceanic base)");
    println!("     • Mid-ocean ridges: Lighter blues (elevated seafloor)");
    println!("   □ Province Widths:");
    println!("     • Trenches: Narrow (~50 km)");
    println!("     • Volcanic arcs: Medium (~50-100 km)");
    println!("     • Accretionary wedges: ~100-200 km");
    println!("     • Forearc basins: ~100-200 km");
    println!("     • Backarc basins: ~200-400 km");
    println!("     • Mid-ocean ridges: ~60-150 km (varies by spreading rate)");
    println!("     • Collision orogens: Wide (~500-2000 km)");
    println!("   □ Province Positioning:");
    println!("     • Trenches: On subducting plate side");
    println!("     • Wedges/arcs/basins: On overriding plate side");
    println!("     • Ridges: At divergent boundaries");
    println!("     • Hotspots: In deep plate interiors (far from edges)");
    println!("     • Shields: Near plate centers");
    println!("   □ No Visual Artifacts:");
    println!("     • No gaps between provinces");
    println!("     • No unexpected overlaps");
    println!("     • No rendering glitches at map edges");
    println!("     • Consistent widths at similar latitudes");

    Ok(())
}

fn generate_test_world(
    seed: u64,
    width: usize,
    height: usize,
    num_plates: usize,
    filename_prefix: &str,
    description: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    // Create world
    let mut world = WorldMap::new(width, height, seed)?;

    // Complete Stage 1 pipeline (same as main production usage)
    // Stage 1.1: Generate plates
    world.generate_tectonics(num_plates)?;

    // Stage 1.2: Boundary refinement (roughening)
    world.refine_boundaries(None)?; // Uses default config

    // Stage 1.3: Island removal
    world.remove_islands(None)?;

    // Stage 1.4: Motion analysis and boundary classification
    world.analyze_boundaries(None)?;

    let metadata = world.get_tectonic_metadata().unwrap();
    println!("   Generated {} plates", metadata.plate_stats.len());

    // Count plate types
    let continental_count = metadata.plate_stats.values()
        .filter(|s| s.plate_type == geoforge::PlateType::Continental)
        .count();
    let oceanic_count = metadata.plate_stats.len() - continental_count;
    println!("   Continental: {}, Oceanic: {}", continental_count, oceanic_count);

    // Generate geology
    let provinces = world.generate_geology(None)?;
    println!("   Generated {} geologic provinces", provinces.len());

    // Count province types
    use std::collections::HashMap;
    let mut type_counts: HashMap<String, usize> = HashMap::new();
    for province in &provinces {
        let type_name = format!("{:?}", province.characteristics.province_type);
        *type_counts.entry(type_name).or_insert(0) += 1;
    }

    // Show province distribution
    let mut sorted_types: Vec<_> = type_counts.iter().collect();
    sorted_types.sort_by_key(|(_, count)| std::cmp::Reverse(**count));

    print!("   Province types: ");
    for (i, (type_name, count)) in sorted_types.iter().enumerate() {
        if i > 0 { print!(", "); }
        print!("{}: {}", type_name, count);
    }
    println!();

    #[cfg(feature = "export-png")]
    {
        // Export tectonics visualization
        world.export_tectonics_png(
            "outputs/geology_qa",
            &format!("{}_tectonics.png", filename_prefix)
        )?;

        // Export geology visualization
        world.export_geology_png(
            "outputs/geology_qa",
            &format!("{}_geology.png", filename_prefix)
        )?;

        println!("   Exported: outputs/geology_qa/{}_*.png", filename_prefix);
        println!("   Description: {}", description);
    }

    #[cfg(not(feature = "export-png"))]
    {
        println!("   ⚠️  PNG export disabled - run with --features export-png");
        let _ = description; // Suppress unused warning
    }

    Ok(())
}
