//! Example demonstrating Stage 1.3: Island Removal / Plate Contiguity
//!
//! This example applies the same boundary refinement levels as the boundary_refinement
//! example, but adds island removal to ensure all plates are contiguous. Each refinement
//! level shows the result after islands have been cleaned up.
//!
//! Compare with boundary_refinement example to see the before/after of island removal.
//!
//! Run with: cargo run --example island_removal --features export-png

use geoforge::{WorldMap, BoundaryRefinementConfig, MapExporter};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🌍 Geoforge - Island Removal Example");
    println!("=====================================\n");
    println!("Same refinement levels as boundary_refinement example,");
    println!("but with island removal applied for contiguous plates.\n");

    let seed = 12345;  // Same seed as boundary_refinement example
    let width = 1800;
    let height = 900;
    let num_plates = 15;

    // Create output directory
    std::fs::create_dir_all("outputs/examples/island_removal")?;

    // ========================================
    // Part 1: Mild refinement + island removal
    // ========================================
    println!("🎨 Part 1: Mild refinement with island removal...");
    let mut world_mild = WorldMap::new(width, height, seed)?;
    world_mild.generate_tectonics(num_plates)?;

    let mild_config = BoundaryRefinementConfig::with_seed(seed)
        .with_noise(0.01, 8.0, 3)
        .with_smoothing(1);

    world_mild.refine_boundaries(Some(mild_config))?;
    let mild_stats = world_mild.remove_islands(None)?;

    println!("  Islands removed: {}", mild_stats.islands_removed);
    println!("  Pixels reassigned: {}", mild_stats.pixels_reassigned);

    #[cfg(feature = "export-png")]
    {
        world_mild.export_tectonics_png(
            "outputs/examples/island_removal",
            "01_mild_refinement_contiguous.png"
        )?;
        println!("✅ Exported: outputs/examples/island_removal/01_mild_refinement_contiguous.png");
    }

    // ========================================
    // Part 2: Moderate refinement + island removal
    // ========================================
    println!("\n🎨 Part 2: Moderate refinement with island removal...");
    let mut world_moderate = WorldMap::new(width, height, seed)?;
    world_moderate.generate_tectonics(num_plates)?;

    let moderate_config = BoundaryRefinementConfig::with_seed(seed)
        .with_noise(0.020, 80.0, 5)
        .with_smoothing(1);

    world_moderate.refine_boundaries(Some(moderate_config))?;
    let moderate_stats = world_moderate.remove_islands(None)?;

    println!("  Islands removed: {}", moderate_stats.islands_removed);
    println!("  Pixels reassigned: {}", moderate_stats.pixels_reassigned);

    #[cfg(feature = "export-png")]
    {
        world_moderate.export_tectonics_png(
            "outputs/examples/island_removal",
            "02_moderate_refinement_contiguous.png"
        )?;
        println!("✅ Exported: outputs/examples/island_removal/02_moderate_refinement_contiguous.png");
    }

    // ========================================
    // Part 3: Strong refinement + island removal
    // ========================================
    println!("\n🎨 Part 3: Strong refinement with island removal...");
    let mut world_strong = WorldMap::new(width, height, seed)?;
    world_strong.generate_tectonics(num_plates)?;

    let strong_config = BoundaryRefinementConfig::with_seed(seed)
        .with_noise(0.025, 100.0, 5)
        .with_smoothing(2);

    world_strong.refine_boundaries(Some(strong_config))?;
    let strong_stats = world_strong.remove_islands(None)?;

    println!("  Islands removed: {}", strong_stats.islands_removed);
    println!("  Pixels reassigned: {}", strong_stats.pixels_reassigned);

    #[cfg(feature = "export-png")]
    {
        world_strong.export_tectonics_png(
            "outputs/examples/island_removal",
            "03_strong_refinement_contiguous.png"
        )?;
        println!("✅ Exported: outputs/examples/island_removal/03_strong_refinement_contiguous.png");
    }

    // ========================================
    // Part 4: Extreme refinement + island removal
    // ========================================
    println!("\n🎨 Part 4: Extreme refinement with island removal...");
    let mut world_extreme = WorldMap::new(width, height, seed)?;
    world_extreme.generate_tectonics(num_plates)?;

    let extreme_config = BoundaryRefinementConfig::with_seed(seed)
        .with_noise(0.030, 120.0, 6)
        .with_smoothing(3);

    world_extreme.refine_boundaries(Some(extreme_config))?;
    let extreme_stats = world_extreme.remove_islands(None)?;

    println!("  Islands removed: {}", extreme_stats.islands_removed);
    println!("  Pixels reassigned: {}", extreme_stats.pixels_reassigned);

    #[cfg(feature = "export-png")]
    {
        world_extreme.export_tectonics_png(
            "outputs/examples/island_removal",
            "04_extreme_refinement_contiguous.png"
        )?;
        println!("✅ Exported: outputs/examples/island_removal/04_extreme_refinement_contiguous.png");
    }

    // ========================================
    // Summary
    // ========================================
    println!("\n🎉 ISLAND REMOVAL EXAMPLE COMPLETE!");
    println!("\n📊 Island removal statistics by refinement level:");
    println!("  Mild:     {} islands removed, {} pixels reassigned",
        mild_stats.islands_removed, mild_stats.pixels_reassigned);
    println!("  Moderate: {} islands removed, {} pixels reassigned",
        moderate_stats.islands_removed, moderate_stats.pixels_reassigned);
    println!("  Strong:   {} islands removed, {} pixels reassigned",
        strong_stats.islands_removed, strong_stats.pixels_reassigned);
    println!("  Extreme:  {} islands removed, {} pixels reassigned",
        extreme_stats.islands_removed, extreme_stats.pixels_reassigned);

    #[cfg(feature = "export-png")]
    println!("\n💡 Compare with outputs/examples/boundary_refinement/ to see the difference!");

    #[cfg(not(feature = "export-png"))]
    println!("\n💡 Run with --features export-png to generate visualization images");

    println!("\n📚 Key insights:");
    println!("  • More aggressive refinement creates more islands");
    println!("  • Island removal maintains jagged boundaries while ensuring contiguity");
    println!("  • Each plate becomes a single contiguous region");
    println!("  • Compare with boundary_refinement outputs to see islands removed");

    Ok(())
}
