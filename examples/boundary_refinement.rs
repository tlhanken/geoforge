//! Example demonstrating Stage 1.2: Boundary Refinement
//!
//! This example generates tectonic plates and then applies boundary refinement
//! to add realistic irregularity to the plate edges. It creates two versions:
//! - One with the original Voronoi boundaries
//! - One with refined, more natural boundaries
//!
//! Run with: cargo run --example boundary_refinement --features export-png

use geoforge::{WorldMap, BoundaryRefinementConfig};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🌍 Geoforge - Boundary Refinement Example");
    println!("==========================================\n");

    let seed = 12345;
    let width = 1800;
    let height = 900;
    let num_plates = 15;

    // Create output directory
    std::fs::create_dir_all("outputs/examples/boundary_refinement")?;

    // ========================================
    // Part 1: Generate with standard boundaries
    // ========================================
    println!("📐 Part 1: Generating plates with standard Voronoi boundaries...");
    let mut world_standard = WorldMap::new(width, height, seed)?;
    world_standard.generate_tectonics(num_plates, true)?;

    #[cfg(feature = "export-png")]
    {
        world_standard.export_tectonics_png(
            "outputs/examples/boundary_refinement",
            "01_standard_boundaries.png"
        )?;
        println!("✅ Standard boundaries exported to: outputs/examples/boundary_refinement/01_standard_boundaries.png");
    }

    // ========================================
    // Part 2: Generate with mild refinement
    // ========================================
    println!("\n🎨 Part 2: Applying mild boundary refinement...");
    let mut world_mild = WorldMap::new(width, height, seed)?;
    world_mild.generate_tectonics(num_plates, true)?;

    let mild_config = BoundaryRefinementConfig::with_seed(seed)
        .with_noise(0.01, 8.0, 3)    // Larger scale, modest warping
        .with_smoothing(1);

    world_mild.refine_boundaries(Some(mild_config))?;

    #[cfg(feature = "export-png")]
    {
        world_mild.export_tectonics_png(
            "outputs/examples/boundary_refinement",
            "02_mild_refinement.png"
        )?;
        println!("✅ Mild refinement exported to: outputs/examples/boundary_refinement/02_mild_refinement.png");
    }

    // ========================================
    // Part 3: Generate with moderate refinement
    // ========================================
    println!("\n🎨 Part 3: Applying moderate boundary refinement...");
    let mut world_moderate = WorldMap::new(width, height, seed)?;
    world_moderate.generate_tectonics(num_plates, true)?;

    let moderate_config = BoundaryRefinementConfig::with_seed(seed)
        .with_noise(0.020, 80.0, 5)  // Medium features, solid warping (was strong)
        .with_smoothing(1);

    world_moderate.refine_boundaries(Some(moderate_config))?;

    #[cfg(feature = "export-png")]
    {
        world_moderate.export_tectonics_png(
            "outputs/examples/boundary_refinement",
            "03_moderate_refinement.png"
        )?;
        println!("✅ Moderate refinement exported to: outputs/examples/boundary_refinement/03_moderate_refinement.png");
    }

    // ========================================
    // Part 4: Generate with strong refinement
    // ========================================
    println!("\n🎨 Part 4: Applying strong boundary refinement...");
    let mut world_strong = WorldMap::new(width, height, seed)?;
    world_strong.generate_tectonics(num_plates, true)?;

    let strong_config = BoundaryRefinementConfig::with_seed(seed)
        .with_noise(0.025, 100.0, 5)  // Tighter features, strong warping (new intermediate)
        .with_smoothing(2);

    world_strong.refine_boundaries(Some(strong_config))?;

    #[cfg(feature = "export-png")]
    {
        world_strong.export_tectonics_png(
            "outputs/examples/boundary_refinement",
            "04_strong_refinement.png"
        )?;
        println!("✅ Strong refinement exported to: outputs/examples/boundary_refinement/04_strong_refinement.png");
    }

    // ========================================
    // Part 5: Generate with extreme refinement
    // ========================================
    println!("\n🎨 Part 5: Applying extreme boundary refinement...");
    let mut world_extreme = WorldMap::new(width, height, seed)?;
    world_extreme.generate_tectonics(num_plates, true)?;

    let extreme_config = BoundaryRefinementConfig::with_seed(seed)
        .with_noise(0.030, 120.0, 6)  // Tight, frequent features with massive warping
        .with_smoothing(3);

    world_extreme.refine_boundaries(Some(extreme_config))?;

    #[cfg(feature = "export-png")]
    {
        world_extreme.export_tectonics_png(
            "outputs/examples/boundary_refinement",
            "05_extreme_refinement.png"
        )?;
        println!("✅ Extreme refinement exported to: outputs/examples/boundary_refinement/05_extreme_refinement.png");
    }

    // ========================================
    // Summary
    // ========================================
    println!("\n🎉 BOUNDARY REFINEMENT EXAMPLE COMPLETE!");
    println!("\n📊 Generated {} variations:", 5);
    println!("  1. Standard Voronoi boundaries (baseline)");
    println!("  2. Mild refinement (subtle irregularity)");
    println!("  3. Moderate refinement (default, realistic)");
    println!("  4. Strong refinement (more pronounced features)");
    println!("  5. Extreme refinement (very irregular, artistic)");

    #[cfg(feature = "export-png")]
    println!("\n💡 Compare the PNG files in outputs/examples/boundary_refinement/");

    #[cfg(not(feature = "export-png"))]
    println!("\n💡 Run with --features export-png to generate visualization images");

    println!("\n📚 Key parameters:");
    println!("  • noise_scale: Controls feature size (0.01-1.0)");
    println!("  • noise_amplitude: How far boundaries shift (1-10 pixels)");
    println!("  • octaves: Number of noise layers (more = more detail)");
    println!("  • smoothing_iterations: Post-refinement smoothing (0-3)");

    Ok(())
}
