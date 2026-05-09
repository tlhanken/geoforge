use clap::Parser;
#[cfg(feature = "export-png")]
use geoforge::MapExporter;
use geoforge::WorldMap;
use std::time::Instant;

/// Realistic geological and climate modeling for procedural world generation
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Cli {
    /// Import an existing PNG file to create tectonic plates
    #[arg(long, value_name = "FILE")]
    pub import_png: Option<String>,

    /// Provide a specific seed for generation
    #[arg(long, short)]
    pub seed: Option<u64>,

    /// Enable full export including extra visualizations
    #[arg(long)]
    pub export_full: bool,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    println!("🌍 Geoforge - Realistic World Generation");
    println!("=====================================");

    if let Some(png_path) = &cli.import_png {
        return import_png_mode(png_path);
    }

    // Seed references:
    // General Test, good set of oceanic interactions: 097243067,
    // Monocontinent with all boundary interactions: 2837,
    // Good continent Continent collision: 3487130930717999446

    let seed = cli.seed.unwrap_or_else(rand::random::<u64>);
    println!("World Seed: {}", seed);

    println!("\n🗺️ Creating new world map (1800x900, seed: {})...", seed);
    let mut world = WorldMap::new(1800, 900, seed)?;

    let start_time = Instant::now();

    // Run complete Stage 1 pipeline using new API
    println!("\n🌍 Running complete Stage 1: Tectonic Foundation pipeline...");
    println!("\n⚡ Stage 1.1: Generating tectonic plates...");
    let t1_start = Instant::now();
    world.tectonics().generate_plates(20)?;
    let t1_duration = t1_start.elapsed();

    println!("🎨 Stage 1.2: Refining plate boundaries...");
    let t2_start = Instant::now();
    world.tectonics().roughen_boundaries(None)?;
    let t2_duration = t2_start.elapsed();

    println!("🏝️  Stage 1.3: Removing plate islands...");
    let t3_start = Instant::now();
    let island_stats = world.tectonics().deisland(None)?;
    let t3_duration = t3_start.elapsed();
    if island_stats.islands_removed > 0 {
        println!(
            "   Removed {} islands ({} pixels reassigned)",
            island_stats.islands_removed, island_stats.pixels_reassigned
        );
    }

    println!("🔍 Stage 1.4: Analyzing plate boundaries...");
    let t4_start = Instant::now();
    let boundary_stats = world.tectonics().analyze(None)?;
    let t4_duration = t4_start.elapsed();
    println!(
        "   Found {} plate boundaries:",
        boundary_stats.total_boundaries
    );
    println!(
        "   • Convergent (colliding):    {}",
        boundary_stats.convergent_count
    );
    println!(
        "   • Divergent (spreading):     {}",
        boundary_stats.divergent_count
    );
    println!(
        "   • Transform (sliding):       {}",
        boundary_stats.transform_count
    );
    println!(
        "   • Total length: {:.0} km",
        boundary_stats.total_length_km
    );
    println!(
        "   • Avg relative velocity: {:.2} cm/year",
        boundary_stats.average_relative_velocity
    );

    // Show statistics
    if let Some(metadata) = world.get_tectonic_metadata() {
        println!(
            "\n📊 Plate Statistics ({} plates total):",
            metadata.plate_stats.len()
        );

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
            };

            // Find motion info
            if let Some(seed) = metadata.plate_seeds.iter().find(|s| s.id == **plate_id) {
                println!(
                    "  {} {} {}: {:.1}% ({:.0} km²) - moving {:.0}° at {:.1} cm/yr",
                    plate_type_icon,
                    category,
                    plate_id,
                    stat.percentage,
                    stat.area_km2,
                    seed.motion_direction,
                    seed.motion_speed
                );
            }
        }

        if sorted_stats.len() > 5 {
            println!("  ... and {} smaller plates", sorted_stats.len() - 5);
        }

        // Show plate type distribution
        let mut oceanic = 0;
        let mut continental = 0;
        for stat in metadata.plate_stats.values() {
            match stat.plate_type {
                geoforge::PlateType::Oceanic => oceanic += 1,
                geoforge::PlateType::Continental => continental += 1,
            }
        }
        println!(
            "\n  Plate Types: {} continental, {} oceanic",
            continental, oceanic
        );
    }

    // Stage 2: Geological Provinces (Orogenic Belts)
    println!("\n🏔️  Stage 2.1: Generating orogenic belts...");
    let t5_start = Instant::now();
    let orogens = world.generate_geology(None)?;
    let t5_duration = t5_start.elapsed();

    // Show geology statistics
    let mut counts = std::collections::HashMap::new();
    for region in &orogens {
        *counts
            .entry(region.characteristics.province_type)
            .or_insert(0) += 1;
    }

    println!("   Generated {} geological provinces:", orogens.len());
    for (province_type, count) in counts.iter() {
        println!("   • {}: {}", province_type.name(), count);
    }

    // Export all tectonic data using new API
    println!("\n💾 Exporting visualizations...");
    let export_start = Instant::now();
    let export_full = cli.export_full;
    world.tectonics().export("outputs", export_full)?;

    #[cfg(feature = "export-png")]
    if export_full {
        println!("✅ Plate boundaries exported: outputs/tectonics.png");
        println!("✅ Boundary types exported: outputs/tectonics_boundaries.png");
        println!("   (Red=convergent, Blue=divergent, Green=transform)");
        println!("✅ Plate motion exported: outputs/tectonics_motion.png");
        println!("   (Color=direction, Brightness=speed)");

        // Export motion color reference
        world.export_motion_reference_png("outputs", "motion_reference.png")?;
        println!("✅ Motion color reference: outputs/motion_reference.png");
        println!("   (Color wheel showing direction → color mapping)");

        // Export plate types
        world.export_plate_types_png("outputs", "tectonics_types.png")?;
        println!("✅ Plate types exported: outputs/tectonics_types.png");
        println!("   (COOL Blue/Cyan=oceanic, WARM Red/Orange=continental)");

        // Export geology
        world.export_geology_png("outputs", "geology.png")?;
        println!("✅ Geological provinces exported: outputs/geology.png");

        // Export geology with boundary overlays
        world.export_geology_with_boundaries_png("outputs", "geology_boundaries.png")?;
        println!("✅ Geology+boundaries exported: outputs/geology_boundaries.png");
        println!("   (Provinces in color + Red/Blue/Green boundary overlays)");
    }
    let export_duration = export_start.elapsed();
    let total_duration = start_time.elapsed();

    println!("✅ Complete world data saved: outputs/world.map");

    println!("\n🎉 STAGES 1-2: TECTONIC & GEOLOGICAL FOUNDATION COMPLETE!");
    println!("\nPipeline executed in {:.2?}:", total_duration);
    println!(
        "  ✅ Stage 1.1: Core Plate Generation   ({:.2?})",
        t1_duration
    );
    println!(
        "  ✅ Stage 1.2: Boundary Refinement     ({:.2?})",
        t2_duration
    );
    println!(
        "  ✅ Stage 1.3: Island Removal          ({:.2?})",
        t3_duration
    );
    println!(
        "  ✅ Stage 1.4: Boundary Analysis       ({:.2?})",
        t4_duration
    );
    println!(
        "  ✅ Stage 2.1: Orogenic Belts          ({:.2?})",
        t5_duration
    );
    println!(
        "  ✅ Exports                            ({:.2?})",
        export_duration
    );
    println!("\nFiles created in outputs/ directory:");
    println!("  • world.map - Complete world data (binary)");

    #[cfg(feature = "export-png")]
    {
        println!("  • tectonics.png - Tectonic plates (color-coded)");
        println!("  • tectonics_boundaries.png - Boundary types (red/blue/green)");
        println!("  • tectonics_motion.png - Motion vectors (hue=direction, sat=speed)");
        println!("  • tectonics_types.png - Plate character (oceanic vs continental)");
        println!("  • geology.png - Geological provinces (various colors)");
        println!("  • geology_boundaries.png - Provinces + boundary overlays");
        println!("\n🌊 Plate Types Color Key:");
        println!("  • COOL Blue/Cyan = Oceanic plates (denser crust)");
        println!("  • WARM Red/Orange = Continental plates (lighter crust)");
        println!("  • Purple/Magenta = Mixed plates");
        println!("\n📖 Motion Visualization Color Key:");
        println!("  • Red → Eastward    • Yellow → Northward");
        println!("  • Cyan → Westward   • Blue → Southward");
        println!("  • Brighter = faster, Grayer = slower");
        println!("\n🏔️  Geology Visualization Color Key:");
        println!("  • Red = Collision orogens (continental-continental)");
        println!("  • Orange = Subduction orogens (oceanic-continental)");
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
    use image::GenericImageView;
    use image::io::Reader as ImageReader;
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
            println!(
                "  {} {}: {:.1}% ({} km²)",
                category, plate_id, stat.percentage, stat.area_km2
            );
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn verify_cli() {
        use clap::CommandFactory;
        Cli::command().debug_assert();
    }
}
