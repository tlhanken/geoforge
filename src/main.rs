use geoforge::{TectonicPlateGenerator, PlateError};
use std::fs::File;
use std::io::Write;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🌍 Geoforge - Tectonic Plate Generation");
    println!("=====================================\n");

    // Example 1: Quick generation with random seed
    println!("🎲 Example 1: Random Generation");
    let mut generator_random = TectonicPlateGenerator::new(900, 450, 8)?;
    println!("Random seed: {}", generator_random.get_seed());
    generator_random.generate("region_growing", true)?;

    let stats = generator_random.get_plate_stats();
    println!("Generated {} plates", stats.len());
    for (plate_id, stat) in stats.iter().take(3) {
        println!("  Plate {}: {:.1}% of surface", plate_id, stat.percentage);
    }
    println!();

    // Example 2: Reproducible generation
    println!("🔄 Example 2: Reproducible Generation");
    let mut generator_fixed = TectonicPlateGenerator::with_seed(900, 450, 8, 42)?;
    println!("Fixed seed: {}", generator_fixed.get_seed());
    generator_fixed.generate("region_growing", true)?;
    generator_fixed.validate()?;

    let stats = generator_fixed.get_plate_stats();
    println!("Generated {} plates (reproducible)", stats.len());
    println!();

    // Example 3: Compare different methods
    println!("⚖️  Example 3: Comparing Methods");
    let base_seed = 123456;
    
    let mut voronoi_gen = TectonicPlateGenerator::with_seed(450, 225, 6, base_seed)?;
    let start = std::time::Instant::now();
    voronoi_gen.generate("voronoi", false)?;
    let voronoi_time = start.elapsed();
    
    let mut region_gen = TectonicPlateGenerator::with_seed(450, 225, 6, base_seed)?;
    let start = std::time::Instant::now();
    region_gen.generate("region_growing", true)?;
    let region_time = start.elapsed();

    println!("Voronoi method: {:.2?}", voronoi_time);
    println!("Region growing: {:.2?}", region_time);
    println!();

    // Example 4: Detailed world analysis
    println!("🔍 Example 4: Detailed World Analysis");
    let mut detailed_gen = TectonicPlateGenerator::with_seed(1800, 900, 15, 999)?;
    detailed_gen.generate("region_growing", true)?;
    detailed_gen.validate()?;

    let stats = detailed_gen.get_plate_stats();
    let (width, height, _plate_data, seeds) = detailed_gen.get_plate_data();
    
    println!("World Details:");
    println!("  Resolution: {}×{} pixels", width, height);
    println!("  Coverage: 0.2° per pixel");
    println!("  Total plates: {}", stats.len());
    
    // Find largest and smallest plates
    let mut sorted_stats: Vec<_> = stats.into_iter().collect();
    sorted_stats.sort_by(|a, b| b.1.area_km2.cmp(&a.1.area_km2));
    
    println!("\nLargest plates:");
    for (plate_id, stat) in sorted_stats.iter().take(5) {
        println!("  Plate {:2}: {:6.2}% ({:>10} km²) - Motion: {:5.1}°@{:.1}cm/yr", 
                 plate_id, 
                 stat.percentage, 
                 stat.area_km2,
                 stat.seed.motion_direction,
                 stat.seed.motion_speed);
    }

    println!("\nSmallest plates:");
    for (plate_id, stat) in sorted_stats.iter().rev().take(3) {
        println!("  Plate {:2}: {:6.2}% ({:>10} km²)", 
                 plate_id, stat.percentage, stat.area_km2);
    }

    // Export data
    println!("\n💾 Exporting Data");
    let raw_bytes = detailed_gen.export_raw_u16_le();
    let mut file = File::create("geoforge_plates.bin")?;
    file.write_all(&raw_bytes)?;
    println!("Exported {} bytes to geoforge_plates.bin", raw_bytes.len());

    // Show GeoTIFF metadata
    let transform = detailed_gen.get_geotransform();
    println!("GeoTIFF transform: {:?}", transform);
    
    println!("\n🎉 All examples completed successfully!");
    println!("Seed for reproduction: {}", detailed_gen.get_seed());
    
    Ok(())
}
