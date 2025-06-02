/// Simple example showing basic tectonic plate generation
use geoforge::TectonicPlateGenerator;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Simple Tectonic Plate Generation Example");
    println!("========================================");

    // Create a small world for demonstration
    let mut generator = TectonicPlateGenerator::with_seed(360, 180, 6, 12345)?;
    
    println!("Generating world with seed: {}", generator.get_seed());
    println!("Resolution: 1° per pixel (360×180 grid)");
    println!("Number of plates: 6");
    
    // Generate the plates
    let _plate_map = generator.generate("region_growing", true)?;
    
    // Validate the result
    generator.validate()?;
    
    // Show statistics
    let stats = generator.get_plate_stats();
    println!("\nGenerated {} plates:", stats.len());
    
    let mut sorted_stats: Vec<_> = stats.into_iter().collect();
    sorted_stats.sort_by(|a, b| b.1.area_km2.cmp(&a.1.area_km2));
    
    for (plate_id, stat) in sorted_stats {
        println!(
            "Plate {}: {:5.1}% of Earth's surface ({:>8} km²)",
            plate_id, stat.percentage, stat.area_km2
        );
    }
    
    // Show world dimensions
    let (width, height, plate_data, seeds) = generator.get_plate_data();
    println!("\nWorld data:");
    println!("  Grid size: {}×{} pixels", width, height);
    println!("  Total data points: {}", plate_data.len());
    println!("  Plate seeds: {}", seeds.len());
    println!("  Memory usage: {:.1} KB", plate_data.len() * 2 / 1024);
    
    println!("\n✅ Example completed successfully!");
    
    Ok(())
}
