/// Example demonstrating different planetary parameters
use geoforge::{WorldMap, PlanetaryParams};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🪐 Planetary Parameters Example");
    println!("===============================");
    println!("Demonstrating world generation with different planetary characteristics");
    
    // Example 1: Earth-like world (default)
    println!("\n🌍 Example 1: Earth-like World");
    println!("------------------------------");
    
    let mut earth_world = WorldMap::new(90, 45, 42)?;
    earth_world.generate_tectonics(6, true)?;
    
    println!("🌍 Earth Parameters:");
    println!("   Radius: {:.0} km", earth_world.planetary_params.radius_km);
    println!("   Gravity: {:.2} m/s²", earth_world.planetary_params.gravity_ms2);
    println!("   Day length: {:.1} hours", earth_world.planetary_params.rotation_period_hours);
    println!("   Axial tilt: {:.1}°", earth_world.planetary_params.axial_tilt_degrees);
    println!("   Year length: {:.1} days", earth_world.planetary_params.orbital_period_days);
    println!("   Orbital distance: {:.3}-{:.3} AU", earth_world.planetary_params.perihelion_au, earth_world.planetary_params.aphelion_au);
    println!("   Atmospheric pressure: {:.1} kPa", earth_world.planetary_params.atmospheric_pressure_kpa);
    println!("   Greenhouse factor: {:.2}x Earth", earth_world.planetary_params.greenhouse_factor);
    
    if let Some(stats) = earth_world.get_tectonic_stats() {
        let total_area: u64 = stats.values().map(|s| s.area_km2).sum();
        println!("   🗺️ Total surface area: {} km²", total_area);
        println!("   🧩 Number of plates: {}", stats.len());
    }
    
    // Example 2: Mars-like world
    println!("\n🔴 Example 2: Mars-like World");
    println!("-----------------------------");
    
    let mut mars_world = WorldMap::new_mars(90, 45, 123)?;
    mars_world.generate_tectonics(5, true)?;
    
    println!("🔴 Mars Parameters:");
    println!("   Radius: {:.0} km", mars_world.planetary_params.radius_km);
    println!("   Gravity: {:.2} m/s²", mars_world.planetary_params.gravity_ms2);
    println!("   Day length: {:.1} hours", mars_world.planetary_params.rotation_period_hours);
    println!("   Axial tilt: {:.1}°", mars_world.planetary_params.axial_tilt_degrees);
    println!("   Year length: {:.0} days ({:.1} Earth years)", mars_world.planetary_params.orbital_period_days, mars_world.planetary_params.year_length_factor());
    println!("   Orbital distance: {:.3}-{:.3} AU", mars_world.planetary_params.perihelion_au, mars_world.planetary_params.aphelion_au);
    println!("   Atmospheric pressure: {:.3} kPa", mars_world.planetary_params.atmospheric_pressure_kpa);
    println!("   Greenhouse factor: {:.2}x Earth", mars_world.planetary_params.greenhouse_factor);
    
    if let Some(stats) = mars_world.get_tectonic_stats() {
        let total_area: u64 = stats.values().map(|s| s.area_km2).sum();
        println!("   🗺️ Total surface area: {} km²", total_area);
        println!("   🧩 Number of plates: {}", stats.len());
    }
    
    // Example 3: Venus-like world
    println!("\n🟡 Example 3: Venus-like World");
    println!("------------------------------");
    
    let mut venus_world = WorldMap::new_venus(90, 45, 456)?;
    venus_world.generate_tectonics(4, true)?;
    
    println!("🟡 Venus Parameters:");
    println!("   Radius: {:.0} km", venus_world.planetary_params.radius_km);
    println!("   Gravity: {:.2} m/s²", venus_world.planetary_params.gravity_ms2);
    println!("   Day length: {:.0} hours ({:.0} Earth days)", 
             venus_world.planetary_params.rotation_period_hours,
             venus_world.planetary_params.rotation_period_hours / 24.0);
    println!("   Axial tilt: {:.1}° (retrograde)", venus_world.planetary_params.axial_tilt_degrees);
    println!("   Year length: {:.1} days ({:.2} Earth years)", venus_world.planetary_params.orbital_period_days, venus_world.planetary_params.year_length_factor());
    println!("   Orbital distance: {:.3}-{:.3} AU", venus_world.planetary_params.perihelion_au, venus_world.planetary_params.aphelion_au);
    println!("   Atmospheric pressure: {:.0} kPa ({:.0} bar)", 
             venus_world.planetary_params.atmospheric_pressure_kpa,
             venus_world.planetary_params.atmospheric_pressure_kpa / 100.0);
    println!("   Greenhouse factor: {:.1}x Earth", venus_world.planetary_params.greenhouse_factor);
    
    if let Some(stats) = venus_world.get_tectonic_stats() {
        let total_area: u64 = stats.values().map(|s| s.area_km2).sum();
        println!("   🗺️ Total surface area: {} km²", total_area);
        println!("   🧩 Number of plates: {}", stats.len());
    }
    
    // Example 4: Custom super-Earth
    println!("\n🌌 Example 4: Custom Super-Earth");
    println!("---------------------------------");
    
    let custom_params = PlanetaryParams::from_radius_and_density(8000.0, 6000.0); // Larger, denser planet
    let mut custom_params = custom_params;
    custom_params.axial_tilt_degrees = 35.0; // High axial tilt
    custom_params.rotation_period_hours = 18.0; // Shorter day
    custom_params.orbital_period_days = 480.0; // Longer year
    custom_params.orbital_eccentricity = 0.05; // More eccentric orbit
    custom_params.semi_major_axis_au = 1.2; // Farther from star
    custom_params.perihelion_au = 1.14; // 1.2 * (1 - 0.05)
    custom_params.aphelion_au = 1.26; // 1.2 * (1 + 0.05)
    custom_params.atmospheric_pressure_kpa = 150.0; // Thicker atmosphere
    custom_params.greenhouse_factor = 1.36; // 45°C / 33°C = 1.36x Earth
    custom_params.stellar_luminosity = 1.2; // Brighter star (20% more luminous than Sun)
    
    let mut custom_world = WorldMap::new_with_params(90, 45, 789, custom_params)?;
    custom_world.generate_tectonics(8, true)?;
    
    println!("🌌 Super-Earth Parameters:");
    println!("   Radius: {:.0} km", custom_world.planetary_params.radius_km);
    println!("   Gravity: {:.2} m/s²", custom_world.planetary_params.gravity_ms2);
    println!("   Mass: {:.2e} kg", custom_world.planetary_params.mass_kg);
    println!("   Day length: {:.1} hours", custom_world.planetary_params.rotation_period_hours);
    println!("   Axial tilt: {:.1}°", custom_world.planetary_params.axial_tilt_degrees);
    println!("   Year length: {:.0} days ({:.1} Earth years)", custom_world.planetary_params.orbital_period_days, custom_world.planetary_params.year_length_factor());
    println!("   Orbital distance: {:.2}-{:.2} AU", custom_world.planetary_params.perihelion_au, custom_world.planetary_params.aphelion_au);
    println!("   Atmospheric pressure: {:.0} kPa", custom_world.planetary_params.atmospheric_pressure_kpa);
    println!("   Greenhouse factor: {:.2}x Earth", custom_world.planetary_params.greenhouse_factor);
    println!("   Stellar luminosity: {:.1}x solar", custom_world.planetary_params.stellar_luminosity);
    println!("   Escape velocity: {:.1} km/s", custom_world.planetary_params.escape_velocity_kms());
    
    if let Some(stats) = custom_world.get_tectonic_stats() {
        let total_area: u64 = stats.values().map(|s| s.area_km2).sum();
        println!("   🗺️ Total surface area: {} km²", total_area);
        println!("   🧩 Number of plates: {}", stats.len());
    }
    
    // Example 5: Climate factors analysis
    println!("\n🌡️ Example 5: Climate Factor Analysis");
    println!("--------------------------------------");
    
    println!("Seasonal variation factors (based on axial tilt):");
    println!("   🌍 Earth: {:.2}", earth_world.planetary_params.seasonal_variation_factor());
    println!("   🔴 Mars: {:.2}", mars_world.planetary_params.seasonal_variation_factor());
    println!("   🟡 Venus: {:.2}", venus_world.planetary_params.seasonal_variation_factor());
    println!("   🌌 Super-Earth: {:.2}", custom_world.planetary_params.seasonal_variation_factor());
    
    println!("\nDiurnal (day/night) variation factors (based on rotation):");
    println!("   🌍 Earth: {:.2}", earth_world.planetary_params.diurnal_variation_factor());
    println!("   🔴 Mars: {:.2}", mars_world.planetary_params.diurnal_variation_factor());
    println!("   🟡 Venus: {:.2}", venus_world.planetary_params.diurnal_variation_factor());
    println!("   🌌 Super-Earth: {:.2}", custom_world.planetary_params.diurnal_variation_factor());
    
    println!("\nOrbital radiation analysis:");
    println!("   🌍 Earth solar flux: {:.2}x Earth baseline", earth_world.planetary_params.average_solar_flux());
    println!("   🔴 Mars solar flux: {:.2}x Earth baseline", mars_world.planetary_params.average_solar_flux());
    println!("   🟡 Venus solar flux: {:.2}x Earth baseline", venus_world.planetary_params.average_solar_flux());
    println!("   🌌 Super-Earth solar flux: {:.2}x Earth baseline", custom_world.planetary_params.average_solar_flux());
    
    println!("\nOrbital eccentricity effects (seasonal radiation variation):");
    println!("   🌍 Earth: {:.1}% variation", earth_world.planetary_params.orbital_radiation_variation() * 100.0);
    println!("   🔴 Mars: {:.1}% variation", mars_world.planetary_params.orbital_radiation_variation() * 100.0);
    println!("   🟡 Venus: {:.1}% variation", venus_world.planetary_params.orbital_radiation_variation() * 100.0);
    println!("   🌌 Super-Earth: {:.1}% variation", custom_world.planetary_params.orbital_radiation_variation() * 100.0);
    
    // Example 6: Insolation calculations
    println!("\n☀️ Example 6: Actual Insolation Values");
    println!("--------------------------------------");
    
    println!("Average insolation (W/m²):");
    println!("   🌍 Earth: {:.0} W/m²", earth_world.planetary_params.average_insolation_wm2());
    println!("   🔴 Mars: {:.0} W/m²", mars_world.planetary_params.average_insolation_wm2());
    println!("   🟡 Venus: {:.0} W/m²", venus_world.planetary_params.average_insolation_wm2());
    println!("   🌌 Super-Earth: {:.0} W/m² (brighter star compensates for distance)", custom_world.planetary_params.average_insolation_wm2());
    
    println!("\nSeasonal insolation range due to orbital eccentricity:");
    println!("   🌍 Earth: {:.0}-{:.0} W/m²", earth_world.planetary_params.aphelion_insolation_wm2(), earth_world.planetary_params.perihelion_insolation_wm2());
    println!("   🔴 Mars: {:.0}-{:.0} W/m²", mars_world.planetary_params.aphelion_insolation_wm2(), mars_world.planetary_params.perihelion_insolation_wm2());
    println!("   🟡 Venus: {:.0}-{:.0} W/m²", venus_world.planetary_params.aphelion_insolation_wm2(), venus_world.planetary_params.perihelion_insolation_wm2());
    println!("   🌌 Super-Earth: {:.0}-{:.0} W/m²", custom_world.planetary_params.aphelion_insolation_wm2(), custom_world.planetary_params.perihelion_insolation_wm2());
    
    println!("\nInsolation at different orbital positions for Mars (most eccentric):");
    for i in 0..=10 {
        let position = i as f64 / 10.0;
        let insolation = mars_world.planetary_params.insolation_at_position(position);
        let season = match i {
            0 => "Perihelion (closest)",
            3 => "Spring equinox",
            5 => "Summer solstice", 
            8 => "Autumn equinox",
            10 => "Aphelion (farthest)",
            _ => "",
        };
        if !season.is_empty() {
            println!("   Position {:.1}: {:.0} W/m² ({})", position, insolation, season);
        }
    }
    
    // Export examples
    #[cfg(feature = "export-png")]
    {
        println!("\n📤 Exporting planetary examples...");
        
        std::fs::create_dir_all("outputs/examples/planetary_params")?;
        
        earth_world.export_tectonics_png("outputs/examples/planetary_params", "earth_tectonics.png")?;
        mars_world.export_tectonics_png("outputs/examples/planetary_params", "mars_tectonics.png")?;
        venus_world.export_tectonics_png("outputs/examples/planetary_params", "venus_tectonics.png")?;
        custom_world.export_tectonics_png("outputs/examples/planetary_params", "super_earth_tectonics.png")?;
        
        println!("✅ Exported PNG visualizations to outputs/examples/planetary_params/");
    }
    
    #[cfg(not(feature = "export-png"))]
    {
        println!("\n⚠️ PNG export features not enabled. Run with --features export-png for visualizations");
    }
    
    println!("\n🎉 Planetary parameters example completed!");
    println!("\n💡 Key insights:");
    println!("  • Different planetary sizes affect total surface area and plate scaling");
    println!("  • Gravity affects geological processes and mountain building");
    println!("  • Rotation period affects day/night temperature cycles");
    println!("  • Axial tilt controls seasonal variation strength");
    println!("  • Atmospheric parameters will influence climate generation");
    
    Ok(())
}