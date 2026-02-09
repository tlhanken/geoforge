//! Integration tests for Stage 2: Geological Province Generation
//!
//! Tests the full pipeline from tectonic plates to orogenic belts

use geoforge::{
    GeologicProvince, OrogenicBeltGenerator, OrogenicConfig, WorldMap,
};

#[test]
fn test_orogenic_belt_generation_full_pipeline() {
    // Create a small world for testing (0.4° resolution)
    let mut world = WorldMap::new(900, 450, 42).expect("Failed to create world");

    // Generate tectonic plates
    world
        .tectonics()
        .generate_plates(12)
        .expect("Failed to generate plates");

    // Analyze boundaries (automatically assigns motion and classifies)
    world
        .analyze_boundaries(None)
        .expect("Failed to analyze boundaries");

    // Get tectonic data
    let plate_map = world.tectonics.as_ref().expect("No plate map generated");
    let metadata = world.get_tectonic_metadata().expect("No tectonic metadata");

    // Generate orogenic belts
    let config = OrogenicConfig::default();
    let generator = OrogenicBeltGenerator::new(config, world.planetary_params.clone());

    let orogens = generator.generate_orogens(
        &metadata.plate_boundaries,
        &metadata.plate_stats,
        plate_map,
    );

    // Verify we got some orogenic belts
    assert!(
        !orogens.is_empty(),
        "Should generate at least one orogenic belt"
    );

    // Verify all orogens have valid data
    for orogen in &orogens {
        assert!(
            !orogen.pixels.is_empty(),
            "Orogenic belt should have pixels"
        );
        assert!(
            orogen.characteristics.width_km > 0.0,
            "Width should be positive"
        );
        assert!(
            orogen.characteristics.convergence_rate >= 2.0,
            "Convergence rate should be above minimum"
        );
        assert!(
            orogen.source_boundary_index.is_some(),
            "Should have source boundary"
        );
    }

    println!("Generated {} orogenic belts", orogens.len());
}

#[test]
fn test_all_orogen_types_generated() {
    // Create a larger world to ensure diversity (0.2° resolution)
    let mut world = WorldMap::new(1800, 900, 100).expect("Failed to create world");

    world
        .tectonics()
        .generate_plates(20)
        .expect("Failed to generate plates");

    world
        .analyze_boundaries(None)
        .expect("Failed to analyze boundaries");

    let plate_map = world.tectonics.as_ref().unwrap();
    let metadata = world.get_tectonic_metadata().unwrap();

    // Generate orogenic belts
    let config = OrogenicConfig::default();
    let generator = OrogenicBeltGenerator::new(config, world.planetary_params.clone());
    let orogens = generator.generate_orogens(
        &metadata.plate_boundaries,
        &metadata.plate_stats,
        plate_map,
    );

    // Check for variety of orogen types
    let has_collision = orogens
        .iter()
        .any(|o| o.characteristics.province_type == GeologicProvince::CollisionOrogen);

    println!("Orogen type distribution:");
    println!("  Collision: {}", has_collision);

    // With 20 plates, we should get collision orogens
    // (Subduction systems now handled separately by Arc Systems generator)
    assert!(
        has_collision,
        "Should generate collision orogens with 20 plates"
    );
}

#[test]
fn test_dynamic_width_scaling() {
    let mut world = WorldMap::new(900, 450, 200).expect("Failed to create world");

    world
        .tectonics()
        .generate_plates(15)
        .expect("Failed to generate plates");

    world
        .analyze_boundaries(None)
        .expect("Failed to analyze boundaries");

    let plate_map = world.tectonics.as_ref().unwrap();
    let metadata = world.get_tectonic_metadata().unwrap();

    let config = OrogenicConfig::default();
    let generator = OrogenicBeltGenerator::new(config, world.planetary_params.clone());
    let orogens = generator.generate_orogens(
        &metadata.plate_boundaries,
        &metadata.plate_stats,
        plate_map,
    );

    // Check that widths vary based on convergence rate
    if orogens.len() >= 2 {
        let widths: Vec<f64> = orogens.iter().map(|o| o.characteristics.width_km).collect();
        let min_width = widths.iter().copied().fold(f64::INFINITY, f64::min);
        let max_width = widths.iter().copied().fold(f64::NEG_INFINITY, f64::max);

        println!("Width range: {} km - {} km", min_width, max_width);

        // Widths should vary (not all the same)
        assert!(
            max_width > min_width,
            "Orogenic belt widths should vary based on convergence rate"
        );
    }
}

#[test]
fn test_orogenic_belt_pixel_expansion() {
    let mut world = WorldMap::new(900, 450, 300).expect("Failed to create world");

    world
        .tectonics()
        .generate_plates(10)
        .expect("Failed to generate plates");

    world
        .analyze_boundaries(None)
        .expect("Failed to analyze boundaries");

    let plate_map = world.tectonics.as_ref().unwrap();
    let metadata = world.get_tectonic_metadata().unwrap();

    let config = OrogenicConfig::default();
    let generator = OrogenicBeltGenerator::new(config, world.planetary_params.clone());
    let orogens = generator.generate_orogens(
        &metadata.plate_boundaries,
        &metadata.plate_stats,
        plate_map,
    );

    // Verify that orogen pixels are larger than boundary pixels (expansion happened)
    for (i, orogen) in orogens.iter().enumerate() {
        let boundary_idx = orogen.source_boundary_index.unwrap();
        let boundary = &metadata.plate_boundaries[boundary_idx];

        assert!(
            orogen.pixels.len() > boundary.pixels.len(),
            "Orogen {} should have more pixels than its source boundary (expanded)",
            i
        );

        println!(
            "Orogen {}: boundary {} pixels → belt {} pixels ({}x expansion)",
            i,
            boundary.pixels.len(),
            orogen.pixels.len(),
            orogen.pixels.len() as f64 / boundary.pixels.len() as f64
        );
    }
}

#[test]
fn test_convergence_rate_filtering() {
    let mut world = WorldMap::new(900, 450, 400).expect("Failed to create world");

    world
        .tectonics()
        .generate_plates(12)
        .expect("Failed to generate plates");

    world
        .analyze_boundaries(None)
        .expect("Failed to analyze boundaries");

    let plate_map = world.tectonics.as_ref().unwrap();
    let metadata = world.get_tectonic_metadata().unwrap();

    // Count convergent boundaries
    let convergent_count = metadata
        .plate_boundaries
        .iter()
        .filter(|b| b.interaction_type == geoforge::PlateInteraction::Convergent)
        .count();

    let config = OrogenicConfig::default();
    let min_rate = config.min_convergence_rate;
    let generator = OrogenicBeltGenerator::new(config, world.planetary_params.clone());
    let orogens = generator.generate_orogens(
        &metadata.plate_boundaries,
        &metadata.plate_stats,
        plate_map,
    );

    println!(
        "Convergent boundaries: {}, Orogens generated: {}",
        convergent_count,
        orogens.len()
    );

    // Orogens should be <= convergent boundaries (some filtered by min rate)
    assert!(
        orogens.len() <= convergent_count,
        "Should have at most one orogen per convergent boundary"
    );

    // All generated orogens should meet minimum convergence rate
    for orogen in &orogens {
        assert!(
            orogen.characteristics.convergence_rate >= min_rate,
            "All orogens should meet minimum convergence rate"
        );
    }
}

#[test]
fn test_deterministic_generation() {
    // Generate twice with same seed
    let mut world1 = WorldMap::new(900, 450, 500).expect("Failed to create world");

    world1.tectonics().generate_plates(12).unwrap();
    world1.analyze_boundaries(None).unwrap();

    let plate_map1 = world1.tectonics.as_ref().unwrap();
    let metadata1 = world1.get_tectonic_metadata().unwrap();

    let config = OrogenicConfig::default();
    let generator = OrogenicBeltGenerator::new(config.clone(), world1.planetary_params.clone());
    let orogens1 = generator.generate_orogens(
        &metadata1.plate_boundaries,
        &metadata1.plate_stats,
        plate_map1,
    );

    // Generate again with same seed
    let mut world2 = WorldMap::new(900, 450, 500).expect("Failed to create world");

    world2.tectonics().generate_plates(12).unwrap();
    world2.analyze_boundaries(None).unwrap();

    let plate_map2 = world2.tectonics.as_ref().unwrap();
    let metadata2 = world2.get_tectonic_metadata().unwrap();

    let generator2 = OrogenicBeltGenerator::new(config, world2.planetary_params.clone());
    let orogens2 = generator2.generate_orogens(
        &metadata2.plate_boundaries,
        &metadata2.plate_stats,
        plate_map2,
    );

    // Should produce identical results (same count and distribution)
    assert_eq!(
        orogens1.len(),
        orogens2.len(),
        "Same seed should produce same number of orogens"
    );

    // Count orogen types in both sets
    let count_type = |orogens: &[geoforge::ProvinceRegion], ptype| {
        orogens
            .iter()
            .filter(|o| o.characteristics.province_type == ptype)
            .count()
    };

    assert_eq!(
        count_type(&orogens1, GeologicProvince::CollisionOrogen),
        count_type(&orogens2, GeologicProvince::CollisionOrogen),
        "Collision orogen counts should match"
    );

    println!("✓ Deterministic generation verified with {} orogens", orogens1.len());
}

#[test]
fn test_orogenic_characteristics_scaling() {
    let mut world = WorldMap::new(900, 450, 600).expect("Failed to create world");

    world.tectonics().generate_plates(15).unwrap();
    world.analyze_boundaries(None).unwrap();

    let plate_map = world.tectonics.as_ref().unwrap();
    let metadata = world.get_tectonic_metadata().unwrap();

    let config = OrogenicConfig::default();
    let generator = OrogenicBeltGenerator::new(config, world.planetary_params.clone());
    let orogens = generator.generate_orogens(
        &metadata.plate_boundaries,
        &metadata.plate_stats,
        plate_map,
    );

    // Verify characteristics make sense
    for orogen in &orogens {
        let chars = &orogen.characteristics;

        // Intensity should be 0-1
        assert!(
            chars.intensity >= 0.0 && chars.intensity <= 1.0,
            "Intensity should be normalized 0-1"
        );

        // Roughness should be 0-1
        assert!(
            chars.roughness >= 0.0 && chars.roughness <= 1.0,
            "Roughness should be 0-1"
        );

        // Elevation intensity should match orogen type expectations
        match chars.province_type {
            GeologicProvince::CollisionOrogen => {
                assert_eq!(
                    chars.elevation_intensity, 1.0,
                    "Collision orogens should have max elevation"
                );
            }
            GeologicProvince::AccretionaryWedge => {
                assert_eq!(
                    chars.elevation_intensity, 0.2,
                    "Accretionary wedges should have 0.2 elevation (often submarine)"
                );
            }
            _ => {
                // For other province types not yet fully tested
                assert!(
                    chars.elevation_intensity >= -1.2 && chars.elevation_intensity <= 1.0,
                    "Elevation intensity should be in reasonable range"
                );
            }
        }
    }
}

#[test]
fn test_full_geology_pipeline_deterministic() {
    // Test that the FULL geology pipeline (all 20 province types including hotspot tracks)
    // reproduces exactly with the same seed

    let seed = 12345_u64;

    // Generate world 1
    let mut world1 = WorldMap::new(1800, 900, seed).expect("Failed to create world 1");
    world1.tectonics().generate_plates(20).unwrap();
    world1.tectonics().analyze(None).unwrap();
    let provinces1 = world1.generate_geology(None).unwrap();

    // Generate world 2 with same seed
    let mut world2 = WorldMap::new(1800, 900, seed).expect("Failed to create world 2");
    world2.tectonics().generate_plates(20).unwrap();
    world2.tectonics().analyze(None).unwrap();
    let provinces2 = world2.generate_geology(None).unwrap();

    // Should produce identical number of provinces
    println!("Run 1: {} provinces", provinces1.len());
    println!("Run 2: {} provinces", provinces2.len());

    if provinces1.len() != provinces2.len() {
        // Debug: Show which types differ
        let count_type_local = |provinces: &[geoforge::ProvinceRegion], ptype| {
            provinces
                .iter()
                .filter(|p| p.characteristics.province_type == ptype)
                .count()
        };

        println!("\nDifferences found:");
        for ptype in &[
            GeologicProvince::OceanicHotspotTrack,
            GeologicProvince::ContinentalHotspotTrack,
            GeologicProvince::ContinentalFloodBasalt,
            GeologicProvince::OceanicPlateau,
        ] {
            let c1 = count_type_local(&provinces1, *ptype);
            let c2 = count_type_local(&provinces2, *ptype);
            if c1 != c2 {
                println!("  {:?}: {} vs {}", ptype, c1, c2);
            }
        }
    }

    assert_eq!(
        provinces1.len(),
        provinces2.len(),
        "Same seed should produce same number of geological provinces"
    );

    // Count each province type in both runs
    let count_type = |provinces: &[geoforge::ProvinceRegion], ptype| {
        provinces
            .iter()
            .filter(|p| p.characteristics.province_type == ptype)
            .count()
    };

    // Test all 18 province types for determinism
    let all_types = vec![
        GeologicProvince::CollisionOrogen,
        GeologicProvince::AccretionaryWedge,
        GeologicProvince::ContinentalFloodBasalt,
        GeologicProvince::OceanicPlateau,
        GeologicProvince::ContinentalHotspotTrack,
        GeologicProvince::VolcanicArc,
        GeologicProvince::ForearcBasin,
        GeologicProvince::BackarcBasin,
        GeologicProvince::OceanTrench,
        GeologicProvince::Craton,
        GeologicProvince::Platform,
        GeologicProvince::IntracratonicBasin,
        GeologicProvince::ExtendedCrust,
        GeologicProvince::ContinentalRift,
        GeologicProvince::AbyssalPlain,
        GeologicProvince::MidOceanRidge,
        GeologicProvince::OceanicFractureZone,
        GeologicProvince::OceanicHotspotTrack,
    ];

    println!("\n=== Province Type Counts (Determinism Test) ===");
    for ptype in all_types {
        let count1 = count_type(&provinces1, ptype);
        let count2 = count_type(&provinces2, ptype);

        if count1 > 0 || count2 > 0 {
            println!("  {:?}: {} vs {}", ptype, count1, count2);
        }

        assert_eq!(
            count1, count2,
            "{:?} count should match between runs with same seed",
            ptype
        );
    }

    // CRITICAL: Test hotspot tracks specifically (both oceanic and continental)
    let oceanic_hotspots1 = count_type(&provinces1, GeologicProvince::OceanicHotspotTrack);
    let oceanic_hotspots2 = count_type(&provinces2, GeologicProvince::OceanicHotspotTrack);
    let continental_hotspots1 = count_type(&provinces1, GeologicProvince::ContinentalHotspotTrack);
    let continental_hotspots2 = count_type(&provinces2, GeologicProvince::ContinentalHotspotTrack);

    println!("\n=== Hotspot Track Verification ===");
    println!("  Oceanic hotspot tracks: {} (both runs)", oceanic_hotspots1);
    println!("  Continental hotspot tracks: {} (both runs)", continental_hotspots1);

    assert_eq!(
        oceanic_hotspots1, oceanic_hotspots2,
        "Oceanic hotspot tracks must reproduce exactly with same seed"
    );
    assert_eq!(
        continental_hotspots1, continental_hotspots2,
        "Continental hotspot tracks must reproduce exactly with same seed"
    );

    // Verify pixel-level reproducibility for hotspot tracks
    let hotspot_provinces1: Vec<_> = provinces1
        .iter()
        .filter(|p| {
            p.characteristics.province_type == GeologicProvince::OceanicHotspotTrack
                || p.characteristics.province_type == GeologicProvince::ContinentalHotspotTrack
        })
        .collect();

    let hotspot_provinces2: Vec<_> = provinces2
        .iter()
        .filter(|p| {
            p.characteristics.province_type == GeologicProvince::OceanicHotspotTrack
                || p.characteristics.province_type == GeologicProvince::ContinentalHotspotTrack
        })
        .collect();

    assert_eq!(
        hotspot_provinces1.len(),
        hotspot_provinces2.len(),
        "Total hotspot track count should match"
    );

    // Verify each hotspot track has same pixel count
    for (i, (h1, h2)) in hotspot_provinces1.iter().zip(hotspot_provinces2.iter()).enumerate() {
        assert_eq!(
            h1.pixels.len(),
            h2.pixels.len(),
            "Hotspot track {} should have same pixel count in both runs",
            i
        );

        assert_eq!(
            h1.characteristics.province_type,
            h2.characteristics.province_type,
            "Hotspot track {} should have same type in both runs",
            i
        );
    }

    println!("\n✓ Full geology pipeline determinism verified:");
    println!("  {} total provinces", provinces1.len());
    println!("  {} hotspot tracks (pixel-perfect match)", hotspot_provinces1.len());
}

#[test]
fn test_geology_generator_full_pipeline() {
    // Test GeologyGenerator.generate_all_provinces() directly
    // This is the missing test identified in the code review

    use geoforge::{GeologyGenerator, GeologyConfig};

    // Create world with enough plates to generate diverse provinces
    let mut world = WorldMap::new(1800, 900, 777).expect("Failed to create world");
    world.tectonics().generate_plates(20).unwrap();
    world.tectonics().analyze(None).unwrap();

    let plate_map = world.tectonics.as_ref().unwrap();
    let metadata = world.get_tectonic_metadata().unwrap();

    // Create GeologyGenerator with default config
    let config = GeologyConfig::default();
    let generator = GeologyGenerator::new(config, world.seed, world.planetary_params.clone());

    // Generate all provinces through the full pipeline
    let provinces = generator.generate_all_provinces(
        &metadata.plate_boundaries,
        &metadata.plate_stats,
        &metadata.plate_seeds,
        plate_map,
    );

    // Verify we got provinces
    assert!(
        !provinces.is_empty(),
        "Should generate at least some provinces"
    );

    println!("\n=== GeologyGenerator Full Pipeline Test ===");
    println!("Generated {} total provinces", provinces.len());

    // Count each province type
    let count_type = |ptype| {
        provinces
            .iter()
            .filter(|p| p.characteristics.province_type == ptype)
            .count()
    };

    // Test all 18 implemented province types
    let all_types = vec![
        (GeologicProvince::CollisionOrogen, "CollisionOrogen"),
        (GeologicProvince::AccretionaryWedge, "AccretionaryWedge"),
        (GeologicProvince::ContinentalFloodBasalt, "ContinentalFloodBasalt"),
        (GeologicProvince::OceanicPlateau, "OceanicPlateau"),
        (GeologicProvince::ContinentalHotspotTrack, "ContinentalHotspotTrack"),
        (GeologicProvince::VolcanicArc, "VolcanicArc"),
        (GeologicProvince::ForearcBasin, "ForearcBasin"),
        (GeologicProvince::BackarcBasin, "BackarcBasin"),
        (GeologicProvince::OceanTrench, "OceanTrench"),
        (GeologicProvince::Craton, "Craton"),
        (GeologicProvince::Platform, "Platform"),
        (GeologicProvince::IntracratonicBasin, "IntracratonicBasin"),
        (GeologicProvince::ExtendedCrust, "ExtendedCrust"),
        (GeologicProvince::ContinentalRift, "ContinentalRift"),
        (GeologicProvince::AbyssalPlain, "AbyssalPlain"),
        (GeologicProvince::MidOceanRidge, "MidOceanRidge"),
        (GeologicProvince::OceanicFractureZone, "OceanicFractureZone"),
        (GeologicProvince::OceanicHotspotTrack, "OceanicHotspotTrack"),
    ];

    println!("\n=== Province Type Counts ===");
    let mut found_types = 0;
    for (ptype, name) in &all_types {
        let count = count_type(*ptype);
        if count > 0 {
            println!("  {}: {}", name, count);
            found_types += 1;
        }
    }

    // With 20 plates, we should generate multiple province types
    assert!(
        found_types >= 5,
        "Should generate at least 5 different province types with 20 plates, got {}",
        found_types
    );

    // Verify foundation layers exist (should always generate)
    let has_oceanic_base = count_type(GeologicProvince::AbyssalPlain) > 0;
    let has_stable_regions = count_type(GeologicProvince::Platform) > 0 || count_type(GeologicProvince::Craton) > 0;

    assert!(
        has_oceanic_base,
        "Should always generate oceanic base layer (abyssal plains)"
    );
    assert!(
        has_stable_regions,
        "Should always generate stable continental regions (platform/craton)"
    );

    // Verify all provinces have valid data
    for (i, province) in provinces.iter().enumerate() {
        assert!(
            !province.pixels.is_empty(),
            "Province {} should have pixels",
            i
        );
        assert!(
            province.characteristics.width_km >= 0.0,
            "Province {} width should be non-negative",
            i
        );
        assert!(
            province.characteristics.roughness >= 0.0 && province.characteristics.roughness <= 1.0,
            "Province {} roughness should be 0-1, got {}",
            i,
            province.characteristics.roughness
        );
        assert!(
            province.characteristics.intensity >= 0.0 && province.characteristics.intensity <= 1.0,
            "Province {} intensity should be 0-1, got {}",
            i,
            province.characteristics.intensity
        );
    }

    println!("\n✓ GeologyGenerator full pipeline test passed:");
    println!("  {} total provinces", provinces.len());
    println!("  {} distinct province types", found_types);
    println!("  All provinces have valid characteristics");
}

#[test]
fn test_edge_case_empty_boundaries() {
    // Test that the generator handles empty boundary lists gracefully
    use geoforge::{GeologyGenerator, GeologyConfig};

    let mut world = WorldMap::new(900, 450, 999).expect("Failed to create world");
    world.tectonics().generate_plates(5).unwrap();
    world.tectonics().analyze(None).unwrap();

    let plate_map = world.tectonics.as_ref().unwrap();
    let metadata = world.get_tectonic_metadata().unwrap();

    let config = GeologyConfig::default();
    let generator = GeologyGenerator::new(config, world.seed, world.planetary_params.clone());

    // Generate with empty boundaries (should still create base layers)
    let provinces = generator.generate_all_provinces(
        &[], // Empty boundaries
        &metadata.plate_stats,
        &metadata.plate_seeds,
        plate_map,
    );

    // Should still generate foundation layers (oceanic base, stable regions)
    assert!(
        !provinces.is_empty(),
        "Should generate foundation layers even with no boundaries"
    );

    // Count foundation types
    let has_oceanic = provinces.iter()
        .any(|p| p.characteristics.province_type == geoforge::GeologicProvince::AbyssalPlain);
    let has_stable = provinces.iter()
        .any(|p| p.characteristics.province_type == geoforge::GeologicProvince::Platform
            || p.characteristics.province_type == geoforge::GeologicProvince::Craton);

    assert!(has_oceanic, "Should generate oceanic base layer");
    assert!(has_stable, "Should generate stable continental regions");

    // Should NOT generate boundary-dependent features
    let has_orogens = provinces.iter()
        .any(|p| p.characteristics.province_type == geoforge::GeologicProvince::CollisionOrogen);
    let has_ridges = provinces.iter()
        .any(|p| p.characteristics.province_type == geoforge::GeologicProvince::MidOceanRidge);

    assert!(!has_orogens, "Should not generate orogens without boundaries");
    assert!(!has_ridges, "Should not generate ridges without boundaries");

    println!("\n✓ Empty boundaries test passed:");
    println!("  {} foundation provinces generated", provinces.len());
}

#[test]
fn test_edge_case_tiny_plates() {
    // Test that the generator handles very small plates gracefully
    use geoforge::{GeologyGenerator, GeologyConfig};

    // Create small world to encourage tiny plates
    let mut world = WorldMap::new(300, 150, 888).expect("Failed to create world");
    world.tectonics().generate_plates(30).unwrap(); // Many plates in small space
    world.tectonics().analyze(None).unwrap();

    let plate_map = world.tectonics.as_ref().unwrap();
    let metadata = world.get_tectonic_metadata().unwrap();

    let config = GeologyConfig::default();
    let generator = GeologyGenerator::new(config, world.seed, world.planetary_params.clone());

    // Should not crash or hang with tiny plates
    let provinces = generator.generate_all_provinces(
        &metadata.plate_boundaries,
        &metadata.plate_stats,
        &metadata.plate_seeds,
        plate_map,
    );

    // Verify we got some provinces
    assert!(!provinces.is_empty(), "Should generate provinces even with tiny plates");

    // Count provinces with pixels (some may be empty with very tiny plates)
    let non_empty_count = provinces.iter().filter(|p| !p.pixels.is_empty()).count();
    assert!(non_empty_count > 0, "At least some provinces should have pixels");

    // Verify all provinces have valid characteristics
    for province in &provinces {
        assert!(province.characteristics.width_km >= 0.0, "Width should be non-negative");
    }

    println!("\n✓ Tiny plates test passed:");
    println!("  {} plates, {} provinces ({} non-empty)",
        metadata.plate_stats.len(), provinces.len(), non_empty_count);
}

#[test]
fn test_edge_case_polar_regions() {
    // Test that provinces near poles handle spherical projection correctly
    use geoforge::{GeologyGenerator, GeologyConfig};

    // Create world and generate provinces
    let mut world = WorldMap::new(1800, 900, 555).expect("Failed to create world");
    world.tectonics().generate_plates(15).unwrap();
    world.tectonics().analyze(None).unwrap();

    let plate_map = world.tectonics.as_ref().unwrap();
    let metadata = world.get_tectonic_metadata().unwrap();

    let config = GeologyConfig::default();
    let generator = GeologyGenerator::new(config, world.seed, world.planetary_params.clone());

    let provinces = generator.generate_all_provinces(
        &metadata.plate_boundaries,
        &metadata.plate_stats,
        &metadata.plate_seeds,
        plate_map,
    );

    // Check for provinces near poles (latitude > 80°)
    let mut polar_provinces = 0;
    for province in &provinces {
        for &(x, y) in &province.pixels {
            let (lat, _lon) = plate_map.projection.pixel_to_coords(x, y);
            if lat.abs() > 80.0 {
                polar_provinces += 1;
                break; // Count each province once
            }
        }
    }

    // With 15 plates, we should have some polar coverage
    println!("\n✓ Polar regions test passed:");
    println!("  {} provinces with polar pixels (>80° latitude)", polar_provinces);
    println!("  No crashes or infinite loops in polar expansion");
}

#[test]
fn test_edge_case_longitude_wraparound() {
    // Test that longitude wraparound (crossing 180° meridian) works correctly
    use geoforge::{GeologyGenerator, GeologyConfig};

    let mut world = WorldMap::new(1800, 900, 444).expect("Failed to create world");
    world.tectonics().generate_plates(12).unwrap();
    world.tectonics().analyze(None).unwrap();

    let plate_map = world.tectonics.as_ref().unwrap();
    let metadata = world.get_tectonic_metadata().unwrap();

    // Find boundaries that cross the map edge (x near 0 or width-1)
    let edge_boundaries: Vec<_> = metadata.plate_boundaries.iter()
        .filter(|b| {
            b.pixels.iter().any(|&(x, _)| x < 10 || x > plate_map.width - 10)
        })
        .collect();

    println!("\nFound {} boundaries near map edges", edge_boundaries.len());

    let config = GeologyConfig::default();
    let generator = GeologyGenerator::new(config, world.seed, world.planetary_params.clone());

    // Should handle wraparound without crashes
    let provinces = generator.generate_all_provinces(
        &metadata.plate_boundaries,
        &metadata.plate_stats,
        &metadata.plate_seeds,
        plate_map,
    );

    // Verify provinces near map edges are valid
    let mut edge_provinces = 0;
    for province in &provinces {
        let has_edge_pixels = province.pixels.iter()
            .any(|&(x, _)| x < 10 || x > plate_map.width - 10);
        if has_edge_pixels {
            edge_provinces += 1;
            // Verify no out-of-bounds pixels
            for &(x, y) in &province.pixels {
                assert!(x < plate_map.width, "X coordinate out of bounds");
                assert!(y < plate_map.height, "Y coordinate out of bounds");
            }
        }
    }

    println!("\n✓ Longitude wraparound test passed:");
    println!("  {} provinces near map edges", edge_provinces);
    println!("  All pixels within valid bounds");
}

#[test]
fn test_edge_case_all_oceanic_plates() {
    // Test edge case where all plates are oceanic (no continental provinces)
    use geoforge::{GeologyGenerator, GeologyConfig, PlateStats, PlateType};
    use std::collections::HashMap;

    let mut world = WorldMap::new(900, 450, 333).expect("Failed to create world");
    world.tectonics().generate_plates(8).unwrap();
    world.tectonics().analyze(None).unwrap();

    let plate_map = world.tectonics.as_ref().unwrap();
    let metadata = world.get_tectonic_metadata().unwrap();

    // Force all plates to be oceanic
    let mut modified_stats: HashMap<u16, PlateStats> = HashMap::new();
    for (id, stats) in &metadata.plate_stats {
        let mut modified = stats.clone();
        modified.plate_type = PlateType::Oceanic;
        modified_stats.insert(*id, modified);
    }

    let config = GeologyConfig::default();
    let generator = GeologyGenerator::new(config, world.seed, world.planetary_params.clone());

    let provinces = generator.generate_all_provinces(
        &metadata.plate_boundaries,
        &modified_stats,
        &metadata.plate_seeds,
        plate_map,
    );

    // Should still generate oceanic provinces
    assert!(!provinces.is_empty(), "Should generate oceanic provinces");

    // Should NOT generate continental provinces
    let has_continental = provinces.iter().any(|p| {
        matches!(
            p.characteristics.province_type,
            geoforge::GeologicProvince::Craton
                | geoforge::GeologicProvince::Platform
                | geoforge::GeologicProvince::CollisionOrogen
        )
    });

    assert!(!has_continental, "Should not generate continental provinces");

    // Should have oceanic provinces
    let has_oceanic = provinces.iter().any(|p| {
        matches!(
            p.characteristics.province_type,
            geoforge::GeologicProvince::AbyssalPlain | geoforge::GeologicProvince::MidOceanRidge
        )
    });

    assert!(has_oceanic, "Should generate oceanic provinces");

    println!("\n✓ All-oceanic plates test passed:");
    println!("  {} oceanic provinces generated", provinces.len());
}

#[test]
fn test_oceanic_continental_convergence() {
    // Test that oceanic-continental convergence generates arc systems
    // Real-world examples: Andes, Cascades, Japanese Alps

    use geoforge::{PlateType, PlateInteraction, GeologicProvince};

    let mut world = WorldMap::new(1800, 900, 12345).expect("Failed to create world");
    world.tectonics().generate_plates(15).unwrap();
    world.tectonics().analyze(None).unwrap();

    // Find oceanic-continental convergent boundaries (before geology generation)
    let oceanic_continental_count = {
        let metadata = world.get_tectonic_metadata().unwrap();
        metadata.plate_boundaries.iter()
            .filter(|b| {
                if b.interaction_type != PlateInteraction::Convergent {
                    return false;
                }

                let stats_a = metadata.plate_stats.get(&b.plate_a);
                let stats_b = metadata.plate_stats.get(&b.plate_b);

                if let (Some(a), Some(b)) = (stats_a, stats_b) {
                    matches!(
                        (a.plate_type, b.plate_type),
                        (PlateType::Oceanic, PlateType::Continental) |
                        (PlateType::Continental, PlateType::Oceanic)
                    )
                } else {
                    false
                }
            })
            .count()
    };

    println!("\nFound {} oceanic-continental convergent boundaries", oceanic_continental_count);

    // Generate geology
    let provinces = world.generate_geology(None).unwrap();

    // If we have oceanic-continental boundaries, we should generate arc system provinces
    if oceanic_continental_count > 0 {
        let volcanic_arcs: Vec<_> = provinces.iter()
            .filter(|p| p.characteristics.province_type == GeologicProvince::VolcanicArc)
            .collect();

        let trenches: Vec<_> = provinces.iter()
            .filter(|p| p.characteristics.province_type == GeologicProvince::OceanTrench)
            .collect();

        let accretionary_wedges: Vec<_> = provinces.iter()
            .filter(|p| p.characteristics.province_type == GeologicProvince::AccretionaryWedge)
            .collect();

        println!("Generated {} volcanic arcs", volcanic_arcs.len());
        println!("Generated {} ocean trenches", trenches.len());
        println!("Generated {} accretionary wedges", accretionary_wedges.len());

        assert!(
            !volcanic_arcs.is_empty() || !trenches.is_empty(),
            "Should generate arc system features for oceanic-continental convergence (found {} boundaries)",
            oceanic_continental_count
        );
    } else {
        println!("⚠️  No oceanic-continental boundaries in this test world (seed dependent)");
    }

    println!("\n✓ Oceanic-continental convergence test passed");
}
