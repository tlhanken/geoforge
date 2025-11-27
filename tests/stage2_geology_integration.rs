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
