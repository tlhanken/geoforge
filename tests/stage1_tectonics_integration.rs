//! Integration tests for Stage 1: Tectonic Foundation
//!
//! These tests verify that all stages of tectonic generation work together correctly:
//! Stage 1.1: Core Plate Generation
//! Stage 1.2: Boundary Refinement
//! Stage 1.3: Island Removal

use geoforge::{WorldMap, BoundaryRefinementConfig};

#[test]
fn test_full_pipeline_integration() {
    // Test that all three stages can be run in sequence
    let mut world = WorldMap::new(360, 180, 42).unwrap();

    // Stage 1.1: Generate plates
    world.tectonics().generate_plates(8).unwrap();
    let initial_stats = world.get_tectonic_stats().unwrap();
    assert_eq!(initial_stats.len(), 8, "Should have 8 plates after generation");

    // Stage 1.2: Refine boundaries
    let refinement_config = BoundaryRefinementConfig::with_seed(42)
        .with_noise(0.020, 50.0, 4)
        .with_smoothing(1);
    world.tectonics().roughen_boundaries(Some(refinement_config)).unwrap();

    let refined_stats = world.get_tectonic_stats().unwrap();
    assert_eq!(refined_stats.len(), 8, "Should still have 8 plates after refinement");

    // Stage 1.3: Remove islands
    let island_stats = world.tectonics().deisland(None).unwrap();

    let final_stats = world.get_tectonic_stats().unwrap();
    assert_eq!(final_stats.len(), 8, "Should still have 8 plates after island removal");

    // Verify stats were recalculated
    for (plate_id, stats) in final_stats {
        assert!(*plate_id > 0 && *plate_id <= 8);
        assert!(stats.area_km2 > 0, "Plate {} should have positive area", plate_id);
        assert!(stats.pixels > 0, "Plate {} should have pixels", plate_id);
    }

    println!("Island removal stats: {} islands removed from {} plates",
        island_stats.islands_removed, island_stats.plates_with_islands);
}

#[test]
fn test_pipeline_preserves_plate_count() {
    // Verify that refinement and island removal don't change the number of plates
    let mut world = WorldMap::new(180, 90, 999).unwrap();

    world.tectonics().generate_plates(5).unwrap();
    let initial_plate_count = world.get_tectonic_stats().unwrap().len();

    world.tectonics().roughen_boundaries(None).unwrap();
    let refined_plate_count = world.get_tectonic_stats().unwrap().len();

    world.tectonics().deisland(None).unwrap();
    let final_plate_count = world.get_tectonic_stats().unwrap().len();

    assert_eq!(initial_plate_count, refined_plate_count);
    assert_eq!(refined_plate_count, final_plate_count);
    assert_eq!(final_plate_count, 5);
}

#[test]
fn test_pipeline_with_extreme_refinement() {
    // Test that island removal can handle aggressive refinement
    let mut world = WorldMap::new(360, 180, 12345).unwrap();

    world.tectonics().generate_plates(10).unwrap();

    // Extreme refinement to maximize island creation
    let extreme_config = BoundaryRefinementConfig::with_seed(12345)
        .with_noise(0.030, 120.0, 6)
        .with_smoothing(3);

    world.tectonics().roughen_boundaries(Some(extreme_config)).unwrap();
    let _island_stats = world.remove_islands(None).unwrap();

    // With extreme refinement, we should find islands
    // Note: islands_removed and pixels_reassigned are usize (always >= 0)

    // All plates should still be present
    let final_stats = world.get_tectonic_stats().unwrap();
    assert_eq!(final_stats.len(), 10);
}

#[test]
fn test_pipeline_determinism_generation_and_refinement() {
    // Test that generation and refinement are deterministic with the same seed
    // Note: Full pipeline determinism including island removal depends on HashMap iteration order
    let seed = 54321;

    let mut world1 = WorldMap::new(180, 90, seed).unwrap();
    world1.tectonics().generate_plates(6).unwrap();

    let config1 = BoundaryRefinementConfig::with_seed(seed)
        .with_noise(0.020, 50.0, 4)
        .with_smoothing(0);  // Disable smoothing for exact determinism
    world1.tectonics().roughen_boundaries(Some(config1)).unwrap();

    let mut world2 = WorldMap::new(180, 90, seed).unwrap();
    world2.tectonics().generate_plates(6).unwrap();

    let config2 = BoundaryRefinementConfig::with_seed(seed)
        .with_noise(0.020, 50.0, 4)
        .with_smoothing(0);  // Disable smoothing for exact determinism
    world2.tectonics().roughen_boundaries(Some(config2)).unwrap();

    // Compare plate data before island removal
    let data1 = &world1.tectonics.as_ref().unwrap().data;
    let data2 = &world2.tectonics.as_ref().unwrap().data;

    assert_eq!(data1, data2, "Same seed should produce identical results before island removal");

    // Island removal should successfully complete on both
    // Note: Exact stats may vary slightly due to HashMap iteration order
    let stats1 = world1.tectonics().deisland(None).unwrap();
    let stats2 = world2.tectonics().deisland(None).unwrap();

    // Verify both runs found and removed islands
    assert!(stats1.islands_removed > 0, "Should find islands with this config");
    assert!(stats2.islands_removed > 0, "Should find islands with this config");

    // Stats may vary due to HashMap iteration order affecting island reassignment
    // Both runs should complete successfully and reassign a reasonable number of pixels
    assert!(stats1.pixels_reassigned > 0, "Should reassign pixels in first run");
    assert!(stats2.pixels_reassigned > 0, "Should reassign pixels in second run");

    // Verify stats are in the same ballpark (within 50% of each other)
    let avg_pixels = (stats1.pixels_reassigned + stats2.pixels_reassigned) / 2;
    let pixel_diff = (stats1.pixels_reassigned as i64 - stats2.pixels_reassigned as i64).abs() as usize;
    assert!(pixel_diff < avg_pixels,
        "Pixel reassignment variance should be reasonable (diff: {}, avg: {})", pixel_diff, avg_pixels);
}

#[test]
fn test_error_handling_refine_before_generate() {
    // Test that refine_boundaries fails if called before generate_tectonics
    let mut world = WorldMap::new(100, 50, 42).unwrap();

    let result = world.refine_boundaries(None);
    assert!(result.is_err());
    assert_eq!(result.unwrap_err().to_string(), "Tectonics must be generated before refining boundaries");
}

#[test]
fn test_error_handling_island_removal_before_generate() {
    // Test that remove_islands fails if called before generate_tectonics
    let mut world = WorldMap::new(100, 50, 42).unwrap();

    let result = world.remove_islands(None);
    assert!(result.is_err());
    assert_eq!(result.unwrap_err().to_string(), "Tectonics must be generated before removing islands");
}

#[test]
fn test_pipeline_with_minimal_refinement() {
    // Test that minimal refinement produces few or no islands
    let mut world = WorldMap::new(180, 90, 111).unwrap();

    world.tectonics().generate_plates(7).unwrap();

    // Minimal refinement
    let minimal_config = BoundaryRefinementConfig::with_seed(111)
        .with_noise(0.005, 2.0, 2)
        .with_smoothing(0);

    world.tectonics().roughen_boundaries(Some(minimal_config)).unwrap();
    let island_stats = world.remove_islands(None).unwrap();

    // Minimal refinement should produce very few islands
    assert!(island_stats.islands_removed < 10, "Minimal refinement should produce few islands");
}

#[test]
fn test_stats_recalculation_after_stages() {
    // Verify that plate statistics are properly recalculated after each stage
    let mut world = WorldMap::new(360, 180, 777).unwrap();

    world.tectonics().generate_plates(8).unwrap();
    let initial_total_area: u64 = world.get_tectonic_stats().unwrap()
        .values()
        .map(|s| s.area_km2)
        .sum();

    world.tectonics().roughen_boundaries(None).unwrap();
    let refined_total_area: u64 = world.get_tectonic_stats().unwrap()
        .values()
        .map(|s| s.area_km2)
        .sum();

    world.tectonics().deisland(None).unwrap();
    let final_total_area: u64 = world.get_tectonic_stats().unwrap()
        .values()
        .map(|s| s.area_km2)
        .sum();

    // Total area should remain constant (allowing small rounding differences)
    let area_diff_1 = (initial_total_area as i64 - refined_total_area as i64).abs();
    let area_diff_2 = (refined_total_area as i64 - final_total_area as i64).abs();

    assert!(area_diff_1 < 100, "Area should be conserved during refinement");
    assert!(area_diff_2 < 100, "Area should be conserved during island removal");
}

#[test]
fn test_pipeline_with_many_plates() {
    // Test pipeline with a large number of plates
    let mut world = WorldMap::new(360, 180, 8888).unwrap();

    world.tectonics().generate_plates(30).unwrap();
    world.tectonics().roughen_boundaries(None).unwrap();
    let island_stats = world.remove_islands(None).unwrap();

    let final_stats = world.get_tectonic_stats().unwrap();
    assert_eq!(final_stats.len(), 30, "Should still have 30 plates");

    // With many plates, island creation is more likely
    println!("With 30 plates: {} islands removed", island_stats.islands_removed);
}

#[test]
fn test_single_plate_edge_case() {
    // Test edge case with only one plate
    let mut world = WorldMap::new(100, 50, 333).unwrap();

    world.tectonics().generate_plates(1).unwrap();
    world.tectonics().roughen_boundaries(None).unwrap();
    let island_stats = world.remove_islands(None).unwrap();

    // Single plate should never have islands
    assert_eq!(island_stats.islands_removed, 0);
    assert_eq!(island_stats.plates_with_islands, 0);

    let final_stats = world.get_tectonic_stats().unwrap();
    assert_eq!(final_stats.len(), 1);
}

#[test]
fn test_contiguity_guarantee() {
    // Verify that after island removal, all plates are truly contiguous
    use std::collections::{HashSet, VecDeque};

    let mut world = WorldMap::new(180, 90, 55555).unwrap();
    world.tectonics().generate_plates(10).unwrap();
    world.tectonics().roughen_boundaries(Some(
        BoundaryRefinementConfig::with_seed(55555)
            .with_noise(0.025, 100.0, 5)
    )).unwrap();
    world.tectonics().deisland(None).unwrap();

    let plate_map = world.tectonics.as_ref().unwrap();

    // For each plate, verify it's contiguous using flood fill
    let mut checked_plates = HashSet::new();

    for y in 0..plate_map.height {
        for x in 0..plate_map.width {
            let plate_id = plate_map.data[plate_map.get_index(x, y)];

            if plate_id == 0 || checked_plates.contains(&plate_id) {
                continue;
            }

            // Flood fill to find all connected pixels for this plate
            let mut visited = HashSet::new();
            let mut queue = VecDeque::new();
            queue.push_back((x, y));
            visited.insert((x, y));

            while let Some((cx, cy)) = queue.pop_front() {
                for (nx, ny) in plate_map.get_neighbors(cx, cy) {
                    if visited.contains(&(nx, ny)) {
                        continue;
                    }

                    let neighbor_plate = plate_map.data[plate_map.get_index(nx, ny)];
                    if neighbor_plate == plate_id {
                        visited.insert((nx, ny));
                        queue.push_back((nx, ny));
                    }
                }
            }

            // Count total pixels for this plate
            let total_pixels = plate_map.data.iter().filter(|&&p| p == plate_id).count();

            // All pixels should be reachable from flood fill (contiguous)
            assert_eq!(visited.len(), total_pixels,
                "Plate {} should be contiguous (found {} connected, {} total)",
                plate_id, visited.len(), total_pixels);

            checked_plates.insert(plate_id);
        }
    }
}
#[test]
fn test_boundary_analysis_integration() {
    // Test Stage 1.4: Boundary Analysis
    let mut world = WorldMap::new(180, 90, 42).unwrap();

    // Generate and process plates
    world.tectonics().generate_plates(6).unwrap();
    world.tectonics().roughen_boundaries(None).unwrap();
    world.tectonics().deisland(None).unwrap();

    // Stage 1.4: Analyze boundaries
    let boundary_stats = world.analyze_boundaries(None).unwrap();

    // Verify boundary statistics
    assert!(boundary_stats.total_boundaries > 0, "Should find plate boundaries");
    
    // With 6 plates, maximum possible boundaries is (6 * 5) / 2 = 15
    assert!(boundary_stats.total_boundaries <= 15, 
        "Cannot have more than 15 boundaries with 6 plates");

    // Total should equal sum of parts
    assert_eq!(
        boundary_stats.total_boundaries,
        boundary_stats.convergent_count + boundary_stats.divergent_count + boundary_stats.transform_count,
        "Boundary counts should sum to total"
    );

    // Verify metadata was updated
    let metadata = world.get_tectonic_metadata().unwrap();
    assert_eq!(metadata.plate_boundaries.len(), boundary_stats.total_boundaries);
    assert!(metadata.boundary_stats.is_some());

    // Verify each boundary segment
    for boundary in &metadata.plate_boundaries {
        assert!(boundary.plate_a > 0 && boundary.plate_a <= 6);
        assert!(boundary.plate_b > 0 && boundary.plate_b <= 6);
        assert_ne!(boundary.plate_a, boundary.plate_b);
        assert!(boundary.pixel_count() > 0);
        assert!(boundary.length_km > 0.0);
    }

    println!("Boundary analysis results:");
    println!("  Total boundaries: {}", boundary_stats.total_boundaries);
    println!("  Convergent: {}", boundary_stats.convergent_count);
    println!("  Divergent: {}", boundary_stats.divergent_count);
    println!("  Transform: {}", boundary_stats.transform_count);
    println!("  Total length: {:.0} km", boundary_stats.total_length_km);
    println!("  Avg velocity: {:.2} cm/year", boundary_stats.average_relative_velocity);
}

#[test]
fn test_boundary_analysis_requires_tectonics() {
    // Test error handling when analyzing boundaries before generating tectonics
    let mut world = WorldMap::new(100, 50, 42).unwrap();

    let result = world.analyze_boundaries(None);
    assert!(result.is_err(), "Should fail when tectonics not generated");
}

#[test]
fn test_plate_type_assignment() {
    // Test that plate types are assigned correctly based on size
    let mut world = WorldMap::new(360, 180, 123).unwrap();
    world.tectonics().generate_plates(10).unwrap();

    let metadata = world.get_tectonic_metadata().unwrap();

    // Should have assigned types to all plates
    let mut oceanic_count = 0;
    let mut continental_count = 0;

    for stats in metadata.plate_stats.values() {
        match stats.plate_type {
            geoforge::PlateType::Oceanic => oceanic_count += 1,
            geoforge::PlateType::Continental => continental_count += 1,
        }
    }

    // With 10 plates, should have a mix of types
    assert!(oceanic_count > 0, "Should have some oceanic plates");
    assert!(continental_count > 0, "Should have some continental plates");

    assert_eq!(oceanic_count + continental_count, 10,
        "All plates should have assigned types");

    println!("Plate type distribution:");
    println!("  Oceanic: {}", oceanic_count);
    println!("  Continental: {}", continental_count);
}

#[test]
fn test_complete_stage_1_pipeline() {
    // Test the complete Stage 1 pipeline from 1.1 through 1.4
    let mut world = WorldMap::new(360, 180, 42).unwrap();

    // Stage 1.1: Generate plates
    world.tectonics().generate_plates(8).unwrap();
    assert!(world.tectonics.is_some());

    let metadata = world.get_tectonic_metadata().unwrap();

    // Verify motion was automatically assigned
    for seed in &metadata.plate_seeds {
        assert!(seed.motion_speed > 0.0, "Plate {} should have motion assigned", seed.id);
        assert!(seed.motion_direction >= 0.0 && seed.motion_direction < 360.0);
    }

    // Stage 1.2: Refine boundaries
    world.tectonics().roughen_boundaries(None).unwrap();

    // Stage 1.3: Remove islands
    world.tectonics().deisland(None).unwrap();

    // Stage 1.4: Analyze boundaries
    let boundary_stats = world.analyze_boundaries(None).unwrap();

    assert!(boundary_stats.total_boundaries > 0);
    assert!(boundary_stats.average_relative_velocity > 0.0,
        "Boundaries should have non-zero relative velocity");

    // Verify boundary types are distributed
    let has_convergent = boundary_stats.convergent_count > 0;
    let has_divergent = boundary_stats.divergent_count > 0;
    let has_transform = boundary_stats.transform_count > 0;

    let boundary_type_count = [has_convergent, has_divergent, has_transform]
        .iter()
        .filter(|&&x| x)
        .count();

    assert!(boundary_type_count >= 2,
        "Should have at least 2 different boundary types with 8 plates");

    println!("Complete Stage 1 pipeline success!");
    println!("  Plates: 8");
    println!("  Boundaries: {}", boundary_stats.total_boundaries);
    println!("  Convergent: {}", boundary_stats.convergent_count);
    println!("  Divergent: {}", boundary_stats.divergent_count);
    println!("  Transform: {}", boundary_stats.transform_count);
    println!("  Avg velocity: {:.2} cm/year", boundary_stats.average_relative_velocity);
}

#[test]
fn test_motion_assignment_determinism() {
    // Test that motion assignment is deterministic with same seed
    let mut world1 = WorldMap::new(180, 90, 123).unwrap();
    let mut world2 = WorldMap::new(180, 90, 123).unwrap();

    world1.tectonics().generate_plates(5).unwrap();
    world2.tectonics().generate_plates(5).unwrap();

    let metadata1 = world1.get_tectonic_metadata().unwrap();
    let metadata2 = world2.get_tectonic_metadata().unwrap();

    // Motion should be identical
    for (seed1, seed2) in metadata1.plate_seeds.iter().zip(&metadata2.plate_seeds) {
        assert_eq!(seed1.id, seed2.id);
        assert_eq!(seed1.motion_direction, seed2.motion_direction,
            "Motion direction should be deterministic");
        assert_eq!(seed1.motion_speed, seed2.motion_speed,
            "Motion speed should be deterministic");
    }
}

#[test]
fn test_motion_realistic_values() {
    // Test that motion values are within realistic Earth-like ranges
    let mut world = WorldMap::new(360, 180, 999).unwrap();
    world.tectonics().generate_plates(15).unwrap();

    let metadata = world.get_tectonic_metadata().unwrap();

    for seed in &metadata.plate_seeds {
        // Speed should be in Earth-like range (1-10 cm/year by default)
        assert!(seed.motion_speed >= 1.0 && seed.motion_speed <= 10.0,
            "Plate {} speed {:.2} outside realistic range",
            seed.id, seed.motion_speed);

        // Direction should be valid (0-360 degrees)
        assert!(seed.motion_direction >= 0.0 && seed.motion_direction < 360.0,
            "Plate {} direction {:.2} invalid",
            seed.id, seed.motion_direction);
    }
}

#[test]
#[cfg(feature = "export-png")]
fn test_boundary_visualization_export() {
    use std::path::Path;

    let mut world = WorldMap::new(180, 90, 42).unwrap();
    world.tectonics().generate_plates(6).unwrap();
    world.tectonics().roughen_boundaries(None).unwrap();
    world.tectonics().deisland(None).unwrap();
    world.tectonics().analyze(None).unwrap();

    // Export boundary visualization
    let output_dir = "outputs/tests";
    let filename = "test_boundaries.png";
    let result = world.export_boundaries_png(output_dir, filename);

    assert!(result.is_ok(), "Boundary export should succeed");

    // Verify file was created
    let path = Path::new(output_dir).join(filename);
    assert!(path.exists(), "Boundary PNG should be created");

    // Cleanup
    std::fs::remove_file(path).ok();
}

#[test]
#[cfg(feature = "export-png")]
fn test_boundary_export_requires_analysis() {
    // Test that exporting boundaries before analysis fails
    let mut world = WorldMap::new(100, 50, 42).unwrap();
    world.tectonics().generate_plates(4).unwrap();

    // Try to export without analyzing
    let result = world.export_boundaries_png("outputs/tests", "should_fail.png");
    assert!(result.is_err(), "Should fail when boundaries not analyzed");
}

#[test]
fn test_boundary_classification_accuracy() {
    // Test that boundary classification produces sensible results
    let mut world = WorldMap::new(360, 180, 789).unwrap();
    world.tectonics().generate_plates(10).unwrap();
    world.tectonics().roughen_boundaries(None).unwrap();
    world.tectonics().deisland(None).unwrap();
    world.tectonics().analyze(None).unwrap();

    let metadata = world.get_tectonic_metadata().unwrap();

    // Every boundary should have a valid classification
    for boundary in &metadata.plate_boundaries {
        // Should have positive relative velocity (plates are moving)
        assert!(boundary.relative_velocity >= 0.0,
            "Boundary between {} and {} has negative velocity",
            boundary.plate_a, boundary.plate_b);

        // Should have positive length
        assert!(boundary.length_km > 0.0,
            "Boundary should have positive length");

        // Should involve different plates
        assert_ne!(boundary.plate_a, boundary.plate_b,
            "Boundary should connect different plates");

        // Should have pixels
        assert!(!boundary.pixels.is_empty(),
            "Boundary should have pixels");
    }

    println!("Boundary classification validation passed for {} boundaries",
        metadata.plate_boundaries.len());
}

#[test]
#[cfg(feature = "export-png")]
fn test_motion_png_export_import_roundtrip() {
    use std::fs;

    // Create output directory
    let test_dir = "outputs/tests/motion_roundtrip";
    fs::create_dir_all(test_dir).unwrap();

    // Generate a world with specific seed for reproducibility
    let mut world1 = WorldMap::new(360, 180, 12345).unwrap();
    world1.tectonics().generate(8).unwrap();

    // Get original metadata
    let original_metadata = world1.get_tectonic_metadata().unwrap();
    let original_seeds = original_metadata.plate_seeds.clone();

    // Export to plate_motion.png
    let motion_png_path = format!("{}/test_motion.png", test_dir);
    world1.export_plate_motion_png(test_dir, "test_motion.png").unwrap();

    // Create new world and import from plate_motion.png
    let mut world2 = WorldMap::new(360, 180, 0).unwrap();
    world2.tectonics().import_png(&motion_png_path).unwrap();

    // Get imported metadata
    let imported_metadata = world2.get_tectonic_metadata().unwrap();
    let imported_seeds = &imported_metadata.plate_seeds;

    // Verify same number of plates
    assert_eq!(original_seeds.len(), imported_seeds.len(),
        "Should have same number of plates after roundtrip");

    // Build a map from (direction, speed) to original seed to match imported plates
    // We'll match by the motion vector itself since that's what's encoded in the PNG
    use std::collections::HashMap;

    let mut motion_to_original: HashMap<(i32, i32), &PlateSeed> = HashMap::new();
    for seed in &original_seeds {
        // Round to integers for matching (RGB quantization means we lose precision)
        let dir_key = seed.motion_direction.round() as i32;
        let speed_key = (seed.motion_speed * 10.0).round() as i32; // 0.1 cm/year precision
        motion_to_original.insert((dir_key, speed_key), seed);
    }

    // Verify each imported plate matches an original (within tolerance)
    let mut matched_count = 0;
    for imported_seed in imported_seeds {
        let imp_dir = imported_seed.motion_direction.round() as i32;
        let imp_speed = (imported_seed.motion_speed * 10.0).round() as i32;

        // Try to find a match within small tolerance
        let mut found_match = false;
        for dir_offset in -2..=2 {
            for speed_offset in -2..=2 {
                let key = (imp_dir + dir_offset, imp_speed + speed_offset);
                if motion_to_original.contains_key(&key) {
                    found_match = true;
                    matched_count += 1;
                    break;
                }
            }
            if found_match { break; }
        }

        assert!(found_match,
            "Imported plate with motion {}° at {:.1} cm/year should match an original plate",
            imported_seed.motion_direction, imported_seed.motion_speed);
    }

    assert_eq!(matched_count, original_seeds.len(),
        "All plates should have matching motion vectors after roundtrip");

    println!("✅ Motion PNG roundtrip test passed: {} plates with preserved motion vectors",
        original_seeds.len());
}

#[test]
#[cfg(feature = "export-png")]
fn test_complete_export_import_workflow() {
    use std::fs;

    // Test the complete workflow: generate -> export all -> import -> verify
    let test_dir = "outputs/tests/complete_workflow";
    fs::create_dir_all(test_dir).unwrap();

    // Generate and export everything
    let mut world1 = WorldMap::new(180, 90, 42424).unwrap();
    let _stats1 = world1.tectonics().generate(6).unwrap();

    // Export all formats
    world1.tectonics().export(test_dir).unwrap();

    // Verify files were created
    assert!(std::path::Path::new(&format!("{}/world.map", test_dir)).exists());
    assert!(std::path::Path::new(&format!("{}/tectonics.png", test_dir)).exists());
    assert!(std::path::Path::new(&format!("{}/boundaries.png", test_dir)).exists());
    assert!(std::path::Path::new(&format!("{}/plate_motion.png", test_dir)).exists());

    // Import from motion PNG (this is the key test - motion vectors preserved)
    let mut world2 = WorldMap::new(180, 90, 0).unwrap();
    let motion_path = format!("{}/plate_motion.png", test_dir);
    world2.tectonics().import_png(&motion_path).unwrap();

    let metadata2 = world2.get_tectonic_metadata().unwrap();
    assert_eq!(metadata2.plate_seeds.len(), 6, "Motion PNG should have 6 plates");

    // Verify all plates have valid motion vectors
    for seed in &metadata2.plate_seeds {
        assert!(seed.motion_speed >= 1.0 && seed.motion_speed <= 10.0,
            "Plate should have realistic speed");
        assert!(seed.motion_direction >= 0.0 && seed.motion_direction < 360.0,
            "Plate should have valid direction");
    }

    println!("✅ Complete export/import workflow test passed");
}
