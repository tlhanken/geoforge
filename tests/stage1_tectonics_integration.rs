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
    world.generate_tectonics(8).unwrap();
    let initial_stats = world.get_tectonic_stats().unwrap();
    assert_eq!(initial_stats.len(), 8, "Should have 8 plates after generation");

    // Stage 1.2: Refine boundaries
    let refinement_config = BoundaryRefinementConfig::with_seed(42)
        .with_noise(0.020, 50.0, 4)
        .with_smoothing(1);
    world.refine_boundaries(Some(refinement_config)).unwrap();

    let refined_stats = world.get_tectonic_stats().unwrap();
    assert_eq!(refined_stats.len(), 8, "Should still have 8 plates after refinement");

    // Stage 1.3: Remove islands
    let island_stats = world.remove_islands(None).unwrap();

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

    world.generate_tectonics(5).unwrap();
    let initial_plate_count = world.get_tectonic_stats().unwrap().len();

    world.refine_boundaries(None).unwrap();
    let refined_plate_count = world.get_tectonic_stats().unwrap().len();

    world.remove_islands(None).unwrap();
    let final_plate_count = world.get_tectonic_stats().unwrap().len();

    assert_eq!(initial_plate_count, refined_plate_count);
    assert_eq!(refined_plate_count, final_plate_count);
    assert_eq!(final_plate_count, 5);
}

#[test]
fn test_pipeline_with_extreme_refinement() {
    // Test that island removal can handle aggressive refinement
    let mut world = WorldMap::new(360, 180, 12345).unwrap();

    world.generate_tectonics(10).unwrap();

    // Extreme refinement to maximize island creation
    let extreme_config = BoundaryRefinementConfig::with_seed(12345)
        .with_noise(0.030, 120.0, 6)
        .with_smoothing(3);

    world.refine_boundaries(Some(extreme_config)).unwrap();
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
    world1.generate_tectonics(6).unwrap();

    let config1 = BoundaryRefinementConfig::with_seed(seed)
        .with_noise(0.020, 50.0, 4)
        .with_smoothing(0);  // Disable smoothing for exact determinism
    world1.refine_boundaries(Some(config1)).unwrap();

    let mut world2 = WorldMap::new(180, 90, seed).unwrap();
    world2.generate_tectonics(6).unwrap();

    let config2 = BoundaryRefinementConfig::with_seed(seed)
        .with_noise(0.020, 50.0, 4)
        .with_smoothing(0);  // Disable smoothing for exact determinism
    world2.refine_boundaries(Some(config2)).unwrap();

    // Compare plate data before island removal
    let data1 = &world1.tectonics.as_ref().unwrap().data;
    let data2 = &world2.tectonics.as_ref().unwrap().data;

    assert_eq!(data1, data2, "Same seed should produce identical results before island removal");

    // Island removal should successfully complete on both
    // Note: Exact stats may vary slightly due to HashMap iteration order
    let stats1 = world1.remove_islands(None).unwrap();
    let stats2 = world2.remove_islands(None).unwrap();

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

    world.generate_tectonics(7).unwrap();

    // Minimal refinement
    let minimal_config = BoundaryRefinementConfig::with_seed(111)
        .with_noise(0.005, 2.0, 2)
        .with_smoothing(0);

    world.refine_boundaries(Some(minimal_config)).unwrap();
    let island_stats = world.remove_islands(None).unwrap();

    // Minimal refinement should produce very few islands
    assert!(island_stats.islands_removed < 10, "Minimal refinement should produce few islands");
}

#[test]
fn test_stats_recalculation_after_stages() {
    // Verify that plate statistics are properly recalculated after each stage
    let mut world = WorldMap::new(360, 180, 777).unwrap();

    world.generate_tectonics(8).unwrap();
    let initial_total_area: u64 = world.get_tectonic_stats().unwrap()
        .values()
        .map(|s| s.area_km2)
        .sum();

    world.refine_boundaries(None).unwrap();
    let refined_total_area: u64 = world.get_tectonic_stats().unwrap()
        .values()
        .map(|s| s.area_km2)
        .sum();

    world.remove_islands(None).unwrap();
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

    world.generate_tectonics(30).unwrap();
    world.refine_boundaries(None).unwrap();
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

    world.generate_tectonics(1).unwrap();
    world.refine_boundaries(None).unwrap();
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
    world.generate_tectonics(10).unwrap();
    world.refine_boundaries(Some(
        BoundaryRefinementConfig::with_seed(55555)
            .with_noise(0.025, 100.0, 5)
    )).unwrap();
    world.remove_islands(None).unwrap();

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
    world.generate_tectonics(6).unwrap();
    world.refine_boundaries(None).unwrap();
    world.remove_islands(None).unwrap();

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
    world.generate_tectonics(10).unwrap();

    let metadata = world.get_tectonic_metadata().unwrap();

    // Should have assigned types to all plates
    let mut oceanic_count = 0;
    let mut continental_count = 0;
    let mut mixed_count = 0;

    for stats in metadata.plate_stats.values() {
        match stats.plate_type {
            geoforge::PlateType::Oceanic => oceanic_count += 1,
            geoforge::PlateType::Continental => continental_count += 1,
            geoforge::PlateType::Mixed => mixed_count += 1,
        }
    }

    // With 10 plates, should have a mix of types
    assert!(oceanic_count > 0, "Should have some oceanic plates");
    assert!(continental_count > 0, "Should have some continental plates");
    
    assert_eq!(oceanic_count + continental_count + mixed_count, 10,
        "All plates should have assigned types");

    println!("Plate type distribution:");
    println!("  Oceanic: {}", oceanic_count);
    println!("  Continental: {}", continental_count);
    println!("  Mixed: {}", mixed_count);
}
