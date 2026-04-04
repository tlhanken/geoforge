//! Island removal and plate contiguity enforcement
//!
//! This module provides post-processing to ensure all tectonic plates are contiguous.
//! After boundary refinement, plates may develop isolated "islands" - small disconnected
//! fragments separated from the main plate body. This module detects such fragments and
//! reassigns them to surrounding plates.

use crate::map::terrain::TerrainMap;
use std::collections::{HashMap, HashSet, VecDeque};

/// Configuration for island removal
#[derive(Debug, Clone)]
pub struct IslandRemovalConfig {
    /// Minimum size (in pixels) for a plate fragment to be considered the "main" body
    /// Fragments smaller than this are always considered islands
    pub min_main_body_size: usize,

    /// Verbose output for debugging
    pub verbose: bool,
}

impl Default for IslandRemovalConfig {
    fn default() -> Self {
        Self {
            min_main_body_size: 100,
            verbose: true,
        }
    }
}

impl IslandRemovalConfig {
    /// Create a new configuration with verbose output enabled
    pub fn new() -> Self {
        Self::default()
    }

    /// Set minimum main body size
    pub fn with_min_body_size(mut self, size: usize) -> Self {
        self.min_main_body_size = size;
        self
    }

    /// Set verbose output
    pub fn with_verbose(mut self, verbose: bool) -> Self {
        self.verbose = verbose;
        self
    }
}

/// Island removal processor
pub struct IslandRemover {
    config: IslandRemovalConfig,
}

impl IslandRemover {
    /// Create a new island remover with configuration
    pub fn new(config: IslandRemovalConfig) -> Self {
        Self { config }
    }

    /// Remove islands from the plate map, ensuring all plates are contiguous
    pub fn remove_islands(&mut self, plate_map: &mut TerrainMap<u16>) -> IslandRemovalStats {
        if self.config.verbose {
            println!("Removing plate islands to ensure contiguity...");
        }

        let mut stats = IslandRemovalStats::default();

        // Get all unique plate IDs
        let mut plate_ids: HashSet<u16> = plate_map.data.iter().copied().collect();
        plate_ids.remove(&0); // Remove "no plate" if present

        if self.config.verbose {
            println!("  Analyzing {} plates...", plate_ids.len());
        }

        // For each plate, find all connected components
        for plate_id in plate_ids {
            let components = self.find_connected_components(plate_map, plate_id);

            if components.len() > 1 {
                stats.plates_with_islands += 1;

                if self.config.verbose {
                    println!("  Plate {} has {} fragments", plate_id, components.len());
                }

                // Find the largest component (main body)
                let main_component_idx = components
                    .iter()
                    .enumerate()
                    .max_by_key(|(_, component)| component.len())
                    .map(|(idx, _)| idx)
                    .unwrap();

                let main_body_size = components[main_component_idx].len();

                // Reassign all smaller components (islands)
                for (idx, component) in components.iter().enumerate() {
                    if idx != main_component_idx {
                        stats.islands_removed += 1;
                        stats.pixels_reassigned += component.len();

                        if self.config.verbose {
                            println!("    Removing island of {} pixels ({:.1}% of main body)",
                                component.len(),
                                (component.len() as f64 / main_body_size as f64) * 100.0);
                        }

                        // Find the most common neighboring plate and reassign to it
                        let new_plate_id = self.find_surrounding_plate(plate_map, component, plate_id);

                        for &(x, y) in component {
                            let idx = plate_map.get_index(x, y);
                            plate_map.data[idx] = new_plate_id;
                        }
                    }
                }
            }
        }

        if self.config.verbose {
            println!("Island removal complete!");
            println!("  {} plates had islands", stats.plates_with_islands);
            println!("  {} islands removed", stats.islands_removed);
            println!("  {} pixels reassigned", stats.pixels_reassigned);
        }

        stats
    }

    /// Find all connected components for a given plate ID using flood fill
    fn find_connected_components(
        &self,
        plate_map: &TerrainMap<u16>,
        plate_id: u16,
    ) -> Vec<Vec<(usize, usize)>> {
        let mut visited = HashSet::new();
        let mut components = Vec::new();

        // Find all pixels belonging to this plate
        for y in 0..plate_map.height {
            for x in 0..plate_map.width {
                if visited.contains(&(x, y)) {
                    continue;
                }

                let idx = plate_map.get_index(x, y);
                if plate_map.data[idx] == plate_id {
                    // Start a new connected component with flood fill
                    let component = self.flood_fill(plate_map, x, y, plate_id, &mut visited);
                    if !component.is_empty() {
                        components.push(component);
                    }
                }
            }
        }

        components
    }

    /// Perform flood fill to find a connected component
    fn flood_fill(
        &self,
        plate_map: &TerrainMap<u16>,
        start_x: usize,
        start_y: usize,
        plate_id: u16,
        visited: &mut HashSet<(usize, usize)>,
    ) -> Vec<(usize, usize)> {
        let mut component = Vec::new();
        let mut queue = VecDeque::new();

        queue.push_back((start_x, start_y));
        visited.insert((start_x, start_y));

        while let Some((x, y)) = queue.pop_front() {
            component.push((x, y));

            // Check all neighbors
            for (nx, ny) in plate_map.get_neighbors(x, y) {
                if visited.contains(&(nx, ny)) {
                    continue;
                }

                let idx = plate_map.get_index(nx, ny);
                if plate_map.data[idx] == plate_id {
                    visited.insert((nx, ny));
                    queue.push_back((nx, ny));
                }
            }
        }

        component
    }

    /// Find the plate that surrounds an island by examining neighbors
    fn find_surrounding_plate(
        &self,
        plate_map: &TerrainMap<u16>,
        island: &[(usize, usize)],
        original_plate_id: u16,
    ) -> u16 {
        let mut neighbor_counts: HashMap<u16, usize> = HashMap::new();

        // Count neighboring plates for all pixels in the island
        for &(x, y) in island {
            for (nx, ny) in plate_map.get_neighbors(x, y) {
                let idx = plate_map.get_index(nx, ny);
                let neighbor_plate = plate_map.data[idx];

                // Don't count the original plate or "no plate"
                if neighbor_plate != original_plate_id && neighbor_plate != 0 {
                    *neighbor_counts.entry(neighbor_plate).or_insert(0) += 1;
                }
            }
        }

        // Return the most common neighboring plate
        neighbor_counts
            .into_iter()
            .max_by_key(|(_, count)| *count)
            .map(|(plate_id, _)| plate_id)
            .unwrap_or(original_plate_id) // Fallback to original if no neighbors found
    }
}

/// Statistics about island removal
#[derive(Debug, Clone, Default)]
pub struct IslandRemovalStats {
    /// Number of plates that had islands
    pub plates_with_islands: usize,

    /// Total number of islands removed
    pub islands_removed: usize,

    /// Total number of pixels reassigned
    pub pixels_reassigned: usize,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_island_removal() {
        // Create a simple map with an island
        // Plate 1: main body on left side (x=0-2)
        // Plate 2: fills most of the map
        // Plate 1 island: small fragment in the middle (x=5, y=4-5), completely surrounded by plate 2
        let mut plate_map = TerrainMap::new(10, 10, 2u16);

        // Make left section plate 1 (main body)
        for y in 0..10 {
            for x in 0..3 {
                let idx = plate_map.get_index(x, y);
                plate_map.data[idx] = 1;
            }
        }

        // Make small island of plate 1 in the middle, surrounded by plate 2
        // This is definitely an island because plate 2 separates it from the main body
        for y in 4..6 {
            let idx = plate_map.get_index(5, y);
            plate_map.data[idx] = 1;
        }

        let config = IslandRemovalConfig::default().with_verbose(false);
        let mut remover = IslandRemover::new(config);
        let stats = remover.remove_islands(&mut plate_map);

        // Verify the island was removed
        assert_eq!(stats.islands_removed, 1, "Expected 1 island to be removed");
        assert_eq!(stats.pixels_reassigned, 2, "Expected 2 pixels to be reassigned");

        // Verify the island was reassigned to plate 2 (the surrounding plate)
        for y in 4..6 {
            let idx = plate_map.get_index(5, y);
            assert_eq!(plate_map.data[idx], 2, "Island pixel at (5, {}) should be reassigned to plate 2", y);
        }
    }

    #[test]
    fn test_no_islands() {
        // Create a map with no islands
        let mut plate_map = TerrainMap::new(10, 10, 1u16);

        // Make right half plate 2
        for y in 0..10 {
            for x in 5..10 {
                let idx = plate_map.get_index(x, y);
                plate_map.data[idx] = 2;
            }
        }

        let config = IslandRemovalConfig::default().with_verbose(false);
        let mut remover = IslandRemover::new(config);
        let stats = remover.remove_islands(&mut plate_map);

        // No islands should be found
        assert_eq!(stats.islands_removed, 0);
        assert_eq!(stats.pixels_reassigned, 0);
    }

    #[test]
    fn test_multiple_islands() {
        // Create a map with multiple islands for the same plate
        let mut plate_map = TerrainMap::new(20, 20, 2u16);

        // Plate 1: main body in center
        for y in 5..15 {
            for x in 5..15 {
                let idx = plate_map.get_index(x, y);
                plate_map.data[idx] = 1;
            }
        }

        // Plate 1: island in top-left corner
        for y in 0..2 {
            for x in 0..2 {
                let idx = plate_map.get_index(x, y);
                plate_map.data[idx] = 1;
            }
        }

        // Plate 1: island in bottom-right corner
        for y in 18..20 {
            for x in 18..20 {
                let idx = plate_map.get_index(x, y);
                plate_map.data[idx] = 1;
            }
        }

        let config = IslandRemovalConfig::default().with_verbose(false);
        let mut remover = IslandRemover::new(config);
        let stats = remover.remove_islands(&mut plate_map);

        // Should find 2 islands
        assert_eq!(stats.islands_removed, 2);
        assert_eq!(stats.plates_with_islands, 1);
    }
}
