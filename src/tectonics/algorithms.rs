//! Tectonic plate generation algorithms

use crate::map::TerrainMap;
use crate::spherical::SphericalPoint;
use crate::tectonics::plates::PlateSeed;
use crate::tectonics::PlateError;
use rand::prelude::*;
use rand::rngs::StdRng;
use std::collections::{HashMap, HashSet, BinaryHeap};
use std::cmp::Ordering;

const MIN_SEED_DISTANCE: f64 = 100.0; // km
const MAX_SEED_ATTEMPTS: usize = 10000;

/// Pixel-based growth frontier item for region growing
#[derive(Clone)]
pub struct GrowthPixel {
    pub x: usize,
    pub y: usize,
    pub plate_id: u16,
    pub distance: f64,  // Geodesic distance from seed
}

impl PartialEq for GrowthPixel {
    fn eq(&self, other: &Self) -> bool {
        self.distance == other.distance
    }
}

impl Eq for GrowthPixel {}

impl PartialOrd for GrowthPixel {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        // Reverse ordering for min-heap behavior
        other.distance.partial_cmp(&self.distance)
    }
}

impl Ord for GrowthPixel {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap_or(Ordering::Equal)
    }
}

/// Generate random seed points for plates with improved distribution
pub fn generate_seeds(
    width: usize,
    height: usize,
    num_plates: usize,
    rng: &mut StdRng,
    plate_size_targets: &[f64],
) -> Result<Vec<PlateSeed>, PlateError> {
    let mut plate_seeds = Vec::with_capacity(num_plates);
    
    for i in 0..num_plates {
        let mut attempts = 0;
        let (x, y, point_3d) = loop {
            attempts += 1;
            
            if attempts > MAX_SEED_ATTEMPTS {
                return Err(PlateError::SeedPlacementFailed);
            }

            // Generate random point on sphere with uniform distribution
            let u = rng.gen::<f64>();
            let v = rng.gen::<f64>();
            
            let theta = 2.0 * std::f64::consts::PI * u;
            let phi = (2.0 * v - 1.0).acos();
            
            let x_3d = phi.sin() * theta.cos();
            let y_3d = phi.sin() * theta.sin();
            let z_3d = phi.cos();
            
            let point_3d = SphericalPoint::from_lat_lon(
                z_3d.asin().to_degrees(),
                y_3d.atan2(x_3d).to_degrees()
            );
            let (lat, lon) = point_3d.to_lat_lon();
            
            // Convert to pixel coordinates
            let x = ((lon + 180.0) / 360.0 * width as f64) as usize;
            let y = ((90.0 - lat) / 180.0 * height as f64) as usize;
            
            // Clamp to valid range
            let x = x.min(width - 1);
            let y = y.min(height - 1);
            
            // Check minimum distance to existing seeds
            let temp_seed = PlateSeed::new(0, 0, 0, lat, lon, 0.0, 0.0);
            let too_close = plate_seeds.iter().any(|seed: &PlateSeed| {
                seed.distance_to_km(&temp_seed) < MIN_SEED_DISTANCE
            });
            
            if !too_close {
                break (x, y, point_3d);
            }
        };

        let (lat, lon) = point_3d.to_lat_lon();

        plate_seeds.push(PlateSeed::new(
            (i + 1) as u16,
            x,
            y,
            lat,
            lon,
            rng.gen_range(0.0..360.0),
            rng.gen_range(2.0..10.0),
        ));
    }

    println!("Generated {} plate seeds with varied sizes", plate_seeds.len());
    Ok(plate_seeds)
}

/// Generate plates using spherical Voronoi diagram
pub fn generate_plates_voronoi(
    map: &mut TerrainMap<u16>,
    seeds: &[PlateSeed],
) {
    println!("Generating plates using spherical Voronoi method...");

    for y in 0..map.height {
        for x in 0..map.width {
            let point = map.get_spherical_point(x, y);

            let mut nearest_plate = 1;
            let mut min_distance = f64::INFINITY;

            // Find closest seed using geodesic distance
            for seed in seeds {
                let distance = point.distance_to(seed.spherical_point());
                if distance < min_distance {
                    min_distance = distance;
                    nearest_plate = seed.id;
                }
            }

            map.set(x, y, nearest_plate);
        }

        if y % 100 == 0 {
            println!("Progress: {:.1}%", (y as f64 / map.height as f64) * 100.0);
        }
    }
}

/// Generate plates using spherical region growing
pub fn generate_plates_region_growing(
    map: &mut TerrainMap<u16>,
    seeds: &[PlateSeed],
    plate_size_targets: &[f64],
    rng: &mut StdRng,
) -> Result<(), PlateError> {
    println!("Generating plates using spherical region growing...");

    // Initialize
    map.data.fill(0);
    let mut assigned_pixels = HashSet::new();
    let mut plate_pixel_counts = vec![0usize; seeds.len() + 1];
    let total_pixels = map.width * map.height;
    
    // Priority queue for growth frontier
    let mut frontier = BinaryHeap::new();
    
    // Place seeds and initialize frontier
    for seed in seeds {
        let idx = map.get_index(seed.x, seed.y);
        map.data[idx] = seed.id;
        assigned_pixels.insert(idx);
        plate_pixel_counts[seed.id as usize] = 1;
        
        // Add neighbors to frontier
        add_pixel_neighbors_to_frontier(
            map, seed.x, seed.y, seed.id, &mut frontier, &assigned_pixels, seeds
        );
    }
    
    // Grow regions
    let mut pixels_assigned = seeds.len();
    
    while let Some(growth_pixel) = frontier.pop() {
        let idx = map.get_index(growth_pixel.x, growth_pixel.y);
        
        // Skip if already assigned
        if assigned_pixels.contains(&idx) {
            continue;
        }
        
        // Check if this plate has reached its target size
        let plate_idx = growth_pixel.plate_id as usize;
        let current_fraction = plate_pixel_counts[plate_idx] as f64 / total_pixels as f64;
        let target_fraction = plate_size_targets[plate_idx - 1];
        
        // Add randomness but bias based on target size
        let growth_probability = if current_fraction < target_fraction {
            0.95 // High probability if under target
        } else {
            // Decreasing probability as we exceed target
            0.95 * (target_fraction / current_fraction).min(1.0)
        };
        
        if rng.gen::<f64>() < growth_probability {
            // Assign pixel
            map.data[idx] = growth_pixel.plate_id;
            assigned_pixels.insert(idx);
            plate_pixel_counts[plate_idx] += 1;
            pixels_assigned += 1;
            
            // Add neighbors to frontier
            add_pixel_neighbors_to_frontier(
                map, growth_pixel.x, growth_pixel.y, growth_pixel.plate_id,
                &mut frontier, &assigned_pixels, seeds
            );
        }
        
        // Progress update
        if pixels_assigned % 10000 == 0 {
            println!("Progress: {:.1}%", (pixels_assigned as f64 / total_pixels as f64) * 100.0);
        }
    }
    
    // Fill any remaining unassigned pixels with nearest plate
    fill_unassigned_pixels(map, seeds);
    
    Ok(())
}

/// Add pixel neighbors to growth frontier
fn add_pixel_neighbors_to_frontier(
    map: &TerrainMap<u16>,
    x: usize, 
    y: usize, 
    plate_id: u16, 
    frontier: &mut BinaryHeap<GrowthPixel>,
    assigned: &HashSet<usize>,
    seeds: &[PlateSeed],
) {
    let center_point = map.get_spherical_point(x, y);
    let neighbors = map.get_neighbors(x, y);
    
    // Find the seed point for this plate to calculate distances
    let seed_point = seeds.iter()
        .find(|s| s.id == plate_id)
        .map(|s| s.spherical_point())
        .unwrap_or(&center_point);
    
    for (nx, ny) in neighbors {
        let idx = map.get_index(nx, ny);
        
        // Skip if already assigned
        if assigned.contains(&idx) {
            continue;
        }
        
        let neighbor_point = map.get_spherical_point(nx, ny);
        let distance = seed_point.distance_to(&neighbor_point);
        
        frontier.push(GrowthPixel {
            x: nx,
            y: ny,
            plate_id,
            distance,
        });
    }
}

/// Fill any remaining unassigned pixels
fn fill_unassigned_pixels(map: &mut TerrainMap<u16>, seeds: &[PlateSeed]) {
    for y in 0..map.height {
        for x in 0..map.width {
            let idx = map.get_index(x, y);
            
            if map.data[idx] == 0 {
                let point = map.get_spherical_point(x, y);
                
                // Find nearest plate
                let mut nearest_plate = 1;
                let mut min_distance = f64::INFINITY;
                
                for seed in seeds {
                    let distance = point.distance_to(seed.spherical_point());
                    if distance < min_distance {
                        min_distance = distance;
                        nearest_plate = seed.id;
                    }
                }
                
                map.data[idx] = nearest_plate;
            }
        }
    }
}

/// Smooth boundaries with geodesic-aware smoothing
pub fn smooth_boundaries(map: &mut TerrainMap<u16>, iterations: usize) {
    println!("Smoothing boundaries ({} iterations)...", iterations);

    for iter in 0..iterations {
        let original_data = map.data.clone();
        
        for y in 0..map.height {
            for x in 0..map.width {
                let idx = map.get_index(x, y);
                let point = map.get_spherical_point(x, y);
                
                // Get neighbors and their weights based on geodesic distance
                let neighbors = map.get_neighbors(x, y);
                let mut weighted_counts: HashMap<u16, f64> = HashMap::new();
                let mut total_weight = 0.0;
                
                for (nx, ny) in neighbors {
                    let neighbor_idx = map.get_index(nx, ny);
                    let neighbor_point = map.get_spherical_point(nx, ny);
                    let neighbor_plate = original_data[neighbor_idx];
                    
                    // Weight by inverse geodesic distance
                    let distance = point.distance_to(&neighbor_point);
                    let weight = 1.0 / (distance + 0.01); // Add small epsilon to avoid division by zero
                    
                    *weighted_counts.entry(neighbor_plate).or_insert(0.0) += weight;
                    total_weight += weight;
                }
                
                // Also consider current pixel with some weight
                let current_plate = original_data[idx];
                let self_weight = total_weight * 0.5; // Current pixel has weight equal to half of all neighbors
                *weighted_counts.entry(current_plate).or_insert(0.0) += self_weight;
                total_weight += self_weight;
                
                // Find plate with highest weighted count
                if let Some((&most_common, &weight)) = weighted_counts
                    .iter()
                    .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
                {
                    // Only change if there's significant majority
                    if weight / total_weight > 0.6 {
                        map.data[idx] = most_common;
                    }
                }
            }
        }

        println!("Completed smoothing iteration {}/{}", iter + 1, iterations);
    }
}

/// Generate varied plate sizes using power law distribution
pub fn generate_plate_sizes(num_plates: usize, rng: &mut StdRng) -> Vec<f64> {
    let mut sizes: Vec<f64> = Vec::new();
    let alpha = 1.5; // Power law exponent
    
    for i in 0..num_plates {
        // Power law: larger index = smaller plate
        let size = 1.0 / ((i + 1) as f64).powf(alpha);
        sizes.push(size);
    }
    
    // Shuffle to randomize which plates get which sizes
    sizes.shuffle(rng);
    
    // Normalize so they sum to 1.0
    let sum: f64 = sizes.iter().sum();
    for size in &mut sizes {
        *size /= sum;
    }
    
    sizes
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::map::TerrainMap;

    #[test]
    fn test_plate_size_generation() {
        let mut rng = StdRng::seed_from_u64(42);
        let sizes = generate_plate_sizes(5, &mut rng);
        
        assert_eq!(sizes.len(), 5);
        
        // Should sum to approximately 1.0
        let sum: f64 = sizes.iter().sum();
        assert!((sum - 1.0).abs() < 0.0001);
        
        // All sizes should be positive
        assert!(sizes.iter().all(|&s| s > 0.0));
    }
    
    #[test]
    fn test_seed_generation() {
        let mut rng = StdRng::seed_from_u64(42);
        let sizes = generate_plate_sizes(3, &mut rng);
        let seeds = generate_seeds(100, 50, 3, &mut rng, &sizes);
        
        assert!(seeds.is_ok());
        let seeds = seeds.unwrap();
        assert_eq!(seeds.len(), 3);
        
        // Check that seeds have unique IDs
        for (i, seed) in seeds.iter().enumerate() {
            assert_eq!(seed.id, (i + 1) as u16);
        }
    }
}