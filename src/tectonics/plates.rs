//! Tectonic plate data structures and utilities

use crate::map::spherical::SphericalPoint;

/// Represents a tectonic plate seed point with motion information
#[derive(Debug, Clone)]
pub struct PlateSeed {
    pub id: u16,
    pub x: usize,
    pub y: usize,
    pub lon: f64,
    pub lat: f64,
    pub motion_direction: f64, // degrees
    pub motion_speed: f64,     // cm/year
    pub(crate) point_3d: SphericalPoint,
}

impl PlateSeed {
    /// Create a new plate seed
    pub fn new(
        id: u16, 
        x: usize, 
        y: usize, 
        lat: f64, 
        lon: f64, 
        motion_direction: f64, 
        motion_speed: f64
    ) -> Self {
        Self {
            id,
            x,
            y,
            lon,
            lat,
            motion_direction,
            motion_speed,
            point_3d: SphericalPoint::from_lat_lon(lat, lon),
        }
    }
    
    /// Get the 3D spherical point for this seed
    pub fn spherical_point(&self) -> &SphericalPoint {
        &self.point_3d
    }
    
    /// Calculate geodesic distance to another seed (in radians)
    pub fn distance_to(&self, other: &PlateSeed) -> f64 {
        self.point_3d.distance_to(&other.point_3d)
    }
    
    /// Calculate geodesic distance to another seed (in kilometers)
    pub fn distance_to_km(&self, other: &PlateSeed) -> f64 {
        const EARTH_RADIUS_KM: f64 = 6371.0;
        self.distance_to(other) * EARTH_RADIUS_KM
    }
}

/// Statistics about a generated tectonic plate
#[derive(Debug, Clone)]
pub struct PlateStats {
    pub pixels: usize,
    pub percentage: f64,
    pub area_km2: u64,
    pub seed: PlateSeed,
    pub plate_type: PlateType,
}

impl PlateStats {
    /// Create new plate statistics
    pub fn new(pixels: usize, area_km2: f64, seed: PlateSeed) -> Self {
        const PLANET_SURFACE_AREA_KM2: f64 = 4.0 * std::f64::consts::PI * 6371.0 * 6371.0;

        Self {
            pixels,
            percentage: (area_km2 / PLANET_SURFACE_AREA_KM2) * 100.0,
            area_km2: area_km2 as u64,
            seed,
            plate_type: PlateType::Oceanic, // Default, will be assigned later based on size
        }
    }

    /// Calculate size percentile of this plate among all plates
    ///
    /// Returns a value from 0.0 to 1.0 indicating what fraction of plates are smaller
    pub fn calculate_size_percentile(&self, all_stats: &std::collections::HashMap<u16, PlateStats>) -> f64 {
        let total_plates = all_stats.len();
        if total_plates == 0 {
            return 0.5;
        }

        let smaller_plates = all_stats.values()
            .filter(|s| s.area_km2 < self.area_km2)
            .count();

        smaller_plates as f64 / total_plates as f64
    }

    /// Assign plate type based on size among all plates
    pub fn assign_plate_type(&mut self, all_stats: &std::collections::HashMap<u16, PlateStats>) {
        let percentile = self.calculate_size_percentile(all_stats);
        self.plate_type = PlateType::from_size_percentile(percentile);
    }
}

/// Character/composition of a tectonic plate
///
/// Note: Plates are classified as primarily oceanic or continental based on
/// their dominant crust type. In reality, many plates contain both types
/// (e.g., North American plate has continental interior + Atlantic ocean floor),
/// but we classify based on the dominant character for geological modeling.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PlateType {
    /// Oceanic plate (denser, subducts beneath continental)
    Oceanic,
    /// Continental plate (lighter, does not subduct)
    Continental,
}

impl PlateType {
    /// Assign plate type based on size heuristics
    ///
    /// Larger plates tend to be oceanic (like Earth's Pacific plate), smaller plates
    /// tend to be continental. Uses a 70/30 split to roughly match Earth's ocean/land
    /// surface ratio (~71% ocean).
    ///
    /// # Arguments
    /// * `percentile` - Size percentile of the plate (0.0 = smallest, 1.0 = largest)
    ///
    /// # Returns
    /// PlateType classification based on size percentile
    pub fn from_size_percentile(percentile: f64) -> Self {
        if percentile > 0.3 {
            PlateType::Oceanic      // Top 70% by size = oceanic
        } else {
            PlateType::Continental  // Bottom 30% by size = continental
        }
    }
}

/// Type of plate interaction at boundaries
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PlateInteraction {
    /// Plates moving apart (divergent boundary)
    Divergent,
    /// Plates moving together (convergent boundary)
    Convergent,
    /// Plates sliding past each other (transform boundary)
    Transform,
}

/// Information about plate boundary characteristics
#[derive(Debug, Clone)]
pub struct PlateBoundary {
    pub plate_a: u16,
    pub plate_b: u16,
    pub interaction_type: PlateInteraction,
    pub relative_velocity: f64, // cm/year
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_plate_seed_creation() {
        let seed = PlateSeed::new(1, 100, 50, 45.0, -120.0, 90.0, 5.0);
        assert_eq!(seed.id, 1);
        assert_eq!(seed.lat, 45.0);
        assert_eq!(seed.lon, -120.0);
        assert_eq!(seed.motion_direction, 90.0);
        assert_eq!(seed.motion_speed, 5.0);
    }
    
    #[test]
    fn test_distance_calculation() {
        let seed1 = PlateSeed::new(1, 0, 0, 0.0, 0.0, 0.0, 0.0);
        let seed2 = PlateSeed::new(2, 0, 0, 0.0, 90.0, 0.0, 0.0);
        
        // 90 degrees longitude difference at equator should be PI/2 radians
        let distance = seed1.distance_to(&seed2);
        assert!((distance - std::f64::consts::PI / 2.0).abs() < 0.0001);
        
        // Same point should have zero distance
        let distance_same = seed1.distance_to(&seed1);
        assert!(distance_same < 0.0001);
    }
}