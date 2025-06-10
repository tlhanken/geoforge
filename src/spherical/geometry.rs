//! Spherical geometry utilities for tectonic plate generation

use std::f64::consts::PI;

/// 3D point on unit sphere
#[derive(Debug, Clone, Copy)]
pub struct SphericalPoint {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl SphericalPoint {
    /// Create from latitude/longitude (in degrees)
    pub fn from_lat_lon(lat_deg: f64, lon_deg: f64) -> Self {
        let lat = lat_deg.to_radians();
        let lon = lon_deg.to_radians();
        
        Self {
            x: lat.cos() * lon.cos(),
            y: lat.cos() * lon.sin(),
            z: lat.sin(),
        }
    }
    
    /// Convert back to latitude/longitude (in degrees)
    pub fn to_lat_lon(&self) -> (f64, f64) {
        let lat = self.z.asin().to_degrees();
        let lon = self.y.atan2(self.x).to_degrees();
        (lat, lon)
    }
    
    /// Great circle distance to another point (on unit sphere)
    pub fn distance_to(&self, other: &SphericalPoint) -> f64 {
        // Use dot product for angle between vectors
        let dot = self.x * other.x + self.y * other.y + self.z * other.z;
        // Clamp to avoid numerical errors
        let dot = dot.max(-1.0).min(1.0);
        dot.acos()
    }
    
    /// Normalize to unit sphere
    pub fn normalize(&mut self) {
        let len = (self.x * self.x + self.y * self.y + self.z * self.z).sqrt();
        if len > 0.0 {
            self.x /= len;
            self.y /= len;
            self.z /= len;
        }
    }
    
    /// Linear interpolation between two points on the sphere
    pub fn slerp(&self, other: &SphericalPoint, t: f64) -> SphericalPoint {
        let omega = self.distance_to(other);
        
        if omega < 0.0001 {
            // Points are very close, use linear interpolation
            return *self;
        }
        
        let sin_omega = omega.sin();
        let a = ((1.0 - t) * omega).sin() / sin_omega;
        let b = (t * omega).sin() / sin_omega;
        
        SphericalPoint {
            x: a * self.x + b * other.x,
            y: a * self.y + b * other.y,
            z: a * self.z + b * other.z,
        }
    }
    
    /// Create a random point on the unit sphere
    pub fn random(rng: &mut impl rand::Rng) -> Self {
        // Use rejection sampling for uniform distribution
        let z = rng.gen_range(-1.0..1.0);
        let theta = rng.gen_range(0.0..2.0 * PI);
        let r = (1.0 - z * z).sqrt();
        
        SphericalPoint {
            x: r * theta.cos(),
            y: r * theta.sin(),
            z,
        }
    }
    
    /// Get neighboring points at a given angular distance
    pub fn get_neighbors(&self, angular_distance: f64, num_samples: usize) -> Vec<SphericalPoint> {
        let mut neighbors = Vec::with_capacity(num_samples);
        
        // Convert to lat/lon for easier manipulation
        let (lat, lon) = self.to_lat_lon();
        let lat_rad = lat.to_radians();
        let lon_rad = lon.to_radians();
        
        for i in 0..num_samples {
            let angle = 2.0 * PI * i as f64 / num_samples as f64;
            
            // Calculate new position using small angle approximation
            let dlat = angular_distance * angle.sin();
            let dlon = angular_distance * angle.cos() / lat_rad.cos().max(0.001);
            
            let new_lat = (lat_rad + dlat).max(-PI/2.0).min(PI/2.0);
            let new_lon = lon_rad + dlon;
            
            let neighbor = SphericalPoint {
                x: new_lat.cos() * new_lon.cos(),
                y: new_lat.cos() * new_lon.sin(),
                z: new_lat.sin(),
            };
            
            neighbors.push(neighbor);
        }
        
        neighbors
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_spherical_point_creation() {
        // Test North Pole
        let north_pole = SphericalPoint::from_lat_lon(90.0, 0.0);
        assert!((north_pole.z - 1.0).abs() < 0.0001);
        assert!(north_pole.x.abs() < 0.0001);
        assert!(north_pole.y.abs() < 0.0001);
        
        // Test South Pole
        let south_pole = SphericalPoint::from_lat_lon(-90.0, 0.0);
        assert!((south_pole.z + 1.0).abs() < 0.0001);
        
        // Test Equator
        let equator = SphericalPoint::from_lat_lon(0.0, 0.0);
        assert!(equator.z.abs() < 0.0001);
        assert!((equator.x - 1.0).abs() < 0.0001);
    }
    
    #[test]
    fn test_distance_calculation() {
        let p1 = SphericalPoint::from_lat_lon(0.0, 0.0);
        let p2 = SphericalPoint::from_lat_lon(0.0, 90.0);
        
        // 90 degrees should be PI/2 radians
        let dist = p1.distance_to(&p2);
        assert!((dist - PI / 2.0).abs() < 0.0001);
        
        // Same point should have zero distance
        let dist_same = p1.distance_to(&p1);
        assert!(dist_same < 0.0001);
    }
    
    #[test]
    fn test_conversion_roundtrip() {
        let original_lat = 45.0;
        let original_lon = -120.0;
        
        let point = SphericalPoint::from_lat_lon(original_lat, original_lon);
        let (lat, lon) = point.to_lat_lon();
        
        assert!((lat - original_lat).abs() < 0.0001);
        assert!((lon - original_lon).abs() < 0.0001);
    }
}
