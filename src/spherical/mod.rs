//! Spherical geometry and coordinate system utilities
//! 
//! This module provides mathematical utilities for working with spherical coordinates,
//! great circle distances, and 3D points on the unit sphere.

pub mod geometry;

pub use geometry::SphericalPoint;

/// Earth's radius in kilometers
pub const EARTH_RADIUS_KM: f64 = 6371.0;

/// Earth's surface area in square kilometers  
pub const EARTH_SURFACE_AREA_KM2: f64 = 4.0 * std::f64::consts::PI * EARTH_RADIUS_KM * EARTH_RADIUS_KM;

/// Convert great circle distance in radians to kilometers
pub fn radians_to_km(radians: f64) -> f64 {
    radians * EARTH_RADIUS_KM
}

/// Convert distance in kilometers to radians on Earth's surface
pub fn km_to_radians(km: f64) -> f64 {
    km / EARTH_RADIUS_KM
}

/// Calculate the bearing from one point to another (in degrees)
pub fn bearing(from: &SphericalPoint, to: &SphericalPoint) -> f64 {
    let (lat1, lon1) = from.to_lat_lon();
    let (lat2, lon2) = to.to_lat_lon();
    
    let lat1_rad = lat1.to_radians();
    let lat2_rad = lat2.to_radians();
    let dlon_rad = (lon2 - lon1).to_radians();
    
    let y = dlon_rad.sin() * lat2_rad.cos();
    let x = lat1_rad.cos() * lat2_rad.sin() - 
            lat1_rad.sin() * lat2_rad.cos() * dlon_rad.cos();
    
    let bearing_rad = y.atan2(x);
    (bearing_rad.to_degrees() + 360.0) % 360.0
}

/// Calculate the cross product of two 3D vectors
pub fn cross_product(a: &SphericalPoint, b: &SphericalPoint) -> SphericalPoint {
    SphericalPoint {
        x: a.y * b.z - a.z * b.y,
        y: a.z * b.x - a.x * b.z,  
        z: a.x * b.y - a.y * b.x,
    }
}

/// Calculate the dot product of two 3D vectors
pub fn dot_product(a: &SphericalPoint, b: &SphericalPoint) -> f64 {
    a.x * b.x + a.y * b.y + a.z * b.z
}

/// Get a point at a given distance and bearing from a starting point
pub fn point_at_distance_and_bearing(
    start: &SphericalPoint,
    distance_radians: f64,
    bearing_degrees: f64,
) -> SphericalPoint {
    let (lat1, lon1) = start.to_lat_lon();
    let lat1_rad = lat1.to_radians();
    let lon1_rad = lon1.to_radians();
    let bearing_rad = bearing_degrees.to_radians();
    
    let lat2_rad = (lat1_rad.sin() * distance_radians.cos() + 
                    lat1_rad.cos() * distance_radians.sin() * bearing_rad.cos()).asin();
    
    let lon2_rad = lon1_rad + (bearing_rad.sin() * distance_radians.sin() * lat1_rad.cos())
        .atan2(distance_radians.cos() - lat1_rad.sin() * lat2_rad.sin());
    
    SphericalPoint::from_lat_lon(lat2_rad.to_degrees(), lon2_rad.to_degrees())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_distance_conversions() {
        let distance_radians = std::f64::consts::PI / 4.0; // 45 degrees
        let distance_km = radians_to_km(distance_radians);
        let back_to_radians = km_to_radians(distance_km);
        
        assert!((back_to_radians - distance_radians).abs() < 0.0001);
    }
    
    #[test]
    fn test_bearing_calculation() {
        let north_pole = SphericalPoint::from_lat_lon(90.0, 0.0);
        let equator = SphericalPoint::from_lat_lon(0.0, 0.0);
        
        let bearing_to_north = bearing(&equator, &north_pole);
        assert!((bearing_to_north - 0.0).abs() < 0.1); // Should be approximately north (0 degrees)
    }
    
    #[test]
    fn test_point_at_distance_and_bearing() {
        let start = SphericalPoint::from_lat_lon(0.0, 0.0); // Equator, prime meridian
        let distance = std::f64::consts::PI / 2.0; // 90 degrees
        let bearing = 0.0; // North
        
        let end_point = point_at_distance_and_bearing(&start, distance, bearing);
        let (lat, _lon) = end_point.to_lat_lon();
        
        assert!((lat - 90.0).abs() < 0.1); // Should reach north pole
    }
}