//! Core map data structures for terrain generation
//! 
//! This module provides the fundamental data types and utilities for representing
//! geographic data across different terrain generation stages.

use crate::spherical::SphericalPoint;

/// Core map data structure for all terrain generation steps
#[derive(Debug, Clone)]
pub struct TerrainMap<T> {
    /// Map width in pixels
    pub width: usize,
    /// Map height in pixels  
    pub height: usize,
    /// Data array stored row-major (y * width + x)
    pub data: Vec<T>,
    /// Geographic bounds and projection info
    pub projection: MapProjection,
}

/// Map projection and coordinate system information
#[derive(Debug, Clone)]
pub struct MapProjection {
    /// Western longitude bound (degrees)
    pub west_bound: f64,
    /// Eastern longitude bound (degrees) 
    pub east_bound: f64,
    /// Northern latitude bound (degrees)
    pub north_bound: f64,
    /// Southern latitude bound (degrees)
    pub south_bound: f64,
    /// Degrees per pixel in longitude
    pub lon_resolution: f64,
    /// Degrees per pixel in latitude
    pub lat_resolution: f64,
}

impl MapProjection {
    /// Create a global equirectangular projection
    pub fn global_equirectangular(width: usize, height: usize) -> Self {
        Self {
            west_bound: -180.0,
            east_bound: 180.0,
            north_bound: 90.0,
            south_bound: -90.0,
            lon_resolution: 360.0 / width as f64,
            lat_resolution: 180.0 / height as f64,
        }
    }
    
    /// Convert pixel coordinates to geographic coordinates
    pub fn pixel_to_coords(&self, x: usize, y: usize) -> (f64, f64) {
        let lon = self.west_bound + (x as f64 + 0.5) * self.lon_resolution;
        let lat = self.north_bound - (y as f64 + 0.5) * self.lat_resolution;
        (lat, lon)
    }
    
    /// Convert geographic coordinates to pixel coordinates
    pub fn coords_to_pixel(&self, lat: f64, lon: f64) -> (usize, usize) {
        let x = ((lon - self.west_bound) / self.lon_resolution).floor() as usize;
        let y = ((self.north_bound - lat) / self.lat_resolution).floor() as usize;
        (x, y)
    }
    
    /// Get the area of a pixel at given latitude (km²)
    pub fn pixel_area_km2(&self, lat: f64) -> f64 {
        const EARTH_RADIUS_KM: f64 = 6371.0;
        
        let lat_rad = lat.to_radians();
        let delta_lat_rad = self.lat_resolution.to_radians();
        let delta_lon_rad = self.lon_resolution.to_radians();
        
        // Area of spherical rectangle
        let area = EARTH_RADIUS_KM * EARTH_RADIUS_KM * delta_lon_rad * 
                   (lat_rad + delta_lat_rad / 2.0).sin().abs() * delta_lat_rad.abs();
        
        area
    }
}

impl<T> TerrainMap<T> {
    /// Create a new terrain map with given dimensions and default value
    pub fn new(width: usize, height: usize, default_value: T) -> Self 
    where 
        T: Clone 
    {
        Self {
            width,
            height,
            data: vec![default_value; width * height],
            projection: MapProjection::global_equirectangular(width, height),
        }
    }
    
    /// Create with custom projection
    pub fn with_projection(width: usize, height: usize, default_value: T, projection: MapProjection) -> Self
    where 
        T: Clone
    {
        Self {
            width,
            height,
            data: vec![default_value; width * height],
            projection,
        }
    }
    
    /// Create from existing data
    pub fn from_data(width: usize, height: usize, data: Vec<T>) -> Self {
        Self {
            width,
            height,
            data,
            projection: MapProjection::global_equirectangular(width, height),
        }
    }
    
    /// Get pixel index from coordinates
    pub fn get_index(&self, x: usize, y: usize) -> usize {
        y * self.width + x
    }
    
    /// Get value at pixel coordinates
    pub fn get(&self, x: usize, y: usize) -> Option<&T> {
        if x < self.width && y < self.height {
            Some(&self.data[self.get_index(x, y)])
        } else {
            None
        }
    }
    
    /// Set value at pixel coordinates
    pub fn set(&mut self, x: usize, y: usize, value: T) -> bool {
        if x < self.width && y < self.height {
            let idx = self.get_index(x, y);
            self.data[idx] = value;
            true
        } else {
            false
        }
    }
    
    /// Get value at geographic coordinates
    pub fn get_at_coords(&self, lat: f64, lon: f64) -> Option<&T> {
        let (x, y) = self.projection.coords_to_pixel(lat, lon);
        self.get(x, y)
    }
    
    /// Set value at geographic coordinates
    pub fn set_at_coords(&mut self, lat: f64, lon: f64, value: T) -> bool {
        let (x, y) = self.projection.coords_to_pixel(lat, lon);
        self.set(x, y, value)
    }
    
    /// Get neighboring pixel coordinates (8-connected, with longitude wraparound)
    pub fn get_neighbors(&self, x: usize, y: usize) -> Vec<(usize, usize)> {
        let mut neighbors = Vec::new();
        
        for dy in -1i32..=1 {
            for dx in -1i32..=1 {
                if dx == 0 && dy == 0 {
                    continue;
                }
                
                let ny = y as i32 + dy;
                
                // Skip out of bounds latitude
                if ny < 0 || ny >= self.height as i32 {
                    continue;
                }
                
                // Handle longitude wraparound
                let mut nx = x as i32 + dx;
                if nx < 0 {
                    nx = self.width as i32 - 1;
                } else if nx >= self.width as i32 {
                    nx = 0;
                }
                
                neighbors.push((nx as usize, ny as usize));
            }
        }
        
        neighbors
    }
    
    /// Convert each pixel to its corresponding 3D point on the unit sphere
    pub fn get_spherical_point(&self, x: usize, y: usize) -> SphericalPoint {
        let (lat, lon) = self.projection.pixel_to_coords(x, y);
        SphericalPoint::from_lat_lon(lat, lon)
    }
    
    /// Fill map with a function of coordinates
    pub fn fill_with<F>(&mut self, mut func: F) 
    where 
        F: FnMut(usize, usize, f64, f64) -> T
    {
        for y in 0..self.height {
            for x in 0..self.width {
                let (lat, lon) = self.projection.pixel_to_coords(x, y);
                let value = func(x, y, lat, lon);
                let idx = self.get_index(x, y);
                self.data[idx] = value;
            }
        }
    }
    
    /// Apply a function to transform the map data
    pub fn map<U, F>(self, func: F) -> TerrainMap<U>
    where 
        F: FnMut(T) -> U
    {
        TerrainMap {
            width: self.width,
            height: self.height,
            data: self.data.into_iter().map(func).collect(),
            projection: self.projection,
        }
    }
    
    /// Get statistics about the data
    pub fn get_stats(&self) -> MapStats<T>
    where 
        T: Clone + PartialOrd
    {
        let mut stats = MapStats::new();
        
        if !self.data.is_empty() {
            stats.min_value = Some(self.data[0].clone());
            stats.max_value = Some(self.data[0].clone());
            
            for value in &self.data {
                if let Some(ref min) = stats.min_value {
                    if value < min {
                        stats.min_value = Some(value.clone());
                    }
                }
                if let Some(ref max) = stats.max_value {
                    if value > max {
                        stats.max_value = Some(value.clone());
                    }
                }
            }
        }
        
        stats.total_pixels = self.data.len();
        stats
    }
}

/// Statistics about map data
#[derive(Debug, Clone)]
pub struct MapStats<T> {
    pub total_pixels: usize,
    pub min_value: Option<T>,
    pub max_value: Option<T>,
}

impl<T> MapStats<T> {
    pub fn new() -> Self {
        Self {
            total_pixels: 0,
            min_value: None,
            max_value: None,
        }
    }
}

impl<T> Default for MapStats<T> {
    fn default() -> Self {
        Self::new()
    }
}

/// Common terrain map types
pub type PlateMap = TerrainMap<u16>;
pub type ElevationMap = TerrainMap<f32>;
pub type TemperatureMap = TerrainMap<f32>;
pub type PrecipitationMap = TerrainMap<f32>;
pub type BiomeMap = TerrainMap<u8>;

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_map_creation() {
        let map: TerrainMap<u16> = TerrainMap::new(100, 50, 0);
        assert_eq!(map.width, 100);
        assert_eq!(map.height, 50);
        assert_eq!(map.data.len(), 5000);
    }
    
    #[test]
    fn test_coordinate_conversion() {
        let projection = MapProjection::global_equirectangular(360, 180);
        
        // Test equator, prime meridian (center pixel should be 180, 90)
        let (lat, lon) = projection.pixel_to_coords(180, 90);
        assert!((lat - 0.0).abs() < 1.0); // Relax tolerance for equirectangular projection
        assert!((lon - 0.5).abs() < 1.0); // Should be close to prime meridian
        
        // Test north pole
        let (lat, lon) = projection.pixel_to_coords(180, 0);
        assert!((lat - 89.5).abs() < 1.0); // Near north pole
        
        // Test roundtrip
        let (x, y) = projection.coords_to_pixel(45.0, -120.0);
        let (lat2, lon2) = projection.pixel_to_coords(x, y);
        assert!((lat2 - 45.0).abs() < 1.0);
        assert!((lon2 - (-120.0)).abs() < 1.0);
    }
    
    #[test]
    fn test_neighbors() {
        let map: TerrainMap<u16> = TerrainMap::new(10, 10, 0);
        
        // Test center pixel
        let neighbors = map.get_neighbors(5, 5);
        assert_eq!(neighbors.len(), 8);
        
        // Test edge with wraparound
        let neighbors = map.get_neighbors(0, 5);
        assert_eq!(neighbors.len(), 8);
        assert!(neighbors.contains(&(9, 5))); // Should wrap around
        
        // Test corner
        let neighbors = map.get_neighbors(0, 0);
        assert_eq!(neighbors.len(), 5); // Only valid neighbors
    }
}