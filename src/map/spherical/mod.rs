//! Spherical geometry and coordinate system utilities
//! 
//! This module provides mathematical utilities for working with spherical coordinates,
//! great circle distances, and 3D points on the unit sphere.

pub mod geometry;
mod planetary_params_tests;

pub use geometry::SphericalPoint;

/// Earth's radius in kilometers
pub const EARTH_RADIUS_KM: f64 = 6371.0;

/// Earth's surface area in square kilometers  
pub const EARTH_SURFACE_AREA_KM2: f64 = 4.0 * std::f64::consts::PI * EARTH_RADIUS_KM * EARTH_RADIUS_KM;

/// Planetary parameters for world generation
#[derive(Debug, Clone, PartialEq)]
pub struct PlanetaryParams {
    // Physical properties
    /// Planet radius in kilometers
    pub radius_km: f64,
    /// Planet surface area in square kilometers
    pub surface_area_km2: f64,
    /// Planet gravity in m/s² (affects geological processes)
    pub gravity_ms2: f64,
    /// Planet mass in kg (derived parameter)
    pub mass_kg: f64,
    /// Planet density in kg/m³ (affects internal structure)
    pub density_kgm3: f64,
    
    // Rotational properties
    /// Axial tilt in degrees (affects seasonal variation)
    pub axial_tilt_degrees: f64,
    /// Rotation period in hours (day length)
    pub rotation_period_hours: f64,
    
    // Orbital properties
    /// Orbital period in Earth days (year length)
    pub orbital_period_days: f64,
    /// Orbital eccentricity (0 = perfect circle, 0.99 = highly elliptical)
    pub orbital_eccentricity: f64,
    /// Semi-major axis in AU (average orbital distance)
    pub semi_major_axis_au: f64,
    /// Perihelion distance in AU (closest approach to star)
    pub perihelion_au: f64,
    /// Aphelion distance in AU (farthest distance from star)
    pub aphelion_au: f64,
    /// Orbital inclination in degrees (relative to ecliptic)
    pub orbital_inclination_degrees: f64,
    
    // Atmospheric properties
    /// Atmospheric pressure at sea level in kPa
    pub atmospheric_pressure_kpa: f64,
    /// Normalized greenhouse effect (1.0 = Earth-like, 0.18 = Mars, ~15 = Venus)
    pub greenhouse_factor: f64,
    
    // Stellar properties
    /// Stellar luminosity relative to Sun (1.0 = solar luminosity)
    pub stellar_luminosity: f64,
}

impl PlanetaryParams {
    /// Create Earth-like planetary parameters
    pub fn earth() -> Self {
        Self {
            // Physical properties
            radius_km: EARTH_RADIUS_KM,
            surface_area_km2: EARTH_SURFACE_AREA_KM2,
            gravity_ms2: 9.81,
            mass_kg: 5.972e24,
            density_kgm3: 5515.0,
            
            // Rotational properties
            axial_tilt_degrees: 23.44,
            rotation_period_hours: 24.0,
            
            // Orbital properties
            orbital_period_days: 365.2,
            orbital_eccentricity: 0.017,
            semi_major_axis_au: 1.000,
            perihelion_au: 0.983, // 147.1 million km / 149.6 million km
            aphelion_au: 1.017,   // 152.1 million km / 149.6 million km
            orbital_inclination_degrees: 0.0,
            
            // Atmospheric properties
            atmospheric_pressure_kpa: 101.325,
            greenhouse_factor: 1.0, // Earth baseline
            
            // Stellar properties
            stellar_luminosity: 1.0, // Solar luminosity
        }
    }
    
    /// Create custom planetary parameters from radius and density
    pub fn from_radius_and_density(radius_km: f64, density_kgm3: f64) -> Self {
        let radius_m = radius_km * 1000.0;
        let volume_m3 = 4.0 / 3.0 * std::f64::consts::PI * radius_m.powi(3);
        let mass_kg = volume_m3 * density_kgm3;
        let surface_area_km2 = 4.0 * std::f64::consts::PI * radius_km * radius_km;
        
        // Calculate surface gravity: g = GM/r²
        const G: f64 = 6.67430e-11; // Gravitational constant
        let gravity_ms2 = G * mass_kg / (radius_m * radius_m);
        
        Self {
            radius_km,
            surface_area_km2,
            gravity_ms2,
            mass_kg,
            density_kgm3,
            axial_tilt_degrees: 23.44, // Default to Earth-like
            rotation_period_hours: 24.0,
            orbital_period_days: 365.2,
            orbital_eccentricity: 0.017,
            semi_major_axis_au: 1.000,
            perihelion_au: 0.983,
            aphelion_au: 1.017,
            orbital_inclination_degrees: 0.0,
            atmospheric_pressure_kpa: 101.325,
            greenhouse_factor: 1.0, // Default to Earth-like
            stellar_luminosity: 1.0, // Default to solar luminosity
        }
    }
    
    /// Create planetary parameters from radius (assuming Earth-like density and other properties)
    pub fn from_radius(radius_km: f64) -> Self {
        Self::from_radius_and_density(radius_km, 5515.0)
    }
    
    /// Create a Mars-like planet
    pub fn mars() -> Self {
        Self {
            radius_km: 3389.5,
            surface_area_km2: 4.0 * std::f64::consts::PI * 3389.5 * 3389.5,
            gravity_ms2: 3.71,
            mass_kg: 6.39e23,
            density_kgm3: 3933.0,
            axial_tilt_degrees: 25.19,
            rotation_period_hours: 24.62,
            orbital_period_days: 687.0,
            orbital_eccentricity: 0.094,
            semi_major_axis_au: 1.524, // 228.0 million km / 149.6 million km
            perihelion_au: 1.381,      // 206.7 million km / 149.6 million km
            aphelion_au: 1.666,        // 249.3 million km / 149.6 million km
            orbital_inclination_degrees: 1.85,
            atmospheric_pressure_kpa: 0.636,
            greenhouse_factor: 0.18, // 6°C / 33°C = 0.18 relative to Earth
            stellar_luminosity: 1.0, // Same Sun as Earth
        }
    }
    
    /// Create a Venus-like planet
    pub fn venus() -> Self {
        Self {
            radius_km: 6051.8,
            surface_area_km2: 4.0 * std::f64::consts::PI * 6051.8 * 6051.8,
            gravity_ms2: 8.87,
            mass_kg: 4.87e24,
            density_kgm3: 5243.0,
            axial_tilt_degrees: 177.4, // Retrograde rotation
            rotation_period_hours: 5832.5, // 243 Earth days
            orbital_period_days: 224.7,
            orbital_eccentricity: 0.007,
            semi_major_axis_au: 0.723, // 108.2 million km / 149.6 million km
            perihelion_au: 0.718,      // 107.5 million km / 149.6 million km  
            aphelion_au: 0.728,        // 108.9 million km / 149.6 million km
            orbital_inclination_degrees: 3.39,
            atmospheric_pressure_kpa: 9200.0, // 92 bar
            greenhouse_factor: 15.15, // 500°C / 33°C = 15.15 relative to Earth
            stellar_luminosity: 1.0, // Same Sun as Earth
        }
    }
    
    /// Convert great circle distance in radians to kilometers for this planet
    pub fn radians_to_km(&self, radians: f64) -> f64 {
        radians * self.radius_km
    }
    
    /// Convert distance in kilometers to radians for this planet
    pub fn km_to_radians(&self, km: f64) -> f64 {
        km / self.radius_km
    }
    
    /// Calculate escape velocity in km/s
    pub fn escape_velocity_kms(&self) -> f64 {
        const G: f64 = 6.67430e-11;
        let radius_m = self.radius_km * 1000.0;
        ((2.0 * G * self.mass_kg) / radius_m).sqrt() / 1000.0
    }
    
    /// Get the effective temperature range due to axial tilt (seasonal variation factor)
    pub fn seasonal_variation_factor(&self) -> f64 {
        // Higher axial tilt = more seasonal variation
        (self.axial_tilt_degrees / 90.0).min(1.0)
    }
    
    /// Get day/night temperature variation factor based on rotation period
    pub fn diurnal_variation_factor(&self) -> f64 {
        // Longer days = more extreme day/night temperature differences
        // Venus with its 243-day rotation has extreme temperature uniformity
        // Fast-rotating planets have less temperature variation
        let earth_ratio = self.rotation_period_hours / 24.0;
        (earth_ratio.ln() + 1.0).clamp(0.1, 5.0)
    }
    
    /// Calculate solar radiation variation due to orbital eccentricity
    pub fn orbital_radiation_variation(&self) -> f64 {
        // Solar flux at perihelion vs aphelion
        // Flux is proportional to 1/distance²
        let perihelion_flux = 1.0 / (self.perihelion_au * self.perihelion_au);
        let aphelion_flux = 1.0 / (self.aphelion_au * self.aphelion_au);
        perihelion_flux / aphelion_flux - 1.0 // Return as fractional variation
    }
    
    /// Get the year length relative to Earth
    pub fn year_length_factor(&self) -> f64 {
        self.orbital_period_days / 365.2
    }
    
    /// Calculate average solar flux relative to Earth (1.0 = Earth distance)
    pub fn average_solar_flux(&self) -> f64 {
        // Flux is proportional to luminosity/distance²
        self.stellar_luminosity / (self.semi_major_axis_au * self.semi_major_axis_au)
    }
    
    /// Calculate actual insolation at perihelion in W/m²
    pub fn perihelion_insolation_wm2(&self) -> f64 {
        const SOLAR_CONSTANT: f64 = 1361.0; // W/m² at 1 AU from Sun
        SOLAR_CONSTANT * self.stellar_luminosity / (self.perihelion_au * self.perihelion_au)
    }
    
    /// Calculate actual insolation at aphelion in W/m²
    pub fn aphelion_insolation_wm2(&self) -> f64 {
        const SOLAR_CONSTANT: f64 = 1361.0; // W/m² at 1 AU from Sun
        SOLAR_CONSTANT * self.stellar_luminosity / (self.aphelion_au * self.aphelion_au)
    }
    
    /// Calculate average insolation in W/m²
    pub fn average_insolation_wm2(&self) -> f64 {
        const SOLAR_CONSTANT: f64 = 1361.0; // W/m² at 1 AU from Sun
        SOLAR_CONSTANT * self.stellar_luminosity / (self.semi_major_axis_au * self.semi_major_axis_au)
    }
    
    /// Calculate insolation at a specific orbital position (0.0 = perihelion, 1.0 = aphelion)
    pub fn insolation_at_position(&self, orbital_position: f64) -> f64 {
        // Linear interpolation between perihelion and aphelion for simplicity
        // In reality, this would require solving Kepler's equation
        let position_clamped = orbital_position.clamp(0.0, 1.0);
        let distance_au = self.perihelion_au + position_clamped * (self.aphelion_au - self.perihelion_au);
        
        const SOLAR_CONSTANT: f64 = 1361.0; // W/m² at 1 AU from Sun
        SOLAR_CONSTANT * self.stellar_luminosity / (distance_au * distance_au)
    }
}

impl Default for PlanetaryParams {
    fn default() -> Self {
        Self::earth()
    }
}

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