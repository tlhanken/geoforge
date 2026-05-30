//! Planetary physical, orbital, and stellar parameters.

use serde::{Deserialize, Serialize};

/// Mean Earth radius (km).
pub const EARTH_RADIUS_KM: f64 = 6371.0;

/// Earth surface area (km²).
pub const EARTH_SURFACE_AREA_KM2: f64 =
    4.0 * std::f64::consts::PI * EARTH_RADIUS_KM * EARTH_RADIUS_KM;

/// Physical and orbital parameters for a planetary body.
///
/// Ported from Geoforge v1 `PlanetaryParams` — f64 SI-ish units for simplicity and
/// serde stability. Physical-unit crates may wrap this later.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PlanetaryParams {
    // --- Physical ---
    /// Planet radius (km).
    pub radius_km: f64,
    /// Total surface area (km²).
    pub surface_area_km2: f64,
    /// Surface gravity (m/s²).
    pub gravity_ms2: f64,
    /// Total mass (kg).
    pub mass_kg: f64,
    /// Mean density (kg/m³).
    pub density_kgm3: f64,

    // --- Rotation ---
    /// Axial tilt (degrees).
    pub axial_tilt_degrees: f64,
    /// Sidereal rotation period (hours).
    pub rotation_period_hours: f64,

    // --- Orbit ---
    /// Orbital period (Earth days).
    pub orbital_period_days: f64,
    /// Orbital eccentricity.
    pub orbital_eccentricity: f64,
    /// Semi-major axis (AU).
    pub semi_major_axis_au: f64,
    /// Perihelion (AU).
    pub perihelion_au: f64,
    /// Aphelion (AU).
    pub aphelion_au: f64,
    /// Orbital inclination (degrees).
    pub orbital_inclination_degrees: f64,

    // --- Atmosphere ---
    /// Sea-level pressure (kPa).
    pub atmospheric_pressure_kpa: f64,
    /// Greenhouse factor relative to Earth (1.0 = Earth-like).
    pub greenhouse_factor: f64,

    // --- Stellar ---
    /// Stellar luminosity relative to the Sun.
    pub stellar_luminosity: f64,
}

impl PlanetaryParams {
    /// Earth-like parameters.
    #[must_use]
    pub fn earth() -> Self {
        Self {
            radius_km: EARTH_RADIUS_KM,
            surface_area_km2: EARTH_SURFACE_AREA_KM2,
            gravity_ms2: 9.81,
            mass_kg: 5.972e24,
            density_kgm3: 5515.0,
            axial_tilt_degrees: 23.44,
            rotation_period_hours: 24.0,
            orbital_period_days: 365.25,
            orbital_eccentricity: 0.017,
            semi_major_axis_au: 1.0,
            perihelion_au: 0.983,
            aphelion_au: 1.017,
            orbital_inclination_degrees: 0.0,
            atmospheric_pressure_kpa: 101.325,
            greenhouse_factor: 1.0,
            stellar_luminosity: 1.0,
        }
    }

    /// Mars-like parameters.
    #[must_use]
    pub fn mars() -> Self {
        Self {
            radius_km: 3389.5,
            surface_area_km2: 4.0 * std::f64::consts::PI * 3389.5_f64.powi(2),
            gravity_ms2: 3.71,
            mass_kg: 6.39e23,
            density_kgm3: 3933.0,
            axial_tilt_degrees: 25.19,
            rotation_period_hours: 24.62,
            orbital_period_days: 687.0,
            orbital_eccentricity: 0.094,
            semi_major_axis_au: 1.524,
            perihelion_au: 1.381,
            aphelion_au: 1.666,
            orbital_inclination_degrees: 1.85,
            atmospheric_pressure_kpa: 0.636,
            greenhouse_factor: 0.18,
            stellar_luminosity: 1.0,
        }
    }

    /// Venus-like parameters.
    #[must_use]
    pub fn venus() -> Self {
        Self {
            radius_km: 6051.8,
            surface_area_km2: 4.0 * std::f64::consts::PI * 6051.8_f64.powi(2),
            gravity_ms2: 8.87,
            mass_kg: 4.87e24,
            density_kgm3: 5243.0,
            axial_tilt_degrees: 177.4,
            rotation_period_hours: 5832.5,
            orbital_period_days: 224.7,
            orbital_eccentricity: 0.007,
            semi_major_axis_au: 0.723,
            perihelion_au: 0.718,
            aphelion_au: 0.728,
            orbital_inclination_degrees: 3.39,
            atmospheric_pressure_kpa: 9200.0,
            greenhouse_factor: 15.15,
            stellar_luminosity: 1.0,
        }
    }

    /// Build from radius and density; other fields default to Earth-like.
    #[must_use]
    pub fn from_radius_and_density(radius_km: f64, density_kgm3: f64) -> Self {
        const G: f64 = 6.674_30e-11;
        let radius_m = radius_km * 1000.0;
        let volume_m3 = (4.0 / 3.0) * std::f64::consts::PI * radius_m.powi(3);
        let mass_kg = volume_m3 * density_kgm3;
        let gravity_ms2 = G * mass_kg / radius_m.powi(2);
        let surface_area_km2 = 4.0 * std::f64::consts::PI * radius_km * radius_km;

        let mut params = Self::earth();
        params.radius_km = radius_km;
        params.surface_area_km2 = surface_area_km2;
        params.mass_kg = mass_kg;
        params.density_kgm3 = density_kgm3;
        params.gravity_ms2 = gravity_ms2;
        params
    }

    /// Top-of-atmosphere insolation (W/m²) from inverse-square law.
    #[must_use]
    pub fn insolation_wm2(&self) -> f64 {
        const SOLAR_CONSTANT: f64 = 1361.0;
        SOLAR_CONSTANT * self.stellar_luminosity / self.semi_major_axis_au.powi(2)
    }

    /// Convert great-circle radians to kilometers for this body.
    #[must_use]
    pub fn radians_to_km(&self, radians: f64) -> f64 {
        radians * self.radius_km
    }

    /// Convert kilometers to great-circle radians.
    #[must_use]
    pub fn km_to_radians(&self, km: f64) -> f64 {
        km / self.radius_km
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn earth_insolation_near_solar_constant() {
        let earth = PlanetaryParams::earth();
        assert!((earth.insolation_wm2() - 1361.0).abs() < 1.0);
    }

    #[test]
    fn mars_smaller_than_earth() {
        assert!(PlanetaryParams::mars().radius_km < PlanetaryParams::earth().radius_km);
    }
}
