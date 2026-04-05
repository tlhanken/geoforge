use crate::objects::materials::{
    AtmosphereComposition,
    HydrosphereComposition,
    LithosphereComposition,
};

#[derive(Debug, Clone, PartialEq)]
pub struct PlanetaryBody {
    /// Planet radius in kilometers
    pub radius_km: f64,
    /// Planet surface area in square kilometers
    pub surface_area_km2: f64,
    /// Planet surface gravity in m/s² (affects geological processes)
    pub surface_gravity_ms2: f64,
    /// Planet mass in kg (derived parameter)
    pub mass_kg: f64,
    /// Planet density in kg/m³ (affects internal structure)
    pub density_kgm3: f64,

    /// Axial tilt in degrees (affects seasonal variation)
    pub axial_tilt_degrees: f64,
    /// Rotation period in hours (day length)
    pub rotation_period_hours: f64,

    /// Atmosphere
    pub atmosphere: Option<Atmosphere>,
    /// Hydrosphere
    pub hydrosphere: Option<Hydrosphere>,
    /// Lithosphere
    pub lithosphere: Option<Lithosphere>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Atmosphere {
    /// Atmospheric Composition
    pub composition: AtmosphereComposition,
    /// Atmospheric pressure at sea level in kPa
    pub atmospheric_pressure_kpa: f64,
    /// Normalized greenhouse effect (1.0 = Earth-like, 0.18 = Mars, ~15 = Venus)
    pub greenhouse_factor: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Hydrosphere {
    /// Ocean Liquid
    pub composition: HydrosphereComposition,
    /// Ocean Coverage
    pub ocean_coverage: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Lithosphere {
    /// Crust Composition
    pub composition: LithosphereComposition,
    /// Tectonic Activity
    pub tectonic_activity: Option<TectonicActivity>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum TectonicActivity {
    Low,
    Medium,
    High,
}