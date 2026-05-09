use crate::objects::composition::{
    AtmosphereComposition, HydrosphereComposition, LithosphereComposition,
};
use uom::si::f64::{Acceleration, Angle, Area, Length, Mass, MassDensity, Pressure, Ratio, Time};

#[derive(Debug, Clone, PartialEq)]
pub struct PlanetaryBody {
    /// Planet radius
    pub radius: Length,
    /// Planet surface area
    pub surface_area: Area,
    /// Planet surface gravity (affects geological processes)
    pub surface_gravity: Acceleration,
    /// Planet mass (derived parameter)
    pub mass: Mass,
    /// Planet density (affects internal structure)
    pub density: MassDensity,

    /// Axial tilt (affects seasonal variation)
    pub axial_tilt: Angle,
    /// Rotation period (day length)
    pub rotation_period: Time,

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
    /// Atmospheric pressure at sea level
    pub atmospheric_pressure: Pressure,
    /// Normalized greenhouse effect (1.0 = Earth-like, 0.18 = Mars, ~15 = Venus)
    pub greenhouse_factor: Ratio,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Hydrosphere {
    /// Ocean Liquid
    pub composition: HydrosphereComposition,
    /// Ocean Coverage
    pub ocean_coverage: Ratio,
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
