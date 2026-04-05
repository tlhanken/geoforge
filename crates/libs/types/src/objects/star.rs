use uom::si::f64::{Ratio, Mass, Length, ThermodynamicTemperature, Time, Power};

#[derive(Debug, Clone, PartialEq)]
pub enum SpectralClass {
    O,
    B,
    A,
    F,
    G,
    K,
    M,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Star {
    pub stellar_luminosity: Power,
    pub mass: Mass,
    pub radius: Length,
    pub surface_temperature: ThermodynamicTemperature,
    pub age: Time,
    pub metallicity: Ratio,
    pub spectral_class: SpectralClass,
}
