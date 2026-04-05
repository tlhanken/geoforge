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
    /// Stellar luminosity relative to Sun (1.0 = solar luminosity)
    pub stellar_luminosity: f64,
    pub mass: f64, // solar masses
    pub radius: f64, // solar radii
    pub temperature: f64, // kelvin
    pub age: f64, // years
    pub metallicity: f64, // solar metallicity
    pub spectral_class: SpectralClass,
}
