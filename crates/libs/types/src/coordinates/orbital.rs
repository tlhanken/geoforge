use super::Coordinates;

#[derive(Debug, Clone, PartialEq)]
pub struct OrbitalCoordinates {
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
    /// Current phase of the orbit [0.0, 1.0).
    /// 0.0 = Perihelion, 0.5 = Aphelion, wrapping back to 0.
    pub orbital_phase: f64,
}

impl Eq for OrbitalCoordinates {}
#[allow(clippy::derive_hash_xor_eq)]
impl std::hash::Hash for OrbitalCoordinates {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.orbital_period_days.to_bits().hash(state);
        self.orbital_eccentricity.to_bits().hash(state);
        self.semi_major_axis_au.to_bits().hash(state);
        self.perihelion_au.to_bits().hash(state);
        self.aphelion_au.to_bits().hash(state);
        self.orbital_inclination_degrees.to_bits().hash(state);
        self.orbital_phase.to_bits().hash(state);
    }
}

impl Coordinates for OrbitalCoordinates {
    fn generate_random_coordinates() -> Self where Self: Sized {
        todo!("Generate a new random coordinate");
    }
    /// Calculate the distance between two coordinates within the same coordinate frame of reference
    fn calculate_distance(&self, other: &Self) -> f64 {
        todo!("Calculate the distance between two coordinates");
    }
}
