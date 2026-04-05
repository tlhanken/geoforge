use super::Coordinates;
use uom::si::f64::{Time, Length, Ratio, Angle};

#[derive(Debug, Clone, PartialEq)]
pub struct OrbitalCoordinates {
    /// Orbital period
    pub orbital_period: Time,
    /// Orbital eccentricity (0 = perfect circle, 0.99 = highly elliptical)
    pub orbital_eccentricity: Ratio,
    /// Semi-major axis (average orbital distance)
    pub semi_major_axis: Length,
    /// Perihelion distance (closest approach to star)
    pub perihelion: Length,
    /// Aphelion distance (farthest distance from star)
    pub aphelion: Length,
    /// Orbital inclination (relative to ecliptic)
    pub orbital_inclination: Angle,
    /// Current phase of the orbit [0.0, 1.0).
    /// 0.0 = Perihelion, 0.5 = Aphelion, wrapping back to 0.
    pub orbital_phase: Ratio,
}

impl Eq for OrbitalCoordinates {}
#[allow(clippy::derive_hash_xor_eq)]
impl std::hash::Hash for OrbitalCoordinates {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.orbital_period.value.to_bits().hash(state);
        self.orbital_eccentricity.value.to_bits().hash(state);
        self.semi_major_axis.value.to_bits().hash(state);
        self.perihelion.value.to_bits().hash(state);
        self.aphelion.value.to_bits().hash(state);
        self.orbital_inclination.value.to_bits().hash(state);
        self.orbital_phase.value.to_bits().hash(state);
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
