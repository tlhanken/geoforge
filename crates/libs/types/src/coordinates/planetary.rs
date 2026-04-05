use super::Coordinates;

#[derive(Debug, Clone, PartialEq)]
pub struct PlanetaryCoordinates {
    pub lat: f64,
    pub lon: f64,
}

impl Coordinates for PlanetaryCoordinates {
    fn generate_random_coordinates() -> Self where Self: Sized {
        todo!("Generate a new random coordinate");
    }
    /// Calculate the distance between two coordinates within the same coordinate frame of reference
    fn calculate_distance(&self, other: &Self) -> f64 {
        todo!("Calculate the distance between two coordinates");
    }
}
