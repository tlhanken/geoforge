use super::Coordinates;
use uom::si::f64::Angle;

#[derive(Debug, Clone, PartialEq)]
pub struct PlanetaryCoordinates {
    pub lat: Angle,
    pub lon: Angle,
}

impl Coordinates for PlanetaryCoordinates {
    fn generate_random_coordinates() -> Self
    where
        Self: Sized,
    {
        todo!("Generate a new random coordinate");
    }
    /// Calculate the distance between two coordinates within the same coordinate frame of reference
    fn calculate_distance(&self, _other: &Self) -> f64 {
        todo!("Calculate the distance between two coordinates");
    }
}
