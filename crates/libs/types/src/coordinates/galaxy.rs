use super::Coordinates;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct GalaxyCoordinates(pub crate::coordinates::XYZPosition);

impl Coordinates for GalaxyCoordinates {
    fn generate_random_coordinates() -> Self where Self: Sized {
        todo!("Generate a new random coordinate");
    }
    /// Calculate the distance between two coordinates within the same coordinate frame of reference
    fn calculate_distance(&self, other: &Self) -> f64 {
        todo!("Calculate the distance between two coordinates");
    }
}
