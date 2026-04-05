pub mod universe;
pub mod supercluster;
pub mod galaxy;
pub mod orbital;
pub mod planetary;

use uom::si::f64::Length;

#[derive(Debug, Clone, PartialEq)]
pub struct XYZPosition {
    pub x: Length,
    pub y: Length,
    pub z: Length,
}

impl Eq for XYZPosition {}
#[allow(clippy::derive_hash_xor_eq)]
impl std::hash::Hash for XYZPosition {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.x.value.to_bits().hash(state);
        self.y.value.to_bits().hash(state);
        self.z.value.to_bits().hash(state);
    }
}

pub trait Coordinates {
    fn generate_random_coordinates() -> Self where Self: Sized {
        todo!("Generate a new random coordinate");
    }
    /// Calculate the distance between two coordinates within the same coordinate frame of reference
    fn calculate_distance(&self, other: &Self) -> f64 {
        todo!("Calculate the distance between two coordinates");
    }
}