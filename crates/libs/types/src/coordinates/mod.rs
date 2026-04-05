pub mod universe;
pub mod supercluster;
pub mod galaxy;
pub mod solar;
pub mod planetary;

#[derive(Debug, Clone, PartialEq)]
pub struct XYZPosition {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}
