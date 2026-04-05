#[derive(Debug, Clone, PartialEq)]
pub enum GalaxyType {
    Elliptical,
    Spiral,
    Irregular,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Galaxy {
    pub galaxy_type: GalaxyType,
    pub mass: f64,
    pub radius: f64,
    pub age: f64,
}
