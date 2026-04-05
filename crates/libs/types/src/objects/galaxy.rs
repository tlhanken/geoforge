use uom::si::f64::{Mass, Length, Time};

#[derive(Debug, Clone, PartialEq)]
pub enum GalaxyType {
    Elliptical,
    Spiral,
    Irregular,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Galaxy {
    pub galaxy_type: GalaxyType,
    pub mass: Mass,
    pub radius: Length,
    pub age: Time,
}
