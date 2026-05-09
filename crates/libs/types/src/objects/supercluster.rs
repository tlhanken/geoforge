use uom::si::f64::{Length, Mass, Time};

#[derive(Debug, Clone, PartialEq)]
pub enum SuperclusterType {
    Regular,
    Irregular,
    Dwarf,
}
#[derive(Debug, Clone, PartialEq)]
pub struct Supercluster {
    pub supercluster_type: SuperclusterType,
    pub mass: Mass,
    pub radius: Length,
    pub age: Time,
}
