#[derive(Debug, Clone, PartialEq)]
pub enum SuperclusterType {
    Regular,
    Irregular,
    Dwarf,
}
#[derive(Debug, Clone, PartialEq)]
pub struct Supercluster {
    pub supercluster_type: SuperclusterType,
    pub mass: f64,
    pub radius: f64,
    pub age: f64,
}
