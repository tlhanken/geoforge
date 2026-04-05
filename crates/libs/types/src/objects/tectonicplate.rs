use uom::si::f64::Time;

#[derive(Debug, Clone, PartialEq)]
pub enum TectonicPlateType {
    Continental,
    Oceanic,
    Mixed,
}

#[derive(Debug, Clone, PartialEq)]
pub struct TectonicPlate {
    pub age: Time,
    pub plate_type: TectonicPlateType,
}
