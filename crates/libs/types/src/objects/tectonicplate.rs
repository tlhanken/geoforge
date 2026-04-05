#[derive(Debug, Clone, PartialEq)]
pub enum TectonicPlateType {
    Continental,
    Oceanic,
    Mixed,
}

#[derive(Debug, Clone, PartialEq)]
pub struct TectonicPlate {
    pub age: f64,
    pub plate_type: TectonicPlateType,
}
