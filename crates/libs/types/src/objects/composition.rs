#[derive(Debug, Clone, PartialEq)]
pub enum AtmosphereComposition {
    Nitrogen,
    Oxygen,
    CarbonDioxide,
    Methane,
    Ammonia,
    Other,
}

#[derive(Debug, Clone, PartialEq)]
pub enum HydrosphereComposition {
    Water,
    Methane,
    Ammonia,
    Other,
}

#[derive(Debug, Clone, PartialEq)]
pub enum LithosphereComposition {
    Rock,
    Metal,
    Ice,
    Other,
}