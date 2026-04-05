use std::collections::HashMap;
use super::solarsystems::SolarSystemLayer;
use crate::coordinates::orbital::OrbitalCoordinates;
use crate::objects::planetary_body::PlanetaryBody;
use super::Seed;

#[derive(Debug, Clone, PartialEq)]
pub struct PlanetaryLayer {
    pub seed: PlanetarySeed,
    pub solar_system_layer: SolarSystemLayer,
    pub planets: HashMap<OrbitalCoordinates, PlanetaryBody>, 
}

#[derive(Debug, Clone, PartialEq)]
pub struct PlanetarySeed(pub u64);

impl Seed for PlanetarySeed {
    type Coordinate = crate::coordinates::planetary::PlanetaryCoordinates;
    type Subseed = Self;
    fn generate_subseed(&self, coordinates: &Self::Coordinate) -> Self::Subseed {
        todo!("Planetary seed generation is terminal for now")
    }
    fn generate_random_seed() -> Self {
        todo!("Generate a U64 using a standard random library")
    }
}
