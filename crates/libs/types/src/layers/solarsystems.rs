use super::Seed;
use super::galaxies::GalaxyLayer;
use super::planetary::PlanetarySeed;
use crate::coordinates::galaxy::GalaxyCoordinates;
use crate::coordinates::orbital::OrbitalCoordinates;
use crate::objects::solarsystem::SolarSystem;
use std::collections::HashMap;

#[derive(Debug, Clone, PartialEq)]
pub struct SolarSystemLayer {
    pub seed: SolarSystemSeed,
    pub galaxy_layer: GalaxyLayer,
    pub solar_systems: HashMap<GalaxyCoordinates, SolarSystem>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SolarSystemSeed(pub u64);

impl Seed for SolarSystemSeed {
    type Coordinate = OrbitalCoordinates;
    type Subseed = PlanetarySeed;
    fn generate_subseed(&self, _coordinates: &Self::Coordinate) -> Self::Subseed {
        todo!("Use seed and parent coordinates to generate a new subseed")
    }
    fn generate_random_seed() -> Self {
        todo!("Generate a U64 using a standard random library")
    }
}
