use super::Seed;
use super::solarsystems::SolarSystemSeed;
use super::superclusters::SuperclusterLayer;
use crate::coordinates::galaxy::GalaxyCoordinates;
use crate::coordinates::supercluster::SuperclusterCoordinates;
use crate::objects::Galaxy;
use std::collections::HashMap;

#[derive(Debug, Clone, PartialEq)]
pub struct GalaxyLayer {
    pub seed: GalaxySeed,
    pub supercluster_layer: SuperclusterLayer,
    pub galaxies: HashMap<SuperclusterCoordinates, Galaxy>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct GalaxySeed(pub u64);

impl Seed for GalaxySeed {
    type Coordinate = GalaxyCoordinates;
    type Subseed = SolarSystemSeed;
    fn generate_subseed(&self, _coordinates: &Self::Coordinate) -> Self::Subseed {
        todo!("Use seed and parent coordinates to generate a new subseed")
    }
    fn generate_random_seed() -> Self {
        todo!("Generate a U64 using a standard random library")
    }
}
