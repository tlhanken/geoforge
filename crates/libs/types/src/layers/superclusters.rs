use super::Seed;
use super::galaxies::GalaxySeed;
use super::universe::UniverseLayer;
use crate::coordinates::supercluster::SuperclusterCoordinates;
use crate::coordinates::universe::UniverseCoordinates;
use crate::objects::Supercluster;
use std::collections::HashMap;

#[derive(Debug, Clone, PartialEq)]
pub struct SuperclusterLayer {
    pub seed: SuperclusterSeed,
    pub universe_layer: UniverseLayer,
    pub superclusters: HashMap<UniverseCoordinates, Supercluster>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SuperclusterSeed(pub u64);

impl Seed for SuperclusterSeed {
    type Coordinate = SuperclusterCoordinates;
    type Subseed = GalaxySeed;
    fn generate_subseed(&self, _coordinates: &Self::Coordinate) -> Self::Subseed {
        todo!("Use seed and parent coordinates to generate a new subseed")
    }
    fn generate_random_seed() -> Self {
        todo!("Generate a U64 using a standard random library")
    }
}
