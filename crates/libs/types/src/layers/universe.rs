use super::Seed;
use super::superclusters::SuperclusterSeed;
use crate::{coordinates::universe::UniverseCoordinates, objects::Universe};

#[derive(Debug, Clone, PartialEq)]
pub struct UniverseLayer {
    pub seed: UniverseSeed,
    pub universe: Universe,
}

#[derive(Debug, Clone, PartialEq)]
pub struct UniverseSeed(pub u64);

impl Seed for UniverseSeed {
    type Coordinate = UniverseCoordinates;
    type Subseed = SuperclusterSeed;
    fn generate_subseed(&self, _coordinates: &Self::Coordinate) -> Self::Subseed {
        todo!("use seed and XYZ coordinates to generate a new seed, share with others if possible")
    }
    fn generate_random_seed() -> Self {
        todo!("generate a U64 using a standard random library")
    }
}
