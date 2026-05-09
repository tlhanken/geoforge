pub mod galaxies;
pub mod planetary;
pub mod solarsystems;
pub mod superclusters;
pub mod universe;

use clap::ValueEnum;
pub use galaxies::GalaxyLayer;
pub use planetary::PlanetaryLayer;
pub use solarsystems::SolarSystemLayer;
pub use superclusters::SuperclusterLayer;
pub use universe::UniverseLayer;

pub trait Seed: Clone + PartialEq {
    type Coordinate: crate::coordinates::Coordinates;
    type Subseed: Seed;
    fn generate_subseed(&self, _coordinates: &Self::Coordinate) -> Self::Subseed {
        todo!("Generate a new seed from a previous seed and coordinates");
    }
    fn generate_random_seed() -> Self {
        todo!("Generate a new random seed");
    }
}

#[derive(ValueEnum, Clone, Debug)]
pub enum Layers {
    Universe,
    Supercluster,
    Galaxy,
    SolarSystem,
    PlanetaryBody,
    Tectonics,
    GeologicDomains,
    Heightmap,
    OceanCoverage,
    PrevailingWinds,
    Temperature,
    Precipitation,
    Watersheds,
    RiversAndLakes,
    Biomes,
    Resources,
    NaturalHazards,
    Settlements,
    TransportPaths,
    PoliticalEntities,
}
