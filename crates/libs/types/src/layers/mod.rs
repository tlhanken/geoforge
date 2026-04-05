pub mod universe;
pub mod superclusters;
pub mod galaxies;
pub mod solarsystems;
pub mod planetary;

pub use universe::UniverseLayer;
pub use superclusters::SuperclusterLayer;
pub use galaxies::GalaxyLayer;
pub use solarsystems::SolarSystemLayer;
pub use planetary::PlanetaryLayer;
use clap::ValueEnum;

pub trait Seed: Clone + PartialEq {
    type Coordinate: crate::coordinates::Coordinates;
    type Subseed: Seed;
    fn generate_subseed(&self, coordinates: &Self::Coordinate) -> Self::Subseed {
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
