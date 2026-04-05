pub mod universe;
pub mod supercluster;
pub mod galaxy;
pub mod solarsystem;
pub mod planetary;

pub use universe::UniverseLayer;
pub use supercluster::SuperclusterLayer;
pub use galaxy::GalaxyLayer;
pub use solarsystem::SolarSystemLayer;
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
