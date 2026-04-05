pub mod universe;
pub mod supercluster;
pub mod galaxy;

pub use universe::UniverseLayer;
pub use supercluster::SuperclusterLayer;
pub use galaxy::GalaxyLayer;
use clap::ValueEnum;

struct Seed(u64);

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
