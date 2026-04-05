use clap::{Subcommand, ValueEnum};
use std::path::PathBuf;

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Generate a new procedural world
    Generate {
        /// Seed for the random number generator
        #[arg(short, long)]
        seed: Option<u64>,

        // Begin Layer
        #[arg(long, default_value = "planetary-body")]
        from_layer: Layers,

        // End Layer
        #[arg(long, default_value = "political-entities")]
        to_layer: Layers,

        /// Output directory
        #[arg(short = 'd', long, default_value = "outputs")]
        output_directory: PathBuf,

        /// Output type
        #[arg(short = 'f', long, default_value = "sqlite")]
        output_file_format: FileFormat,

        /// CLI Verbosity
        #[arg(short = 'v', long, action = clap::ArgAction::Count)]
        cli_verbosity: u8,
    },
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

#[derive(ValueEnum, Clone, Debug)]
pub enum FileFormat {
    Bin,
    Sqlite,
    Png,
    BinAndPng,
    SqliteAndPng,
    Geotiff,
}
