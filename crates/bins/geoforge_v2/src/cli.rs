use clap::{Subcommand, ValueEnum};
use geoforge_types::PipelineLayerId;
use std::path::PathBuf;
use std::str::FromStr;

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Generate a new procedural world
    Generate {
        /// Seed for the random number generator
        #[arg(short, long)]
        seed: Option<u64>,

        /// First pipeline layer to generate (inclusive)
        #[arg(long, default_value = "tectonics")]
        from_layer: PipelineLayerArg,

        /// Last pipeline layer to generate (inclusive)
        #[arg(long, default_value = "geology")]
        to_layer: PipelineLayerArg,

        /// Output directory
        #[arg(short = 'd', long, default_value = "outputs")]
        output_directory: PathBuf,

        /// Output type
        #[arg(short = 'f', long, default_value = "bin")]
        output_file_format: FileFormat,

        /// CLI verbosity (repeat `-v` for more)
        #[arg(short = 'v', long, action = clap::ArgAction::Count)]
        cli_verbosity: u8,
    },
}

/// CLI wrapper for [`PipelineLayerId`] with clap parsing.
#[derive(Clone, Debug)]
pub struct PipelineLayerArg(pub PipelineLayerId);

impl ValueEnum for PipelineLayerArg {
    fn value_variants<'a>() -> &'a [Self] {
        static VARIANTS: [PipelineLayerArg; 8] = [
            PipelineLayerArg(PipelineLayerId::PlanetarySystem),
            PipelineLayerArg(PipelineLayerId::Tectonics),
            PipelineLayerArg(PipelineLayerId::Geology),
            PipelineLayerArg(PipelineLayerId::Elevation),
            PipelineLayerArg(PipelineLayerId::Climate),
            PipelineLayerArg(PipelineLayerId::Biomes),
            PipelineLayerArg(PipelineLayerId::Hydrology),
            PipelineLayerArg(PipelineLayerId::Resources),
        ];
        &VARIANTS
    }

    fn to_possible_value(&self) -> Option<clap::builder::PossibleValue> {
        Some(clap::builder::PossibleValue::new(self.0.as_str()))
    }
}

impl FromStr for PipelineLayerArg {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        PipelineLayerId::parse_str(s)
            .map(Self)
            .ok_or_else(|| format!("unknown layer {s:?}"))
    }
}

#[derive(ValueEnum, Clone, Debug)]
pub enum FileFormat {
    Bin,
    Png,
    Geotiff,
}
