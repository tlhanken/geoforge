use clap::{Subcommand, ValueEnum};
use geoforge_types::cosmology::CosmologyPreset;
use geoforge_types::PipelineLayerId;
use std::path::PathBuf;
use std::str::FromStr;

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Generate a new procedural world from seed (no disk persistence required).
    Generate {
        /// Seed for procedural generation
        #[arg(short, long)]
        seed: Option<u64>,

        /// Cosmology scale: minimal (1 system), rich (local region), expansive
        #[arg(long, default_value = "minimal")]
        cosmology_scale: CosmologyScaleArg,

        /// First pipeline layer (inclusive)
        #[arg(long, default_value = "galaxy")]
        from_layer: PipelineLayerArg,

        /// Last pipeline layer (inclusive)
        #[arg(long, default_value = "solar_system")]
        to_layer: PipelineLayerArg,

        /// Output directory for export files
        #[arg(short = 'd', long, default_value = "outputs")]
        output_directory: PathBuf,

        /// Write inspectable export files: json, text, or both (default: json)
        #[arg(long, value_enum, default_value = "json")]
        export: ExportFormatArg,

        /// Skip export files (stdout summary only)
        #[arg(long)]
        no_export: bool,

        /// CLI verbosity (repeat `-v` for more)
        #[arg(short = 'v', long, action = clap::ArgAction::Count)]
        cli_verbosity: u8,
    },
}

/// Export format for inspection.
#[derive(Clone, Copy, Debug, ValueEnum)]
pub enum ExportFormatArg {
    /// Pretty-printed JSON (`cosmology_seed<N>.json`).
    Json,
    /// Human-readable text (`cosmology_seed<N>.txt`).
    Text,
    /// Both JSON and text.
    Both,
}

impl ExportFormatArg {
    #[must_use]
    pub fn into_export(self) -> crate::export::ExportFormat {
        match self {
            Self::Json => crate::export::ExportFormat::Json,
            Self::Text => crate::export::ExportFormat::Text,
            Self::Both => crate::export::ExportFormat::Both,
        }
    }
}

/// CLI wrapper for [`CosmologyPreset`].
#[derive(Clone, Debug)]
pub struct CosmologyScaleArg(pub CosmologyPreset);

impl ValueEnum for CosmologyScaleArg {
    fn value_variants<'a>() -> &'a [Self] {
        static VARIANTS: [CosmologyScaleArg; 3] = [
            CosmologyScaleArg(CosmologyPreset::Minimal),
            CosmologyScaleArg(CosmologyPreset::Rich),
            CosmologyScaleArg(CosmologyPreset::Expansive),
        ];
        &VARIANTS
    }

    fn to_possible_value(&self) -> Option<clap::builder::PossibleValue> {
        let s = match self.0 {
            CosmologyPreset::Minimal => "minimal",
            CosmologyPreset::Rich => "rich",
            CosmologyPreset::Expansive => "expansive",
        };
        Some(clap::builder::PossibleValue::new(s))
    }
}

/// CLI wrapper for [`PipelineLayerId`].
#[derive(Clone, Debug)]
pub struct PipelineLayerArg(pub PipelineLayerId);

impl ValueEnum for PipelineLayerArg {
    fn value_variants<'a>() -> &'a [Self] {
        static VARIANTS: [PipelineLayerArg; 12] = [
            PipelineLayerArg(PipelineLayerId::Universe),
            PipelineLayerArg(PipelineLayerId::Supercluster),
            PipelineLayerArg(PipelineLayerId::Galaxy),
            PipelineLayerArg(PipelineLayerId::SolarSystem),
            PipelineLayerArg(PipelineLayerId::PlanetaryBody),
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
