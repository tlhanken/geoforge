use clap::Parser;
use geoforge_export::{
    build_cosmology_report, format_pre_tectonic_layers, run_pre_tectonic_layers,
    write_cosmology_report, ExportFormat,
};
use geoforge_types::{PipelineLayerId, Seed};

pub mod cli;
use cli::Commands;

/// GeoForge V2 — procedural world generation CLI
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

/// Run the CLI (called from `main`).
pub fn run() {
    let cli = Cli::parse();

    match &cli.command {
        Commands::Generate {
            seed,
            cosmology_scale,
            from_layer,
            to_layer,
            output_directory,
            format,
            no_export,
            cli_verbosity,
        } => run_generate(
            *seed,
            cosmology_scale.0,
            from_layer.0,
            to_layer.0,
            output_directory,
            format.into_export(),
            *no_export,
            *cli_verbosity,
        ),
    }
}

fn run_generate(
    seed: Option<u64>,
    preset: geoforge_types::cosmology::CosmologyPreset,
    from_layer: PipelineLayerId,
    to_layer: PipelineLayerId,
    output_directory: &std::path::Path,
    export_format: ExportFormat,
    no_export: bool,
    verbosity: u8,
) {
    let root = Seed::new(seed.unwrap_or(0x60F0_06E5_5EED));
    println!("geoforge_v2 {}", env!("CARGO_PKG_VERSION"));
    println!("Seed: {}", root.value());
    println!("Preset: {:?}", preset);
    println!(
        "Pipeline layers: {} → {}",
        from_layer.as_str(),
        to_layer.as_str()
    );
    if !no_export {
        println!(
            "Export: {} → {}",
            output_directory.display(),
            format_label(export_format)
        );
    }

    if from_layer > to_layer {
        eprintln!("Error: from_layer must be <= to_layer");
        std::process::exit(1);
    }

    if to_layer >= PipelineLayerId::Universe {
        run_pre_tectonic(
            root,
            preset,
            to_layer,
            output_directory,
            export_format,
            no_export,
            verbosity,
        );
    }

    if to_layer >= PipelineLayerId::Tectonics {
        println!("\n── Layer 6 · tectonics (stage6_tectonics) ──");
        println!("  (not implemented in v2 yet)");
    }
}

fn format_label(f: ExportFormat) -> &'static str {
    match f {
        ExportFormat::Json => "json",
        ExportFormat::Text => "text",
        ExportFormat::Both => "json + text",
    }
}

fn run_pre_tectonic(
    root: Seed,
    preset: geoforge_types::cosmology::CosmologyPreset,
    to_layer: PipelineLayerId,
    output_directory: &std::path::Path,
    export_format: ExportFormat,
    no_export: bool,
    verbosity: u8,
) {
    let layers = run_pre_tectonic_layers(root, preset, to_layer);

    print!("{}", format_pre_tectonic_layers(&layers, preset, to_layer, verbosity));

    if !no_export && to_layer >= PipelineLayerId::Galaxy {
        let report = build_cosmology_report(&layers, preset, to_layer);
        match write_cosmology_report(output_directory, &layers, &report, export_format) {
            Ok(paths) => {
                println!("\n── Export ──");
                for p in paths {
                    println!("  {}", p.display());
                }
            }
            Err(e) => {
                eprintln!("Export failed: {e}");
                std::process::exit(1);
            }
        }
    }
}
