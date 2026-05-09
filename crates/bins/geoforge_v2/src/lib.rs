use clap::Parser;
pub mod cli;
use cli::Commands;

/// GeoForge V2 - Procedural world generation CLI
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

pub fn run() {
    let cli = Cli::parse();

    match &cli.command {
        Commands::Generate {
            seed,
            from_layer,
            to_layer,
            output_directory,
            output_file_format,
            cli_verbosity,
        } => {
            println!("Initializing GeoForge V2 generation...");
            if let Some(s) = seed {
                println!("Using seed: {}", s);
            } else {
                println!("Using random seed");
            }
            println!("Layer scope: {:?} to {:?}", from_layer, to_layer);
            println!("Output directory: {}", output_directory.display());
            println!("Output format: {:?}", output_file_format);
            println!("Verbosity: {}", cli_verbosity);
        }
    }
}
