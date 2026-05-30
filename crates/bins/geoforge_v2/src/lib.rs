use clap::Parser;
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
            from_layer,
            to_layer,
            output_directory,
            output_file_format,
            cli_verbosity,
        } => {
            println!("GeoForge V2 — generation not yet wired to stage crates");
            if let Some(s) = seed {
                println!("Seed: {s}");
            } else {
                println!("Seed: (random — not implemented)");
            }
            println!(
                "Layers: {} → {}",
                from_layer.0.as_str(),
                to_layer.0.as_str()
            );
            println!("Output: {} ({output_file_format:?})", output_directory.display());
            println!("Verbosity: {cli_verbosity}");
        }
    }
}
