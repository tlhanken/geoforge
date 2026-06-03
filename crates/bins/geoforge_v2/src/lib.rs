use clap::Parser;
use geoforge_cosmology::CosmologyContext;
use geoforge_types::{PipelineLayerId, Seed};

pub mod cli;
pub mod export;
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
            export,
            no_export,
            cli_verbosity,
        } => run_generate(
            *seed,
            cosmology_scale.0,
            from_layer.0,
            to_layer.0,
            output_directory,
            export.into_export(),
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
    export_format: export::ExportFormat,
    no_export: bool,
    verbosity: u8,
) {
    let root = Seed::new(seed.unwrap_or(0x60F0_06E5_5EED));
    println!("GeoForge V2 — seed {}", root.value());
    println!("Cosmology preset: {:?}", preset);
    println!("Pipeline: {} → {}", from_layer.as_str(), to_layer.as_str());
    if verbosity > 0 {
        println!("Output dir (optional): {}", output_directory.display());
    }

    if from_layer > to_layer {
        eprintln!("Error: from_layer must be <= to_layer");
        std::process::exit(1);
    }

    if to_layer >= PipelineLayerId::Galaxy {
        run_cosmology(
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
        println!("Tectonics and below: not yet implemented in v2 pipeline.");
    }
}

fn run_cosmology(
    root: Seed,
    preset: geoforge_types::cosmology::CosmologyPreset,
    to_layer: PipelineLayerId,
    output_directory: &std::path::Path,
    export_format: export::ExportFormat,
    no_export: bool,
    verbosity: u8,
) {
    let ctx = CosmologyContext::new(root, preset);
    let g_count = ctx.galaxy_count();

    println!("\n=== Cosmology (regenerated from seed) ===");
    println!("Galaxies: {g_count}");

    for g in 0..g_count {
        let galaxy = ctx.galaxy(g);
        println!(
            "  [{}] {} — {:?}, R={:.0} ly",
            galaxy.index,
            galaxy.label,
            galaxy.profile.morphology,
            galaxy.profile.radius_ly
        );

        if to_layer < PipelineLayerId::Galaxy {
            continue;
        }

        let regions = ctx.regions_near(g);
        if verbosity >= 1 {
            println!("    Regions near focus: {}", regions.len());
        }

        let markers = ctx.markers_near_focus(g);
        println!("    Stellar markers loaded: {}", markers.len());

        if verbosity >= 2 {
            for m in markers.iter().take(10) {
                println!(
                    "      #{}{} @ ({:.1}, {:.1}, {:.1}) ly — {:?} {:?}",
                    m.index_in_region,
                    if m.is_primary { " [primary]" } else { "" },
                    m.position_ly.x_ly,
                    m.position_ly.y_ly,
                    m.position_ly.z_ly,
                    m.phenotype.multiplicity,
                    m.phenotype.primary_class
                );
            }
            if markers.len() > 10 {
                println!("      ...");
            }
        }
    }

    if to_layer >= PipelineLayerId::SolarSystem {
        let primary = ctx.primary_ref();
        let system = ctx.solar_system(primary);
        println!("\n=== Primary solar system (zoom) ===");
        println!(
            "  Position: ({:.1}, {:.1}, {:.1}) ly",
            system.barycenter_ly.x_ly,
            system.barycenter_ly.y_ly,
            system.barycenter_ly.z_ly
        );
        println!("  Stars: {}", system.stars.len());
        for (i, star) in system.stars.iter().enumerate() {
            println!(
                "    Star {i}: {:?}, {:.2} L☉, offset ({:.1}, {:.1}, {:.1}) AU",
                star.spectral_class,
                star.luminosity_solar,
                star.position_au.x_au,
                star.position_au.y_au,
                star.position_au.z_au
            );
        }
        println!("  Planet slots: {}", system.planet_slots.len());
        for slot in &system.planet_slots {
            println!(
                "    Slot {}: {:?} @ {:.2} AU",
                slot.index, slot.mass_class, slot.orbit.semi_major_axis_au
            );
        }
        if !system.belts.is_empty() {
            println!("  Belts: {}", system.belts.len());
        }
    }

    if !no_export && to_layer >= PipelineLayerId::Galaxy {
        let report = export::build_report(&ctx, root, preset, to_layer);
        match export::write_report(output_directory, &report, export_format) {
            Ok(paths) => {
                println!("\n=== Exported ===");
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
