//! Export pipeline outputs to inspectable JSON and text files.

#![deny(clippy::todo, clippy::unimplemented)]

mod layer_output;
mod pipeline;
mod report;

pub use geoforge_types::cosmology::CosmologyLayers;
pub use pipeline::{
    galaxy, galaxy_count, markers_in_region, markers_near_focus, primary_marker, primary_ref,
    primary_solar_system, regions_near, run_cosmology_layers, run_pre_tectonic_layers,
    run_pre_tectonic_layers_with_config, solar_system_at, with_focus,
};
pub use layer_output::format_pre_tectonic_layers;
pub use report::{
    build_cosmology_report, format_cosmology_report_text, write_cosmology_report, CosmologyReport,
    ExportFormat, GalaxyReport,
};

#[cfg(test)]
mod lib_tests;
