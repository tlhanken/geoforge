//! Cosmology report serialization and file output.

use crate::layer_output::format_pre_tectonic_layers;
use crate::pipeline::{
    galaxy, galaxy_count, markers_near_focus, primary_ref, regions_near, solar_system_at,
};
use geoforge_types::cosmology::CosmologyLayers;
use geoforge_types::cosmology::{
    CosmologyPreset, GalaxyInstance, SolarSystem, StellarSystemMarker, StellarSystemRef,
};
use geoforge_types::coordinates::GalacticPoint;
use geoforge_types::planetary::PlanetaryBodies;
use geoforge_types::PipelineLayerId;
use serde::Serialize;
use std::fs;
use std::io;
use std::path::Path;

/// Export format.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExportFormat {
    /// `cosmology_seed<N>.json`
    Json,
    /// `cosmology_seed<N>.txt`
    Text,
    /// Both files.
    Both,
}

/// Full pre-tectonic run for inspection (regenerable from `seed` + `preset`).
#[derive(Debug, Serialize)]
pub struct CosmologyReport {
    /// Root seed value.
    pub seed: u64,
    /// Preset used for generation.
    pub preset: CosmologyPreset,
    /// Last pipeline layer included.
    pub pipeline_to: String,
    /// Load focus (ly).
    pub focus_ly: GalacticPoint,
    /// Per-galaxy summary.
    pub galaxies: Vec<GalaxyReport>,
    /// Primary system reference.
    pub primary_ref: StellarSystemRef,
    /// Stage 4 output when `pipeline_to` includes solar system.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub primary_solar_system: Option<SolarSystem>,
    /// Stage 5 output when `pipeline_to` includes planetary body.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub planetary_bodies: Option<PlanetaryBodies>,
}

/// One galaxy in the export report.
#[derive(Debug, Serialize)]
pub struct GalaxyReport {
    /// Galaxy instance metadata.
    pub galaxy: GalaxyInstance,
    /// Regions intersecting load sphere.
    pub regions_near_focus: usize,
    /// Stellar markers loaded.
    pub markers: Vec<StellarSystemMarker>,
}

/// Build a report from staged output.
#[must_use]
pub fn build_cosmology_report(
    layers: &CosmologyLayers,
    preset: CosmologyPreset,
    to_layer: PipelineLayerId,
) -> CosmologyReport {
    let g_count = galaxy_count(layers);
    let mut galaxies = Vec::new();

    for g in 0..g_count {
        let gal = galaxy(layers, g);
        let regions = regions_near(layers, g);
        let markers = markers_near_focus(layers, g);
        galaxies.push(GalaxyReport {
            galaxy: gal,
            regions_near_focus: regions.len(),
            markers,
        });
    }

    let pref = primary_ref(layers);
    let primary_solar_system = if to_layer >= PipelineLayerId::SolarSystem {
        Some(
            layers
                .solar_system
                .clone()
                .unwrap_or_else(|| solar_system_at(layers, pref)),
        )
    } else {
        None
    };

    let planetary_bodies = if to_layer >= PipelineLayerId::PlanetaryBody {
        layers.planetary.clone()
    } else {
        None
    };

    CosmologyReport {
        seed: layers.root().value(),
        preset,
        pipeline_to: to_layer.as_str().to_string(),
        focus_ly: layers.focus_ly(),
        galaxies,
        primary_ref: pref,
        primary_solar_system,
        planetary_bodies,
    }
}

/// Write report under `dir`.
pub fn write_cosmology_report(
    dir: &Path,
    layers: &CosmologyLayers,
    report: &CosmologyReport,
    format: ExportFormat,
) -> io::Result<Vec<std::path::PathBuf>> {
    fs::create_dir_all(dir)?;
    let stem = format!("cosmology_seed{}", report.seed);
    let mut written = Vec::new();

    if matches!(format, ExportFormat::Json | ExportFormat::Both) {
        let path = dir.join(format!("{stem}.json"));
        let json = serde_json::to_string_pretty(report)?;
        fs::write(&path, json)?;
        written.push(path);
    }

    if matches!(format, ExportFormat::Text | ExportFormat::Both) {
        let path = dir.join(format!("{stem}.txt"));
        let to_layer = PipelineLayerId::parse_str(&report.pipeline_to)
            .unwrap_or(PipelineLayerId::PlanetaryBody);
        fs::write(
            &path,
            format_cosmology_report_text_from_layers(layers, report.preset, to_layer),
        )?;
        written.push(path);
    }

    Ok(written)
}

/// Human-readable report body (layer-by-layer; use [`format_cosmology_report_text_from_layers`] when layers are available).
#[must_use]
pub fn format_cosmology_report_text(report: &CosmologyReport) -> String {
    let mut out = String::new();
    out.push_str("Geoforge V2 — pre-tectonic pipeline report\n");
    out.push_str(&format!("Seed: {}\n", report.seed));
    out.push_str(&format!("Preset: {:?}\n", report.preset));
    out.push_str(&format!("Pipeline to: {}\n", report.pipeline_to));
    out.push_str(&format!(
        "Focus (ly): {:.4}, {:.4}, {:.4}\n\n",
        report.focus_ly.x_ly, report.focus_ly.y_ly, report.focus_ly.z_ly
    ));
    out.push_str("(Re-run export from pipeline layers for full layer sections.)\n");
    out
}

/// Layer-declared text report from staged output (preferred for `.txt` export).
#[must_use]
pub fn format_cosmology_report_text_from_layers(
    layers: &CosmologyLayers,
    preset: CosmologyPreset,
    to_layer: PipelineLayerId,
) -> String {
    let mut out = String::new();
    out.push_str("Geoforge V2 — pre-tectonic pipeline report\n\n");
    out.push_str(&format_pre_tectonic_layers(layers, preset, to_layer, 1));
    out
}
