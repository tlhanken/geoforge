//! Write generation results to inspectable files (JSON / text).

use geoforge_cosmology::CosmologyContext;
use geoforge_types::cosmology::{
    CosmologyPreset, GalaxyInstance, SolarSystem, StellarSystemMarker, StellarSystemRef,
};
use geoforge_types::coordinates::GalacticPoint;
use geoforge_types::cosmology::RegionId;
use geoforge_types::{PipelineLayerId, Seed};
use serde::Serialize;
use std::fs;
use std::io;
use std::path::Path;

/// Export format for CLI.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExportFormat {
    Json,
    Text,
    Both,
}

/// Full cosmology run captured for inspection (regenerable from `seed` + `preset`).
#[derive(Debug, Serialize)]
pub struct CosmologyReport {
    pub seed: u64,
    pub preset: CosmologyPreset,
    pub pipeline_to: String,
    pub focus_ly: GalacticPoint,
    pub galaxies: Vec<GalaxyReport>,
    pub primary_ref: StellarSystemRef,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub primary_solar_system: Option<SolarSystem>,
}

#[derive(Debug, Serialize)]
pub struct GalaxyReport {
    pub galaxy: GalaxyInstance,
    pub regions_near_focus: usize,
    pub markers: Vec<StellarSystemMarker>,
}

/// Build a report from the current context.
#[must_use]
pub fn build_report(
    ctx: &CosmologyContext,
    root: Seed,
    preset: CosmologyPreset,
    to_layer: PipelineLayerId,
) -> CosmologyReport {
    let g_count = ctx.galaxy_count();
    let mut galaxies = Vec::new();

    for g in 0..g_count {
        let galaxy = ctx.galaxy(g);
        let regions = ctx.regions_near(g);
        let markers = ctx.markers_near_focus(g);
        galaxies.push(GalaxyReport {
            galaxy,
            regions_near_focus: regions.len(),
            markers,
        });
    }

    let primary_ref = ctx.primary_ref();
    let primary_solar_system = if to_layer >= PipelineLayerId::SolarSystem {
        Some(ctx.solar_system(primary_ref))
    } else {
        None
    };

    CosmologyReport {
        seed: root.value(),
        preset,
        pipeline_to: to_layer.as_str().to_string(),
        focus_ly: ctx.focus_ly,
        galaxies,
        primary_ref,
        primary_solar_system,
    }
}

/// Write report to `dir` using the given format(s).
pub fn write_report(
    dir: &Path,
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
        fs::write(&path, format_report_text(report))?;
        written.push(path);
    }

    Ok(written)
}

fn format_report_text(report: &CosmologyReport) -> String {
    let mut out = String::new();
    out.push_str(&format!("Geoforge V2 cosmology report\n"));
    out.push_str(&format!("Seed: {}\n", report.seed));
    out.push_str(&format!("Preset: {:?}\n", report.preset));
    out.push_str(&format!("Pipeline to: {}\n", report.pipeline_to));
    out.push_str(&format!(
        "Focus (ly): {:.4}, {:.4}, {:.4}\n\n",
        report.focus_ly.x_ly, report.focus_ly.y_ly, report.focus_ly.z_ly
    ));

    for g in &report.galaxies {
        out.push_str(&format!(
            "=== Galaxy {}: {} ===\n",
            g.galaxy.index, g.galaxy.label
        ));
        out.push_str(&format!(
            "  Morphology: {:?}\n  Radius: {:.1} ly\n  Regions near focus: {}\n  Markers: {}\n\n",
            g.galaxy.profile.morphology,
            g.galaxy.profile.radius_ly,
            g.regions_near_focus,
            g.markers.len()
        ));
        for m in &g.markers {
            out.push_str(&format!(
                "  - #{} galaxy={} region={} @ ({:.2}, {:.2}, {:.2}) ly  {:?}/{:?}{}\n",
                m.index_in_region,
                m.galaxy_index,
                region_summary(m.region),
                m.position_ly.x_ly,
                m.position_ly.y_ly,
                m.position_ly.z_ly,
                m.phenotype.multiplicity,
                m.phenotype.primary_class,
                if m.is_primary { " [PRIMARY]" } else { "" }
            ));
        }
        out.push('\n');
    }

    out.push_str("=== Primary reference ===\n");
    out.push_str(&format!(
        "  galaxy={} region={} index={}\n\n",
        report.primary_ref.galaxy_index,
        region_summary(report.primary_ref.region),
        report.primary_ref.index_in_region
    ));

    if let Some(sys) = &report.primary_solar_system {
        out.push_str("=== Primary solar system ===\n");
        out.push_str(&format!(
            "  Barycenter (ly): {:.4}, {:.4}, {:.4}\n",
            sys.barycenter_ly.x_ly, sys.barycenter_ly.y_ly, sys.barycenter_ly.z_ly
        ));
        out.push_str(&format!("  Stars: {}\n", sys.stars.len()));
        for (i, star) in sys.stars.iter().enumerate() {
            out.push_str(&format!(
                "    [{i}] {:?} {:.3} L☉  offset AU ({:.2}, {:.2}, {:.2})\n",
                star.spectral_class,
                star.luminosity_solar,
                star.position_au.x_au,
                star.position_au.y_au,
                star.position_au.z_au
            ));
        }
        out.push_str(&format!("  Planet slots: {}\n", sys.planet_slots.len()));
        for slot in &sys.planet_slots {
            out.push_str(&format!(
                "    [{}] {:?}  a={:.3} AU  period={:.1} d  e={:.3}\n",
                slot.index,
                slot.mass_class,
                slot.orbit.semi_major_axis_au,
                slot.orbit.period_days,
                slot.orbit.eccentricity
            ));
        }
        if !sys.belts.is_empty() {
            out.push_str(&format!("  Belts: {}\n", sys.belts.len()));
            for b in &sys.belts {
                out.push_str(&format!(
                    "    [{}] {:.2}–{:.2} AU\n",
                    b.index, b.inner_au, b.outer_au
                ));
            }
        }
    }

    out
}

fn region_summary(r: RegionId) -> String {
    match r {
        RegionId::Cylindrical { ring, wedge, layer } => {
            format!("cyl(r={ring},w={wedge},z={layer})")
        }
        RegionId::Spherical {
            shell,
            theta_wedge,
            phi_wedge,
        } => format!("sph(s={shell},θ={theta_wedge},φ={phi_wedge})"),
    }
}
