//! Layer-by-layer human-readable pipeline output.

use crate::pipeline::{galaxy, galaxy_count, markers_near_focus, regions_near};
use geoforge_types::cosmology::CosmologyPreset;
use geoforge_types::cosmology::{CosmologyLayers, RegionId, SolarSystem};
use geoforge_types::planetary::PlanetaryBodies;
use geoforge_types::PipelineLayerId;
use std::fmt::Write as _;

/// Format the full pre-tectonic pipeline as declared layers and their outputs.
#[must_use]
pub fn format_pre_tectonic_layers(
    layers: &CosmologyLayers,
    preset: CosmologyPreset,
    to_layer: PipelineLayerId,
    verbosity: u8,
) -> String {
    let mut out = String::new();
    writeln!(out, "Pre-tectonic pipeline (seed {})", layers.root().value()).ok();
    writeln!(out, "Preset: {preset:?}  |  Generated through: {}", to_layer.as_str()).ok();
    writeln!(out).ok();

    format_layer_universe(&mut out, layers, verbosity);
    format_layer_supercluster(&mut out, layers, verbosity);
    format_layer_galaxy(&mut out, layers, to_layer, verbosity);
    format_layer_solar_system(&mut out, layers, to_layer, verbosity);
    format_layer_planetary_body(&mut out, layers, to_layer, verbosity);

    out
}

fn layer_heading(out: &mut String, stage: u8, layer: PipelineLayerId, crate_name: &str, type_name: &str) {
    writeln!(
        out,
        "── Layer {stage} · {} ({crate_name} → {type_name}) ──",
        layer.as_str()
    )
    .ok();
}

fn format_layer_universe(out: &mut String, layers: &CosmologyLayers, verbosity: u8) {
    layer_heading(out, 1, PipelineLayerId::Universe, "stage1_universe", "Universe");
    let u = &layers.universe;
    let c = &u.config;
    writeln!(out, "  root seed: {}", u.root.value()).ok();
    writeln!(out, "  label: {}", u.label).ok();
    writeln!(out, "  H₀: {:.1} km/s/Mpc", u.hubble_km_s_mpc).ok();
    writeln!(out, "  Ω_m: {:.3}", u.omega_matter).ok();
    if verbosity >= 1 {
        writeln!(out, "  config.preset: {:?}", c.preset).ok();
        writeln!(out, "  config.galaxy_count: {}", c.galaxy_count).ok();
        writeln!(out, "  config.load_radius_ly: {:.1}", c.load_radius_ly).ok();
        writeln!(
            out,
            "  config.max_systems_per_region: {}",
            c.max_systems_per_region
        )
        .ok();
        writeln!(
            out,
            "  config.max_markers_per_query: {}",
            c.max_markers_per_query
        )
        .ok();
    }
    writeln!(out).ok();
}

fn format_layer_supercluster(out: &mut String, layers: &CosmologyLayers, verbosity: u8) {
    layer_heading(
        out,
        2,
        PipelineLayerId::Supercluster,
        "stage2_supercluster",
        "Supercluster",
    );
    let s = &layers.supercluster;
    writeln!(out, "  label: {}", s.label).ok();
    writeln!(out, "  extent_ly: {:.1}", s.extent_ly).ok();
    writeln!(
        out,
        "  focus_ly: ({:.4}, {:.4}, {:.4})",
        s.focus_ly.x_ly, s.focus_ly.y_ly, s.focus_ly.z_ly
    )
    .ok();
    if verbosity >= 1 {
        writeln!(
            out,
            "  origin_ly: ({:.4}, {:.4}, {:.4})",
            s.origin_ly.x_ly, s.origin_ly.y_ly, s.origin_ly.z_ly
        )
        .ok();
        writeln!(out, "  root seed: {}", s.root.value()).ok();
    }
    writeln!(out).ok();
}

fn format_layer_galaxy(
    out: &mut String,
    layers: &CosmologyLayers,
    to_layer: PipelineLayerId,
    verbosity: u8,
) {
    layer_heading(
        out,
        3,
        PipelineLayerId::Galaxy,
        "stage3_galaxy",
        "GalaxyCatalog",
    );
    if to_layer < PipelineLayerId::Galaxy {
        writeln!(out, "  (not requested — pipeline stopped before galaxy)").ok();
        writeln!(out).ok();
        return;
    }

    let cat = &layers.catalog;
    let g_count = galaxy_count(layers);
    writeln!(out, "  catalog.root: {}", cat.root.value()).ok();
    writeln!(out, "  galaxies: {g_count}").ok();
    writeln!(
        out,
        "  primary_ref: galaxy={} region={} index={}",
        cat.primary_ref.galaxy_index,
        region_summary(cat.primary_ref.region),
        cat.primary_ref.index_in_region
    )
    .ok();
    writeln!(
        out,
        "  markers_near_focus (catalog): {}",
        cat.markers_near_focus.len()
    )
    .ok();

    for g in 0..g_count {
        let gal = galaxy(layers, g);
        writeln!(
            out,
            "  galaxy[{}]: {} — {:?}, R={:.0} ly",
            gal.index,
            gal.label,
            gal.profile.morphology,
            gal.profile.radius_ly
        )
        .ok();
        let region_count = regions_near(layers, g).len();
        writeln!(out, "    regions intersecting load: {region_count}").ok();
        let markers = markers_near_focus(layers, g);
        writeln!(out, "    markers loaded: {}", markers.len()).ok();
        if verbosity >= 1 {
            for m in &markers {
                writeln!(
                    out,
                    "      marker #{} @ ({:.2}, {:.2}, {:.2}) ly  {:?}/{:?}{}",
                    m.index_in_region,
                    m.position_ly.x_ly,
                    m.position_ly.y_ly,
                    m.position_ly.z_ly,
                    m.phenotype.multiplicity,
                    m.phenotype.primary_class,
                    if m.is_primary { " [PRIMARY]" } else { "" }
                )
                .ok();
            }
        }
    }
    writeln!(out).ok();
}

fn format_layer_solar_system(
    out: &mut String,
    layers: &CosmologyLayers,
    to_layer: PipelineLayerId,
    verbosity: u8,
) {
    layer_heading(
        out,
        4,
        PipelineLayerId::SolarSystem,
        "stage4_solar_system",
        "SolarSystem",
    );
    if to_layer < PipelineLayerId::SolarSystem {
        writeln!(
            out,
            "  (not generated — pipeline stopped before solar_system)"
        )
        .ok();
        writeln!(out).ok();
        return;
    }

    let Some(sys) = layers.solar_system.as_ref() else {
        writeln!(out, "  (missing — stage 4 did not run)").ok();
        writeln!(out).ok();
        return;
    };

    write_solar_system_body(out, sys, verbosity);
    writeln!(out).ok();
}

fn write_solar_system_body(out: &mut String, sys: &SolarSystem, verbosity: u8) {
    writeln!(out, "  system_seed: {}", sys.system_seed.value()).ok();
    writeln!(
        out,
        "  host: galaxy={} region={} index={}",
        sys.galaxy_index,
        region_summary(sys.region),
        sys.index_in_region
    )
    .ok();
    writeln!(
        out,
        "  barycenter_ly: ({:.4}, {:.4}, {:.4})",
        sys.barycenter_ly.x_ly,
        sys.barycenter_ly.y_ly,
        sys.barycenter_ly.z_ly
    )
    .ok();
    writeln!(
        out,
        "  phenotype: {:?} / {:?}",
        sys.phenotype.multiplicity, sys.phenotype.primary_class
    )
    .ok();
    writeln!(out, "  stars: {}", sys.stars.len()).ok();
    for (i, star) in sys.stars.iter().enumerate() {
        writeln!(
            out,
            "    star[{i}]: {:?}  {:.3} L☉  AU ({:.2}, {:.2}, {:.2})",
            star.spectral_class,
            star.luminosity_solar,
            star.position_au.x_au,
            star.position_au.y_au,
            star.position_au.z_au
        )
        .ok();
    }
    writeln!(out, "  planet_slots: {}", sys.planet_slots.len()).ok();
    let show_slots = verbosity >= 1 || sys.planet_slots.len() <= 12;
    if show_slots {
        for slot in &sys.planet_slots {
            writeln!(
                out,
                "    slot[{}]: {:?}  a={:.3} AU  period={:.1} d  e={:.3}",
                slot.index,
                slot.mass_class,
                slot.orbit.semi_major_axis_au,
                slot.orbit.period_days,
                slot.orbit.eccentricity
            )
            .ok();
        }
    } else {
        writeln!(
            out,
            "    ({} slots — use -v to list)",
            sys.planet_slots.len()
        )
        .ok();
    }
    writeln!(out, "  belts: {}", sys.belts.len()).ok();
    if verbosity >= 1 {
        for b in &sys.belts {
            writeln!(
                out,
                "    belt[{}]: {:.2}–{:.2} AU",
                b.index, b.inner_au, b.outer_au
            )
            .ok();
        }
    }
}

fn format_layer_planetary_body(
    out: &mut String,
    layers: &CosmologyLayers,
    to_layer: PipelineLayerId,
    verbosity: u8,
) {
    layer_heading(
        out,
        5,
        PipelineLayerId::PlanetaryBody,
        "stage5_planetary_body",
        "PlanetaryBodies",
    );
    if to_layer < PipelineLayerId::PlanetaryBody {
        writeln!(
            out,
            "  (not generated — pipeline stopped before planetary_body)"
        )
        .ok();
        writeln!(out).ok();
        return;
    }

    let Some(planets) = layers.planetary.as_ref() else {
        writeln!(out, "  (missing — stage 5 did not run)").ok();
        writeln!(out).ok();
        return;
    };

    write_planetary_bodies(out, planets, verbosity);
    writeln!(out).ok();
}

fn write_planetary_bodies(out: &mut String, planets: &PlanetaryBodies, verbosity: u8) {
    writeln!(out, "  focus_slot_index: {}", planets.focus_slot_index).ok();
    writeln!(out, "  bodies: {}", planets.bodies.len()).ok();
    if let Some(focus) = planets.focus() {
        writeln!(
            out,
            "  focus body: slot {} ({:?})  R={:.0} km  map {}×{}",
            focus.slot_index,
            focus.mass_class,
            focus.params.radius_km,
            focus.map_width,
            focus.map_height
        )
        .ok();
        writeln!(
            out,
            "    insolation: {:.0} W/m²  year: {:.1} d  day: {:.1} h  g: {:.2} m/s²",
            focus.params.insolation_wm2(),
            focus.params.orbital_period_days,
            focus.params.rotation_period_hours,
            focus.params.gravity_ms2
        )
        .ok();
    }
    let show_all = verbosity >= 1 || planets.bodies.len() <= 12;
    if show_all {
        for body in &planets.bodies {
            let tag = if body.slot_index == planets.focus_slot_index {
                " [FOCUS]"
            } else {
                ""
            };
            writeln!(
                out,
                "    body[{}]: {:?}{}  R={:.0} km  a={:.3} AU  insolation={:.0} W/m²",
                body.slot_index,
                body.mass_class,
                tag,
                body.params.radius_km,
                body.params.semi_major_axis_au,
                body.params.insolation_wm2()
            )
            .ok();
        }
    } else {
        writeln!(
            out,
            "    ({} bodies — use -v to list all)",
            planets.bodies.len()
        )
        .ok();
    }
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
