//! Run pre-tectonic stages and query helpers.

use geoforge_types::coordinates::GalacticPoint;
use geoforge_types::cosmology::galaxy_seed;
use geoforge_types::cosmology::{
    regions_intersecting_load, CosmologyConfig, CosmologyLayers, CosmologyPreset, RegionId,
    SolarSystem, StellarSystemMarker, StellarSystemRef,
};
use geoforge_types::{PipelineLayerId, Seed};
use stage3_galaxy::markers_in_region as stage3_markers_in_region;

/// Run stages 1–5 up to `to_layer` (stops before tectonics).
#[must_use]
pub fn run_pre_tectonic_layers(
    root: Seed,
    preset: CosmologyPreset,
    to_layer: PipelineLayerId,
) -> CosmologyLayers {
    run_pre_tectonic_layers_with_config(root, CosmologyConfig::from_preset(preset), to_layer)
}

/// Run with explicit config.
#[must_use]
pub fn run_pre_tectonic_layers_with_config(
    root: Seed,
    config: CosmologyConfig,
    to_layer: PipelineLayerId,
) -> CosmologyLayers {
    let universe = stage1_universe::generate_with_config(root, config);
    let supercluster = stage2_supercluster::generate(&universe);
    let catalog = stage3_galaxy::generate(&supercluster);

    let mut layers = CosmologyLayers {
        universe,
        supercluster,
        catalog,
        solar_system: None,
        planetary: None,
    };
    finish_zoom_stages(&mut layers, to_layer);
    layers
}

/// Rebuild stage 3+ after focus change.
#[must_use]
pub fn with_focus(
    mut layers: CosmologyLayers,
    focus: GalacticPoint,
    to_layer: PipelineLayerId,
) -> CosmologyLayers {
    layers.supercluster.focus_ly = focus;
    layers.catalog = stage3_galaxy::generate(&layers.supercluster);
    layers.solar_system = None;
    layers.planetary = None;
    finish_zoom_stages(&mut layers, to_layer);
    layers
}

fn finish_zoom_stages(layers: &mut CosmologyLayers, to_layer: PipelineLayerId) {
    if to_layer >= PipelineLayerId::SolarSystem {
        let reference = layers.catalog.primary_ref;
        layers.solar_system = Some(stage4_solar_system::generate(&layers.catalog, reference));
    } else {
        layers.solar_system = None;
        layers.planetary = None;
        return;
    }

    if to_layer >= PipelineLayerId::PlanetaryBody {
        if let Some(ref solar) = layers.solar_system {
            layers.planetary = Some(stage5_planetary_body::generate(solar));
        }
    } else {
        layers.planetary = None;
    }
}

/// Number of galaxies.
#[must_use]
pub fn galaxy_count(layers: &CosmologyLayers) -> u32 {
    layers.catalog.galaxies.len() as u32
}

/// Galaxy by index.
#[must_use]
pub fn galaxy(layers: &CosmologyLayers, index: u32) -> geoforge_types::cosmology::GalaxyInstance {
    layers
        .catalog
        .galaxy(index)
        .cloned()
        .unwrap_or_else(|| stage3_galaxy::galaxy_instance(layers.root(), index, layers.config()))
}

/// Load sphere for region queries.
#[must_use]
pub fn load_sphere(layers: &CosmologyLayers) -> geoforge_types::coordinates::LoadSphere {
    geoforge_types::coordinates::LoadSphere {
        center: layers.focus_ly(),
        radius_ly: layers.config().load_radius_ly,
    }
}

/// Regions near focus.
#[must_use]
pub fn regions_near(layers: &CosmologyLayers, galaxy_index: u32) -> Vec<RegionId> {
    let g = galaxy(layers, galaxy_index);
    regions_intersecting_load(&g.profile, load_sphere(layers))
}

/// Primary system reference.
#[must_use]
pub fn primary_ref(layers: &CosmologyLayers) -> StellarSystemRef {
    layers.catalog.primary_ref
}

/// Markers in one region.
#[must_use]
pub fn markers_in_region(
    layers: &CosmologyLayers,
    galaxy_index: u32,
    region: RegionId,
) -> Vec<StellarSystemMarker> {
    let g = galaxy(layers, galaxy_index);
    let gseed = galaxy_seed(layers.root(), galaxy_index);
    let primary = primary_ref(layers);
    let primary_key = Some((primary.region, primary.index_in_region));
    stage3_markers_in_region(
        galaxy_index,
        gseed,
        &g.profile,
        region,
        layers.config(),
        primary_key,
    )
}

/// Markers near focus.
#[must_use]
pub fn markers_near_focus(layers: &CosmologyLayers, galaxy_index: u32) -> Vec<StellarSystemMarker> {
    if galaxy_index == 0 {
        return layers.catalog.markers_near_focus.clone();
    }
    let mut all = Vec::new();
    for region in regions_near(layers, galaxy_index) {
        for m in markers_in_region(layers, galaxy_index, region) {
            all.push(m);
            if all.len() as u32 >= layers.config().max_markers_per_query {
                return all;
            }
        }
    }
    all
}

/// Primary solar system (stage 4 must have been generated).
#[must_use]
pub fn primary_solar_system(layers: &CosmologyLayers) -> &SolarSystem {
    layers
        .solar_system
        .as_ref()
        .expect("solar_system layer not generated")
}

/// Regenerate solar system at a reference.
#[must_use]
pub fn solar_system_at(layers: &CosmologyLayers, reference: StellarSystemRef) -> SolarSystem {
    stage4_solar_system::generate(&layers.catalog, reference)
}

/// Primary marker.
#[must_use]
pub fn primary_marker(layers: &CosmologyLayers) -> StellarSystemMarker {
    layers
        .catalog
        .primary_marker()
        .cloned()
        .unwrap_or_else(|| {
            let pref = primary_ref(layers);
            markers_in_region(layers, pref.galaxy_index, pref.region)
                .into_iter()
                .find(|m| m.is_primary)
                .expect("minimal preset guarantees one marker")
        })
}

/// Alias for older name.
#[must_use]
#[deprecated(note = "use run_pre_tectonic_layers")]
pub fn run_cosmology_layers(
    root: Seed,
    preset: CosmologyPreset,
    to_layer: PipelineLayerId,
) -> CosmologyLayers {
    run_pre_tectonic_layers(root, preset, to_layer)
}
