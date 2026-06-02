//! [`CosmologyContext`] — on-demand cosmology from root seed.

use geoforge_types::coordinates::{GalacticPoint, LoadSphere};
use geoforge_types::cosmology::{
    regions_intersecting_load, CosmologyConfig, CosmologyPreset, GalaxyInstance, RegionId,
    SolarSystem, StellarSystemMarker, StellarSystemRef,
};
use geoforge_types::Seed;

use crate::galaxy::{galaxy_count, galaxy_instance};
use crate::marker::markers_in_region;
use crate::seeds::galaxy_seed;
use crate::solar_system::generate_solar_system;

/// Seed-driven cosmology generator (no disk persistence).
#[derive(Debug, Clone)]
pub struct CosmologyContext {
    /// Root seed.
    pub root: Seed,
    /// Generation parameters.
    pub config: CosmologyConfig,
    /// Focus for local region loading (ly).
    pub focus_ly: GalacticPoint,
}

impl CosmologyContext {
    /// Create with preset defaults.
    #[must_use]
    pub fn new(root: Seed, preset: CosmologyPreset) -> Self {
        Self {
            root,
            config: CosmologyConfig::from_preset(preset),
            focus_ly: GalacticPoint::origin(),
        }
    }

    /// Create with full config.
    #[must_use]
    pub fn with_config(root: Seed, config: CosmologyConfig) -> Self {
        Self {
            root,
            config,
            focus_ly: GalacticPoint::origin(),
        }
    }

    /// Set load focus (e.g. primary system position).
    #[must_use]
    pub fn with_focus(mut self, focus: GalacticPoint) -> Self {
        self.focus_ly = focus;
        self
    }

    /// Number of galaxies in this run.
    #[must_use]
    pub fn galaxy_count(&self) -> u32 {
        galaxy_count(&self.config)
    }

    /// Galaxy metadata by index.
    #[must_use]
    pub fn galaxy(&self, index: u32) -> GalaxyInstance {
        galaxy_instance(self.root, index, &self.config)
    }

    /// Load sphere for region queries.
    #[must_use]
    pub fn load_sphere(&self) -> LoadSphere {
        LoadSphere {
            center: self.focus_ly,
            radius_ly: self.config.load_radius_ly,
        }
    }

    /// Region ids intersecting the current load sphere for a galaxy.
    #[must_use]
    pub fn regions_near(&self, galaxy_index: u32) -> Vec<RegionId> {
        let galaxy = self.galaxy(galaxy_index);
        regions_intersecting_load(&galaxy.profile, self.load_sphere())
    }

    /// Primary system reference (deterministic).
    #[must_use]
    pub fn primary_ref(&self) -> StellarSystemRef {
        let galaxy = self.galaxy(0);
        let region = match galaxy.profile.morphology {
            geoforge_types::cosmology::GalaxyMorphology::Elliptical
            | geoforge_types::cosmology::GalaxyMorphology::Bubble => RegionId::Spherical {
                shell: galaxy.profile.ring_count / 2,
                theta_wedge: 0,
                phi_wedge: 0,
            },
            _ => RegionId::Cylindrical {
                ring: galaxy.profile.ring_count / 3,
                wedge: 0,
                layer: 0_u16,
            },
        };
        StellarSystemRef {
            galaxy_index: 0,
            region,
            index_in_region: 0,
        }
    }

    /// All markers in one region.
    #[must_use]
    pub fn markers_in_region(
        &self,
        galaxy_index: u32,
        region: RegionId,
    ) -> Vec<StellarSystemMarker> {
        let galaxy = self.galaxy(galaxy_index);
        let gseed = galaxy_seed(self.root, galaxy_index);
        let primary = self.primary_ref();
        let primary_key = Some((primary.region, primary.index_in_region));
        markers_in_region(
            self.root,
            galaxy_index,
            gseed,
            &galaxy.profile,
            region,
            &self.config,
            primary_key,
        )
    }

    /// Markers for all regions near focus, capped by config.
    pub fn markers_near_focus(&self, galaxy_index: u32) -> Vec<StellarSystemMarker> {
        let mut all = Vec::new();
        for region in self.regions_near(galaxy_index) {
            for m in self.markers_in_region(galaxy_index, region) {
                all.push(m);
                if all.len() as u32 >= self.config.max_markers_per_query {
                    return all;
                }
            }
        }
        all
    }

    /// Regenerate full solar system at reference.
    #[must_use]
    pub fn solar_system(&self, reference: StellarSystemRef) -> SolarSystem {
        let galaxy = self.galaxy(reference.galaxy_index);
        let gseed = galaxy_seed(self.root, reference.galaxy_index);
        let sys = crate::seeds::system_seed(
            crate::seeds::region_seed(gseed, reference.region),
            reference.index_in_region,
        );
        let force_g = self.config.preset == CosmologyPreset::Minimal;
        let phenotype = crate::phenotype::derive_phenotype(sys, force_g);
        generate_solar_system(
            self.root,
            reference.galaxy_index,
            &galaxy.profile,
            reference,
            phenotype,
        )
    }

    /// Primary marker (first generated in minimal, or flagged primary).
    #[must_use]
    pub fn primary_marker(&self) -> StellarSystemMarker {
        let pref = self.primary_ref();
        let markers = self.markers_in_region(pref.galaxy_index, pref.region);
        markers
            .into_iter()
            .find(|m| m.is_primary)
            .unwrap_or_else(|| {
                self.markers_in_region(0, pref.region)
                    .into_iter()
                    .next()
                    .expect("minimal preset guarantees one marker")
            })
    }
}
