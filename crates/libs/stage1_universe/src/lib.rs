//! Pipeline stage 1 — universe from root seed.

#![deny(clippy::todo, clippy::unimplemented)]

use geoforge_types::cosmology::{
    CosmologyConfig, CosmologyPreset, OBSERVABLE_HUBBLE_KM_S_MPC, OBSERVABLE_OMEGA_MATTER,
    OBSERVABLE_UNIVERSE_LABEL, Universe,
};
use geoforge_types::Seed;

/// Generate the observable universe layer from a preset.
#[must_use]
pub fn generate(root: Seed, preset: CosmologyPreset) -> Universe {
    generate_with_config(root, CosmologyConfig::from_preset(preset))
}

/// Generate the observable universe layer with explicit config.
///
/// Cosmological parameters (H₀, Ω_m) are fixed; the seed drives later stages only.
#[must_use]
pub fn generate_with_config(root: Seed, config: CosmologyConfig) -> Universe {
    Universe {
        root,
        label: OBSERVABLE_UNIVERSE_LABEL.to_string(),
        hubble_km_s_mpc: OBSERVABLE_HUBBLE_KM_S_MPC,
        omega_matter: OBSERVABLE_OMEGA_MATTER,
        config,
    }
}
