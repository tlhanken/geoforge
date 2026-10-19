//! Stage 1 output — universe-scale constants and generation config.

use serde::{Deserialize, Serialize};

use super::preset::CosmologyConfig;
use crate::Seed;

/// Hubble constant for the observable universe (Planck 2018, km/s/Mpc).
pub const OBSERVABLE_HUBBLE_KM_S_MPC: f64 = 67.4;

/// Matter density parameter Ω_m for the observable universe (Planck 2018).
pub const OBSERVABLE_OMEGA_MATTER: f64 = 0.315;

/// Default display label for stage 1 (not seed-randomized).
pub const OBSERVABLE_UNIVERSE_LABEL: &str = "Observable universe (ΛCDM)";

/// Universe-scale layer produced from the root seed (pipeline stage 1).
///
/// Physical cosmological parameters are fixed to observed ΛCDM values; only
/// [`CosmologyConfig`] (preset / load caps) and [`root`](Self::root) vary per run.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Universe {
    /// Root procedural seed (stable identity for the run).
    pub root: Seed,
    /// Human-readable label (fixed for the observable universe).
    pub label: String,
    /// Hubble parameter (km/s/Mpc).
    pub hubble_km_s_mpc: f64,
    /// Matter density parameter Ω_m (dimensionless).
    pub omega_matter: f64,
    /// Cosmology generation config (galaxy count, load caps, preset).
    pub config: CosmologyConfig,
}

impl Universe {
    /// Observed-universe cosmological parameters (no seed randomization).
    #[must_use]
    pub fn observable_constants() -> (f64, f64) {
        (OBSERVABLE_HUBBLE_KM_S_MPC, OBSERVABLE_OMEGA_MATTER)
    }
}
