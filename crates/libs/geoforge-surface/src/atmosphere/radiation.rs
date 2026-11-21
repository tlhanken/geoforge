//! The energy balance by latitude lives in the core crate, so the solar-system layer's CO₂
//! thermostat can use it too.

pub use geoforge_core::energy_balance::{seasonal_swing, Radiation};
