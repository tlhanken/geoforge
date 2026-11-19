//! Gray-atmosphere radiative balance: how much of a body's thermal emission its atmosphere
//! traps, from composition and pressure. Pure functions of the body, like [`crate::derived`].
//!
//! Outgoing longwave is `σT⁴ / (1 + ¾τ)` (Eddington gray atmosphere). The optical depth τ
//! sums three parts:
//!
//! - **CO₂ bands**, growing with the CO₂ partial pressure and broadened by total pressure;
//! - **collision-induced absorption** of dense gas, ∝ P²/g, strongest for cold bodies that
//!   emit in the far infrared (Titan's N₂–CH₄);
//! - **water vapour** where liquid water stands on the surface, rising with temperature (the
//!   water-vapour feedback). Its slope sets Earth's climate sensitivity.
//!
//! The constants are fitted to the presets: Earth (OLR 231 W/m² at 15 °C, slope 2.09 W/m²/K),
//! Mars (≈ −60 °C), Venus (≈ 460 °C) and Titan (≈ −179 °C).

use crate::body::{BodyContext, BodyInterior, Species};
use crate::derived::STEFAN_BOLTZMANN;

pub const KELVIN: f64 = 273.15;

const CO2_K: f64 = 4.71;
const CO2_PARTIAL_EXP: f64 = 0.4445;
const CO2_PRESSURE_EXP: f64 = 0.2857;
const CIA_K: f64 = 0.44;
const WATER_TAU_15C: f64 = 0.7635;
const WATER_PER_C: f64 = 0.0143;

/// A body's atmosphere as a gray absorber.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Greenhouse {
    /// CO₂ band optical depth (independent of temperature).
    pub tau_co2: f64,
    /// Collision-induced optical depth at 100 K; scales as (100 K / T)².
    pub tau_cia_100k: f64,
    /// Water-vapour optical depth at 15 °C; 0 without surface water.
    pub tau_water_15c: f64,
}

impl Greenhouse {
    #[must_use]
    pub fn of(body: &BodyContext) -> Self {
        let i = &body.interior;
        let p = i.surface_pressure_bar();
        let co2 = i.partial_pressure_bar(Species::CarbonDioxide);
        Self {
            tau_co2: if co2 > 0.0 { CO2_K * co2.powf(CO2_PARTIAL_EXP) * p.powf(CO2_PRESSURE_EXP) } else { 0.0 },
            tau_cia_100k: CIA_K * p * p / body.gravity_m_s2(),
            tau_water_15c: if i.has_surface_liquid_water() { WATER_TAU_15C } else { 0.0 },
        }
    }
    #[must_use]
    pub fn optical_depth(&self, t_c: f64) -> f64 {
        let tk = (t_c + KELVIN).max(20.0);
        self.tau_co2 + self.tau_cia_100k * (100.0 / tk).powi(2) + self.tau_water_15c * (WATER_PER_C * (t_c - 15.0).min(300.0)).exp()
    }
    /// Outgoing longwave at surface temperature `t_c`, W/m².
    #[must_use]
    pub fn olr_w_m2(&self, t_c: f64) -> f64 {
        let tk = (t_c + KELVIN).max(1.0);
        STEFAN_BOLTZMANN * tk.powi(4) / (1.0 + 0.75 * self.optical_depth(t_c))
    }
    /// d(OLR)/dT, W/m²/K.
    #[must_use]
    pub fn olr_slope(&self, t_c: f64) -> f64 {
        self.olr_w_m2(t_c + 0.5) - self.olr_w_m2(t_c - 0.5)
    }
}

/// Bond albedo where the surface is free of ice: clouds on wet worlds, the bright cloud decks
/// of very thick atmospheres, and bare regolith on airless bodies.
#[must_use]
pub fn open_albedo(interior: &BodyInterior) -> f64 {
    let p = interior.surface_pressure_bar();
    if !interior.has_atmosphere() {
        return 0.12;
    }
    let clouds = if interior.has_surface_liquid_water() { 0.05 } else { 0.0 };
    0.25 + clouds + 0.52 * p * p / (p * p + 400.0)
}

/// Global mean surface temperature from radiative balance alone (no ice, no transport), °C.
/// The coldest balance: a body past the water-vapour runaway limit reports the hot end.
#[must_use]
pub fn balance_temperature_c(body: &BodyContext) -> f64 {
    let absorbed = body.irradiance.flux_w_m2 / 4.0 * (1.0 - open_albedo(&body.interior));
    let g = Greenhouse::of(body);
    let mut lo = -270.0;
    while lo < 2000.0 {
        let hi = lo + 2.0;
        if g.olr_w_m2(hi) >= absorbed {
            let (a, b) = (g.olr_w_m2(lo), g.olr_w_m2(hi));
            return lo + 2.0 * ((absorbed - a) / (b - a)).clamp(0.0, 1.0);
        }
        lo = hi;
    }
    2000.0
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::body::Archetype;
    use crate::key::Key;

    fn body(a: Archetype) -> BodyContext {
        a.context(Key(1), a.slug())
    }

    #[test]
    fn earth_calibration() {
        let g = Greenhouse::of(&body(Archetype::Earth));
        assert!((g.olr_w_m2(15.0) - 231.4).abs() < 2.0, "OLR {}", g.olr_w_m2(15.0));
        assert!((g.olr_slope(15.0) - 2.09).abs() < 0.1, "slope {}", g.olr_slope(15.0));
    }

    #[test]
    fn other_worlds() {
        let t = |a| balance_temperature_c(&body(a));
        assert!((t(Archetype::Venus) - 464.0).abs() < 30.0, "Venus {}", t(Archetype::Venus));
        assert!((t(Archetype::Mars) + 60.0).abs() < 10.0, "Mars {}", t(Archetype::Mars));
        assert!((t(Archetype::Titan) + 179.0).abs() < 6.0, "Titan {}", t(Archetype::Titan));
        // An airless body has no greenhouse: its balance is the bare equilibrium temperature.
        let moon = body(Archetype::Moon);
        let teq = crate::derived::equilibrium_temperature_k(moon.irradiance.flux_w_m2, 0.12) - KELVIN;
        assert!((t(Archetype::Moon) - teq).abs() < 0.5);
    }

    #[test]
    fn more_co2_is_warmer() {
        let mut b = body(Archetype::Mars);
        let t0 = balance_temperature_c(&b);
        let air = b.interior.shells.iter_mut().find(|s| s.kind == crate::body::ShellKind::Atmosphere).unwrap();
        air.pressure_bar = Some(1.0);
        assert!(balance_temperature_c(&b) > t0 + 20.0);
    }
}
