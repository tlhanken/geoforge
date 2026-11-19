//! The body's own layer: the shell stack, core outward, and how its outer lid behaves.

use super::matter::{Composition, Phase, Species};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ShellKind {
    Core,
    Mantle,
    Crust,
    Ocean,
    IceShell,
    Atmosphere,
}

/// One layer of a body, innermost first in [`BodyInterior::shells`].
///
/// For an ocean, `thickness_km` is the global-equivalent layer: the depth if the water were
/// spread evenly over the whole body. For an atmosphere, `pressure_bar` is surface pressure.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Shell {
    pub kind: ShellKind,
    pub phase: Phase,
    pub thickness_km: f64,
    pub composition: Composition,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pressure_bar: Option<f64>,
}

impl Shell {
    #[must_use]
    pub fn new(kind: ShellKind, phase: Phase, thickness_km: f64, composition: Composition) -> Self {
        Self { kind, phase, thickness_km, composition, pressure_bar: None }
    }
    #[must_use]
    pub fn atmosphere(pressure_bar: f64, thickness_km: f64, composition: Composition) -> Self {
        Self { kind: ShellKind::Atmosphere, phase: Phase::Gas, thickness_km, composition, pressure_bar: Some(pressure_bar) }
    }
}

/// How the outer solid layer moves and renews itself. Picks the surface generators.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum TectonicRegime {
    /// A mobile lid broken into moving plates that recycle crust (Earth).
    Plates,
    /// One rigid lid; heat escapes by conduction and hotspot volcanism (Mars, Mercury, Moon).
    StagnantLid,
    /// A lid that founders and resurfaces the planet in bursts (Venus).
    Episodic,
    /// A floating ice crust over a liquid ocean (Europa, Enceladus, Titan).
    IceShell,
    /// No geologically active solid surface (giants, small bodies).
    Inactive,
}

/// The shell stack, core outward, and the regime of its lid.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BodyInterior {
    pub shells: Vec<Shell>,
    pub regime: TectonicRegime,
}

impl BodyInterior {
    #[must_use]
    pub fn shell(&self, kind: ShellKind) -> Option<&Shell> {
        self.shells.iter().rev().find(|s| s.kind == kind)
    }
    /// The outermost solid layer people could stand on (rock crust or ice shell), if any.
    #[must_use]
    pub fn surface_shell(&self) -> Option<&Shell> {
        self.surface_index().map(|i| &self.shells[i])
    }
    fn surface_index(&self) -> Option<usize> {
        self.shells.iter().rposition(|s| matches!(s.kind, ShellKind::Crust | ShellKind::IceShell))
    }
    #[must_use]
    pub fn has_solid_surface(&self) -> bool {
        // A solid shell must not lie under a deep gas envelope (gas giants).
        let Some(idx) = self.surface_index() else { return false };
        self.shells[idx + 1..].iter().all(|s| match s.kind {
            ShellKind::Atmosphere => s.pressure_bar.unwrap_or(0.0) < 1000.0,
            ShellKind::Ocean => true,
            _ => false,
        })
    }
    /// The surface is rock crust (not ice).
    #[must_use]
    pub fn has_lithosphere(&self) -> bool {
        self.has_solid_surface() && self.surface_shell().is_some_and(|s| s.kind == ShellKind::Crust)
    }
    #[must_use]
    pub fn has_atmosphere(&self) -> bool {
        self.surface_pressure_bar() > 1e-3
    }
    #[must_use]
    pub fn surface_pressure_bar(&self) -> f64 {
        self.shell(ShellKind::Atmosphere).and_then(|s| s.pressure_bar).unwrap_or(0.0)
    }
    /// Partial pressure of `s` at the surface, bar.
    #[must_use]
    pub fn partial_pressure_bar(&self, s: Species) -> f64 {
        self.shell(ShellKind::Atmosphere).map_or(0.0, |a| a.pressure_bar.unwrap_or(0.0) * a.composition.mole_fraction(s))
    }
    /// Global-equivalent depth of oceans lying on the solid surface, km.
    #[must_use]
    pub fn surface_ocean_gel_km(&self) -> f64 {
        let Some(surface) = self.surface_index() else { return 0.0 };
        self.shells[surface + 1..].iter().filter(|s| s.kind == ShellKind::Ocean).map(|s| s.thickness_km).sum()
    }
    /// Liquid water stands on the surface, so rain, rivers and the water cycle apply.
    #[must_use]
    pub fn has_surface_liquid_water(&self) -> bool {
        let Some(surface) = self.surface_index() else { return false };
        self.shells[surface + 1..].iter().any(|s| {
            s.kind == ShellKind::Ocean && s.phase == Phase::Liquid && s.thickness_km > 0.0 && s.composition.fraction(Species::Water) > 0.5
        })
    }
}

#[cfg(test)]
mod tests {
    use super::super::Archetype;
    use super::*;

    #[test]
    fn presence_from_shells() {
        let earth = Archetype::Earth.interior();
        assert!(earth.has_solid_surface() && earth.has_lithosphere() && earth.has_atmosphere());
        assert!(earth.has_surface_liquid_water());
        assert!((earth.surface_ocean_gel_km() - 2.62).abs() < 1e-9);

        let giant = Archetype::GasGiant.interior();
        assert!(!giant.has_solid_surface() && !giant.has_lithosphere());

        // Ice shell over an ocean over rock: a solid surface, but no lithosphere and no
        // surface water.
        let europa = Archetype::Europa.interior();
        assert!(europa.has_solid_surface() && !europa.has_lithosphere());
        assert_eq!(europa.surface_ocean_gel_km(), 0.0, "subsurface ocean is not a surface ocean");
        assert!(!europa.has_surface_liquid_water());

        let mars = Archetype::Mars.interior();
        assert!(mars.has_lithosphere() && mars.has_atmosphere() && !mars.has_surface_liquid_water());
        assert!(!Archetype::Moon.interior().has_atmosphere());
    }

    #[test]
    fn partial_pressures() {
        let venus = Archetype::Venus.interior();
        let co2 = venus.partial_pressure_bar(Species::CarbonDioxide);
        assert!((co2 - 0.965 * 92.0).abs() < 0.5, "Venus pCO₂ {co2}");
        let earth = Archetype::Earth.interior();
        let co2 = earth.partial_pressure_bar(Species::CarbonDioxide);
        assert!((co2 - 4.2e-4).abs() < 1e-4, "Earth pCO₂ {co2}");
        assert_eq!(Archetype::Moon.interior().partial_pressure_bar(Species::Nitrogen), 0.0);
    }
}
