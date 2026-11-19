//! Real solar-system bodies as ready-made contexts: starting points for hand-written worlds
//! and references for realism checks. Values are present-day means.

use super::context::{BodyContext, Irradiance, TidalPartner};
use super::interior::{BodyInterior, Shell, ShellKind, TectonicRegime};
use super::matter::{Composition, Phase, Species};
use super::system::{BodyClass, BodyExterior, Elements, Epoch};
use super::AU_KM;
use crate::key::Key;
use std::fmt;
use std::str::FromStr;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Archetype {
    Earth,
    Mars,
    Venus,
    Mercury,
    Moon,
    Titan,
    Europa,
    GasGiant,
    IceGiant,
}

impl Archetype {
    pub const ALL: [Self; 9] =
        [Self::Earth, Self::Mars, Self::Venus, Self::Mercury, Self::Moon, Self::Titan, Self::Europa, Self::GasGiant, Self::IceGiant];

    #[must_use]
    pub fn slug(self) -> &'static str {
        match self {
            Self::Earth => "earth",
            Self::Mars => "mars",
            Self::Venus => "venus",
            Self::Mercury => "mercury",
            Self::Moon => "moon",
            Self::Titan => "titan",
            Self::Europa => "europa",
            Self::GasGiant => "gas_giant",
            Self::IceGiant => "ice_giant",
        }
    }

    /// The body as a full context named `name`.
    #[must_use]
    pub fn context(self, key: Key, name: &str) -> BodyContext {
        let (exterior, orbit, irradiance, tides) = self.outside();
        BodyContext { key, name: name.into(), exterior, interior: self.interior(), orbit, irradiance, tides, epoch: Epoch { system_age_gyr: 4.57 } }
    }

    /// The shell stack and regime.
    #[must_use]
    pub fn interior(self) -> BodyInterior {
        use Phase::{Liquid, Solid, Supercritical};
        use ShellKind::{Core, Crust, IceShell, Mantle, Ocean};
        use Species::*;
        let w = Composition::from_weights;
        let one = Composition::single;
        let (shells, regime) = match self {
            Self::Earth => (
                vec![
                    Shell::new(Core, Solid, 3480.0, w(&[(Iron, 0.85), (Nickel, 0.05), (Sulfur, 0.10)])),
                    Shell::new(Mantle, Solid, 2856.0, one(Silicate)),
                    Shell::new(Crust, Solid, 35.0, w(&[(Basalt, 0.6), (Granite, 0.4)])),
                    Shell::new(Ocean, Liquid, 2.62, one(Water)),
                    Shell::atmosphere(1.013, 100.0, Composition::from_moles(&[(Nitrogen, 0.781), (Oxygen, 0.209), (Argon, 0.0093), (CarbonDioxide, 0.00042)])),
                ],
                TectonicRegime::Plates,
            ),
            Self::Mars => (
                vec![
                    Shell::new(Core, Liquid, 1830.0, w(&[(Iron, 0.8), (Sulfur, 0.2)])),
                    Shell::new(Mantle, Solid, 1509.5, one(Silicate)),
                    Shell::new(Crust, Solid, 50.0, one(Basalt)),
                    Shell::atmosphere(0.006, 11.0, Composition::from_moles(&[(CarbonDioxide, 0.953), (Nitrogen, 0.027), (Argon, 0.016)])),
                ],
                TectonicRegime::StagnantLid,
            ),
            Self::Venus => (
                vec![
                    Shell::new(Core, Solid, 3000.0, w(&[(Iron, 0.85), (Nickel, 0.05), (Sulfur, 0.10)])),
                    Shell::new(Mantle, Solid, 3021.8, one(Silicate)),
                    Shell::new(Crust, Solid, 30.0, one(Basalt)),
                    Shell::atmosphere(92.0, 250.0, Composition::from_moles(&[(CarbonDioxide, 0.965), (Nitrogen, 0.035)])),
                ],
                TectonicRegime::Episodic,
            ),
            Self::Mercury => (
                vec![
                    Shell::new(Core, Liquid, 2020.0, w(&[(Iron, 0.9), (Nickel, 0.05), (Sulfur, 0.05)])),
                    Shell::new(Mantle, Solid, 384.7, one(Silicate)),
                    Shell::new(Crust, Solid, 35.0, one(Basalt)),
                ],
                TectonicRegime::StagnantLid,
            ),
            Self::Moon => (
                vec![
                    Shell::new(Core, Solid, 350.0, one(Iron)),
                    Shell::new(Mantle, Solid, 1347.4, one(Silicate)),
                    Shell::new(Crust, Solid, 40.0, w(&[(Anorthosite, 0.8), (Basalt, 0.2)])),
                ],
                TectonicRegime::StagnantLid,
            ),
            Self::Titan => (
                vec![
                    Shell::new(Core, Solid, 2050.0, one(Silicate)),
                    Shell::new(Mantle, Solid, 300.0, one(Ice)),
                    Shell::new(Ocean, Liquid, 100.0, w(&[(Water, 0.9), (Ammonia, 0.1)])),
                    Shell::new(IceShell, Solid, 124.73, one(Ice)),
                    Shell::atmosphere(1.467, 600.0, Composition::from_moles(&[(Nitrogen, 0.951), (Methane, 0.049)])),
                ],
                TectonicRegime::IceShell,
            ),
            Self::Europa => (
                vec![
                    Shell::new(Core, Solid, 600.0, w(&[(Iron, 0.7), (Sulfur, 0.3)])),
                    Shell::new(Mantle, Solid, 840.8, one(Silicate)),
                    Shell::new(Ocean, Liquid, 100.0, one(Water)),
                    Shell::new(IceShell, Solid, 20.0, one(Ice)),
                ],
                TectonicRegime::IceShell,
            ),
            Self::GasGiant => (
                vec![
                    Shell::new(Core, Solid, 10_000.0, w(&[(Silicate, 0.5), (Iron, 0.5)])),
                    Shell::new(Mantle, Liquid, 40_000.0, w(&[(Hydrogen, 0.75), (Helium, 0.25)])),
                    Shell::atmosphere(1.0e6, 19_911.0, Composition::from_moles(&[(Hydrogen, 0.898), (Helium, 0.102)])),
                ],
                TectonicRegime::Inactive,
            ),
            Self::IceGiant => (
                vec![
                    Shell::new(Core, Solid, 5000.0, w(&[(Silicate, 0.6), (Iron, 0.4)])),
                    Shell::new(Mantle, Supercritical, 12_000.0, w(&[(Water, 0.65), (Ammonia, 0.15), (Methane, 0.2)])),
                    Shell::atmosphere(1.0e5, 7622.0, Composition::from_moles(&[(Hydrogen, 0.8), (Helium, 0.19), (Methane, 0.015)])),
                ],
                TectonicRegime::Inactive,
            ),
        };
        BodyInterior { shells, regime }
    }

    /// Exterior, own orbit, irradiance and tide raisers.
    fn outside(self) -> (BodyExterior, Elements, Irradiance, Vec<TidalPartner>) {
        let ext = |class, mass_kg, radius_km, rotation_h, obliquity_deg| BodyExterior { class, mass_kg, radius_km, rotation_h, obliquity_deg };
        let orbit = |a_au, e, i_deg, lan_deg, argp_deg| Elements { a_au, e, i_deg, lan_deg, argp_deg, m0_deg: 0.0 };
        let sun = |flux_w_m2, eccentricity, obliquity_deg, year_days| Irradiance { flux_w_m2, star_temperature_k: 5772.0, eccentricity, obliquity_deg, year_days };
        let tide = |name: &str, mass_kg, distance_km| vec![TidalPartner { name: name.into(), mass_kg, distance_km }];
        use BodyClass::{GasGiant, IceGiant, Terrestrial};
        match self {
            Self::Earth => (
                ext(Terrestrial, 5.972e24, 6371.0, 23.934, 23.44),
                orbit(1.0, 0.0167, 0.0, 0.0, 102.9),
                sun(1361.0, 0.0167, 23.44, 365.256),
                tide("Moon", 7.342e22, 384_400.0),
            ),
            Self::Mars => (
                ext(Terrestrial, 6.417e23, 3389.5, 24.623, 25.19),
                orbit(1.5237, 0.0934, 1.85, 49.56, 286.5),
                sun(586.2, 0.0934, 25.19, 686.98),
                Vec::new(),
            ),
            // Retrograde spin: obliquity past 90°.
            Self::Venus => (
                ext(Terrestrial, 4.867e24, 6051.8, 5832.5, 177.36),
                orbit(0.7233, 0.0068, 3.39, 76.68, 54.88),
                sun(2601.0, 0.0068, 177.36, 224.70),
                Vec::new(),
            ),
            Self::Mercury => (
                ext(Terrestrial, 3.301e23, 2439.7, 1407.5, 0.034),
                orbit(0.3871, 0.2056, 7.0, 48.33, 29.12),
                sun(9083.0, 0.2056, 0.034, 87.969),
                Vec::new(),
            ),
            Self::Moon => (
                ext(Terrestrial, 7.342e22, 1737.4, 655.72, 6.68),
                orbit(384_400.0 / AU_KM, 0.0549, 5.145, 125.08, 318.15),
                sun(1361.0, 0.0167, 1.54, 365.256),
                tide("Earth", 5.972e24, 384_400.0),
            ),
            Self::Titan => (
                ext(Terrestrial, 1.3452e23, 2574.73, 382.69, 0.3),
                orbit(1_221_870.0 / AU_KM, 0.0288, 0.348, 0.0, 0.0),
                sun(14.9, 0.0565, 26.73, 10_759.0),
                tide("Saturn", 5.683e26, 1_221_870.0),
            ),
            Self::Europa => (
                ext(Terrestrial, 4.7998e22, 1560.8, 85.23, 0.1),
                orbit(671_034.0 / AU_KM, 0.009, 0.47, 0.0, 0.0),
                sun(50.3, 0.0489, 3.1, 4332.6),
                tide("Jupiter", 1.898e27, 671_034.0),
            ),
            Self::GasGiant => (
                ext(GasGiant, 1.898e27, 69_911.0, 9.925, 3.13),
                orbit(5.2034, 0.0489, 1.30, 100.46, 273.87),
                sun(50.3, 0.0489, 3.13, 4332.6),
                Vec::new(),
            ),
            Self::IceGiant => (
                ext(IceGiant, 1.024e26, 24_622.0, 16.11, 28.32),
                orbit(30.07, 0.0086, 1.77, 131.78, 273.19),
                sun(1.505, 0.0086, 28.32, 60_190.0),
                Vec::new(),
            ),
        }
    }
}

impl fmt::Display for Archetype {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.slug())
    }
}

impl FromStr for Archetype {
    type Err = crate::Error;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let want = s.trim().to_ascii_lowercase().replace('-', "_");
        Self::ALL.into_iter().find(|a| a.slug() == want).ok_or_else(|| {
            let names: Vec<_> = Self::ALL.iter().map(|a| a.slug()).collect();
            crate::Error::Invalid(format!("unknown archetype '{s}' (one of {})", names.join(", ")))
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shells_fill_the_radius() {
        for a in Archetype::ALL {
            let ctx = a.context(Key(1), a.slug());
            // Giants count their envelope; solid bodies stop at the surface.
            let counted = |s: &&Shell| s.kind != ShellKind::Atmosphere || !ctx.interior.has_solid_surface();
            let total: f64 = ctx.interior.shells.iter().filter(counted).map(|s| s.thickness_km).sum();
            let r = ctx.exterior.radius_km;
            assert!((total - r).abs() / r < 0.01, "{a}: shells {total} km vs radius {r} km");
        }
    }

    #[test]
    fn real_values() {
        let g = |a: Archetype| a.context(Key(1), "x").gravity_m_s2();
        assert!((g(Archetype::Earth) - 9.81).abs() < 0.05);
        assert!((g(Archetype::Mars) - 3.73).abs() < 0.05);
        assert!((g(Archetype::Moon) - 1.62).abs() < 0.05);
        assert!((g(Archetype::Titan) - 1.35).abs() < 0.05);
        let regime = |a: Archetype| a.interior().regime;
        assert_eq!(regime(Archetype::Earth), TectonicRegime::Plates);
        assert_eq!(regime(Archetype::Venus), TectonicRegime::Episodic);
        assert_eq!(regime(Archetype::Mars), TectonicRegime::StagnantLid);
        assert_eq!(regime(Archetype::Europa), TectonicRegime::IceShell);
    }

    #[test]
    fn names_parse() {
        for a in Archetype::ALL {
            assert_eq!(a.slug().parse::<Archetype>().unwrap(), a);
        }
        assert_eq!("Gas-Giant".parse::<Archetype>().unwrap(), Archetype::GasGiant);
        assert!("pluto".parse::<Archetype>().is_err());
    }

    #[test]
    fn serde_round_trip() {
        for a in Archetype::ALL {
            let ctx = a.context(Key(7), a.slug());
            let back: BodyContext = serde_json::from_str(&serde_json::to_string(&ctx).unwrap()).unwrap();
            assert_eq!(back, ctx);
        }
    }
}
