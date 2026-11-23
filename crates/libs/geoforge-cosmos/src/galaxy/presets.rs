//! Galaxies by morphology.
//!
//! The Milky Way's numbers are Bland-Hawthorn & Gerhard (2016)'s: thin disk 300 ± 50 pc high
//! and 2.6 ± 0.5 kpc long ((3.5 ± 1)×10¹⁰ M☉); thick disk 900 ± 180 pc and 2.0 ± 0.2 kpc,
//! (6 ± 3)×10⁹ M☉, 4 ± 2% of the thin disk's density at the Sun; halo a power law of slope
//! −2.5 ± 0.3 to a break at 25 ± 10 kpc and −3.7 to −5 beyond, flattened to 0.65, (4–7)×10⁸
//! M☉; a box/peanut bulge and bar of (1.4–1.7)×10¹⁰ M☉ with exponential scale lengths
//! 0.70 : 0.44 : 0.18 kpc, its long axis 27° ± 2° from the Sun–centre line; a nuclear cluster
//! of (1.8 ± 0.3)×10⁷ M☉ with a half-light radius of 4.2 pc; the Sun 8.2 ± 0.1 kpc from the
//! centre and 25 ± 5 pc above the plane. Four logarithmic arms of pitch 13.1° ± 0.6°
//! (Vallée 2015) from 2.2 kpc, turning clockwise seen from the north, the Sun between two of
//! them. The thin disk's stars thicken from 50 pc at birth to ~350 pc (a mass-weighted 270
//! pc), as its older stars are hotter (Bovy et al. 2012).
//!
//! Other morphologies follow the Hubble sequence in their proportions (bulges falling and
//! arms opening from Sa to Sd, star formation over in S0 and E) at masses typical of each;
//! they are illustrations, not fits.

use super::{Arms, Component, ComponentKind, Dust, Galaxy, Height, Imf, Metallicity, Shape, StarFormation};
use serde::{Deserialize, Serialize};

/// A galaxy's place on the Hubble sequence, and the dwarfs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Morphology {
    Spiral { stage: SpiralStage, barred: bool },
    Lenticular,
    Elliptical,
    DwarfSpheroidal,
    DwarfIrregular,
    Starburst,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SpiralStage {
    A,
    B,
    Bc,
    C,
    D,
}

impl Morphology {
    /// Its Hubble type: SBbc, Sa, S0, E, dSph, dIrr.
    #[must_use]
    pub fn label(self) -> String {
        match self {
            Self::Spiral { stage, barred } => {
                let s = match stage {
                    SpiralStage::A => "a",
                    SpiralStage::B => "b",
                    SpiralStage::Bc => "bc",
                    SpiralStage::C => "c",
                    SpiralStage::D => "d",
                };
                format!("{}{s}", if barred { "SB" } else { "S" })
            }
            Self::Lenticular => "S0".into(),
            Self::Elliptical => "E".into(),
            Self::DwarfSpheroidal => "dSph".into(),
            Self::DwarfIrregular => "dIrr".into(),
            Self::Starburst => "starburst".into(),
        }
    }
}

/// Galaxies to choose from by name.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GalaxyPreset {
    /// The Milky Way (SBbc), as the module notes give it.
    MilkyWay,
    /// A spiral of a stage, with or without a bar.
    Spiral(SpiralStage, bool),
    /// A disk and bulge whose star formation is over.
    Lenticular,
    /// A giant elliptical: stars formed in the first 2 Gyr, slightly metal-rich, none forming
    /// now. No blue stars; many white dwarfs and red giants.
    Elliptical,
    /// A faint old metal-poor spheroid, like the Milky Way's satellites.
    DwarfSpheroidal,
    /// A dwarf irregular like the Small Magellanic Cloud: a fifth of solar metals, so few
    /// giant planets.
    Dwarf,
    /// A disk galaxy in a starburst: a quarter of its stars born in the last 100 Myr, so
    /// blue stars and fresh remnants are common.
    Starburst,
}

impl GalaxyPreset {
    pub const ALL: [Self; 14] = [
        Self::MilkyWay,
        Self::Spiral(SpiralStage::A, false),
        Self::Spiral(SpiralStage::B, false),
        Self::Spiral(SpiralStage::C, false),
        Self::Spiral(SpiralStage::D, false),
        Self::Spiral(SpiralStage::A, true),
        Self::Spiral(SpiralStage::B, true),
        Self::Spiral(SpiralStage::C, true),
        Self::Spiral(SpiralStage::D, true),
        Self::Lenticular,
        Self::Elliptical,
        Self::DwarfSpheroidal,
        Self::Dwarf,
        Self::Starburst,
    ];

    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Self::MilkyWay => "milky-way",
            Self::Spiral(SpiralStage::A, false) => "sa",
            Self::Spiral(SpiralStage::B, false) => "sb",
            Self::Spiral(SpiralStage::Bc, false) => "sbc-unbarred",
            Self::Spiral(SpiralStage::C, false) => "sc",
            Self::Spiral(SpiralStage::D, false) => "sd",
            Self::Spiral(SpiralStage::A, true) => "sba",
            Self::Spiral(SpiralStage::B, true) => "sbb",
            Self::Spiral(SpiralStage::Bc, true) => "sbbc",
            Self::Spiral(SpiralStage::C, true) => "sbc",
            Self::Spiral(SpiralStage::D, true) => "sbd",
            Self::Lenticular => "lenticular",
            Self::Elliptical => "elliptical",
            Self::DwarfSpheroidal => "dwarf-spheroidal",
            Self::Dwarf => "dwarf",
            Self::Starburst => "starburst",
        }
    }

    /// The preset named `name` (as [`Self::name`] gives it).
    #[must_use]
    pub fn from_name(name: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|p| p.name().eq_ignore_ascii_case(name))
    }

    #[must_use]
    pub fn galaxy(self) -> Galaxy {
        match self {
            Self::MilkyWay => Galaxy::milky_way(),
            Self::Spiral(stage, barred) => spiral(stage, barred),
            Self::Lenticular => lenticular(),
            Self::Elliptical => elliptical(),
            Self::DwarfSpheroidal => dwarf_spheroidal(),
            Self::Dwarf => dwarf_irregular(),
            Self::Starburst => starburst(),
        }
    }
}

/// The galaxy of the preset named in a system's site, the Milky Way if none matches.
#[must_use]
pub fn by_name(name: &str) -> Galaxy {
    GalaxyPreset::ALL.into_iter().map(GalaxyPreset::galaxy).find(|g| g.name == name).unwrap_or_default()
}

fn halo(mass_msun: f64) -> Component {
    Component {
        kind: ComponentKind::Halo,
        mass_msun,
        shape: Shape::PowerLaw { inner: -2.5, outer: -3.7, break_kpc: 25.0, core_kpc: 1.0, q: 0.65, edge_kpc: 100.0 },
        history: StarFormation::between(13.0, 10.0),
        metals: Metallicity::flat(-1.5, 0.5),
        arm_contrast: 0.0,
    }
}

fn thick_disk(mass_msun: f64, scale_length_kpc: f64, dex: f64) -> Component {
    Component {
        kind: ComponentKind::ThickDisk,
        mass_msun,
        shape: Shape::Disk { scale_length_kpc, height: Height::Fixed { pc: 900.0 } },
        history: StarFormation::between(12.0, 8.0),
        metals: Metallicity::flat(dex - 0.5, 0.25),
        arm_contrast: 0.0,
    }
}

fn nucleus(mass_msun: f64, dex: f64) -> Component {
    Component { kind: ComponentKind::Nucleus, mass_msun, shape: Shape::Plummer { a_pc: 3.2 }, history: StarFormation::between(12.0, 0.0), metals: Metallicity::flat(dex + 0.2, 0.3), arm_contrast: 0.0 }
}

fn thin_disk(mass_msun: f64, scale_length_kpc: f64, decline_gyr: Option<f64>, metals: Metallicity, arm_contrast: f64) -> Component {
    Component {
        kind: ComponentKind::ThinDisk,
        mass_msun,
        shape: Shape::Disk { scale_length_kpc, height: Height::Aging { young_pc: 50.0, old_pc: 350.0, tau_gyr: 3.0 } },
        history: StarFormation { decline_gyr, ..StarFormation::between(11.0, 0.0) },
        metals,
        arm_contrast,
    }
}

impl Galaxy {
    /// The Milky Way, as the module notes give it, but for its thin disk's mass: 2.15×10¹⁰ M☉
    /// sets 0.08 systems per cubic parsec at the Sun (the 10 pc census: 339 systems, Reylé
    /// et al. 2021), where Bland-Hawthorn & Gerhard's 3.5 × 10¹⁰ gives 0.13 with this disk's
    /// thin young layers. With the galactic IMF it gives the Sun's sky: about 9,000 stars to
    /// magnitude 6.5 (9,110 in the Yale Bright Star Catalogue), 20–25 brighter than 1.5.
    #[must_use]
    pub fn milky_way() -> Self {
        let dex = 0.0;
        Self {
            name: "Milky Way".into(),
            morphology: Morphology::Spiral { stage: SpiralStage::Bc, barred: true },
            imf: Imf::galactic(),
            components: vec![
                nucleus(1.8e7, dex),
                Component {
                    kind: ComponentKind::Bulge,
                    mass_msun: 1.55e10,
                    shape: Shape::Ellipsoid { scales_kpc: [0.70, 0.44, 0.18], angle_deg: -27.0 },
                    history: StarFormation::between(12.0, 7.0),
                    metals: Metallicity::flat(-0.1, 0.4),
                    arm_contrast: 0.0,
                },
                thin_disk(2.15e10, 2.6, None, Metallicity { at_reference_dex: 0.0, gradient_dex_per_kpc: -0.06, per_gyr_dex: -0.02, scatter_dex: 0.2 }, 0.8),
                thick_disk(6.0e9, 2.0, dex),
                halo(5.5e8),
            ],
            arms: Some(Arms { count: 4, pitch_deg: 13.1, start_kpc: 2.2, phase_deg: -9.0, width_kpc: 0.4 }),
            dust: Dust { mag_per_kpc: 0.8, height_pc: 100.0, scale_length_kpc: 3.5, cloud_gas: 2.0 },
            reference_kpc: 8.2,
            edge_kpc: 15.0,
        }
    }
}

/// A spiral of `stage`: total stellar mass, bulge share, disk scale length, star formation
/// decline, arms (count, pitch), arm contrast and metals by stage.
fn spiral(stage: SpiralStage, barred: bool) -> Galaxy {
    let (mass, bulge, rd, decline, (count, pitch), contrast, dex) = match stage {
        SpiralStage::A => (8e10, 0.40, 3.5, Some(3.0), (2, 8.0), 0.5, 0.1),
        SpiralStage::B => (6e10, 0.25, 3.0, Some(6.0), (2, 12.0), 0.7, 0.0),
        SpiralStage::Bc => (5e10, 0.15, 2.7, None, (4, 14.0), 0.8, 0.0),
        SpiralStage::C => (3e10, 0.08, 2.5, None, (4, 18.0), 0.8, -0.1),
        SpiralStage::D => (8e9, 0.02, 2.0, Some(-8.0), (4, 25.0), 0.7, -0.3),
    };
    let s = 0.25 * rd;
    let bulge_shape = if barred { Shape::Ellipsoid { scales_kpc: [3.0 * s, s, 0.4 * s], angle_deg: 25.0 } } else { Shape::Ellipsoid { scales_kpc: [s, s, 0.6 * s], angle_deg: 0.0 } };
    let disk = mass * (1.0 - bulge) * 0.92;
    let morphology = Morphology::Spiral { stage, barred };
    Galaxy {
        name: format!("{} spiral", morphology.label()),
        morphology,
        imf: Imf::galactic(),
        components: vec![
            nucleus(5e-4 * mass, dex),
            Component {
                kind: ComponentKind::Bulge,
                mass_msun: mass * bulge,
                shape: bulge_shape,
                history: StarFormation::between(12.0, 6.0),
                metals: Metallicity::flat(dex - 0.1, 0.4),
                arm_contrast: 0.0,
            },
            thin_disk(disk, rd, decline, Metallicity { at_reference_dex: dex, gradient_dex_per_kpc: -0.06 * 2.6 / rd, per_gyr_dex: -0.02, scatter_dex: 0.2 }, contrast),
            thick_disk(mass * (1.0 - bulge) * 0.08, 0.8 * rd, dex),
            halo(0.01 * mass),
        ],
        arms: Some(Arms { count, pitch_deg: pitch, start_kpc: 0.85 * rd, phase_deg: 0.0, width_kpc: 0.15 * rd }),
        dust: Dust { mag_per_kpc: 0.5, height_pc: 100.0, scale_length_kpc: 1.3 * rd, cloud_gas: 2.0 },
        reference_kpc: 3.15 * rd,
        edge_kpc: 6.0 * rd,
    }
}

fn lenticular() -> Galaxy {
    let mass = 6e10;
    Galaxy {
        name: "Lenticular".into(),
        morphology: Morphology::Lenticular,
        imf: Imf::galactic(),
        components: vec![
            Component {
                kind: ComponentKind::Bulge,
                mass_msun: 0.45 * mass,
                shape: Shape::Ellipsoid { scales_kpc: [1.2, 1.2, 0.7], angle_deg: 0.0 },
                history: StarFormation::between(12.0, 7.0),
                metals: Metallicity::flat(0.05, 0.3),
                arm_contrast: 0.0,
            },
            Component {
                kind: ComponentKind::ThinDisk,
                mass_msun: 0.5 * mass,
                shape: Shape::Disk { scale_length_kpc: 3.0, height: Height::Fixed { pc: 500.0 } },
                history: StarFormation { decline_gyr: Some(2.0), ..StarFormation::between(12.0, 4.0) },
                metals: Metallicity { at_reference_dex: 0.05, gradient_dex_per_kpc: -0.04, per_gyr_dex: -0.01, scatter_dex: 0.2 },
                arm_contrast: 0.0,
            },
            thick_disk(0.04 * mass, 2.5, 0.05),
            halo(0.01 * mass),
        ],
        arms: None,
        dust: Dust { mag_per_kpc: 0.1, ..Dust::none() },
        reference_kpc: 9.0,
        edge_kpc: 18.0,
    }
}

fn elliptical() -> Galaxy {
    Galaxy {
        name: "Elliptical".into(),
        morphology: Morphology::Elliptical,
        imf: Imf::galactic(),
        components: vec![Component {
            kind: ComponentKind::Spheroid,
            mass_msun: 1e11,
            shape: Shape::Hernquist { a_kpc: 2.75, q: 0.7 },
            history: StarFormation { age_gyr: 12.0, youngest_gyr: 0.0, decline_gyr: Some(0.7), burst_share: 0.0, burst_gyr: 0.1 },
            metals: Metallicity { at_reference_dex: 0.15, gradient_dex_per_kpc: -0.02, per_gyr_dex: 0.0, scatter_dex: 0.25 },
            arm_contrast: 0.0,
        }],
        arms: None,
        dust: Dust::none(),
        reference_kpc: 4.0,
        edge_kpc: 30.0,
    }
}

fn dwarf_spheroidal() -> Galaxy {
    Galaxy {
        name: "Dwarf spheroidal".into(),
        morphology: Morphology::DwarfSpheroidal,
        imf: Imf::galactic(),
        components: vec![Component {
            kind: ComponentKind::Spheroid,
            mass_msun: 2e7,
            shape: Shape::Plummer { a_pc: 300.0 },
            history: StarFormation::between(12.0, 8.0),
            metals: Metallicity::flat(-1.6, 0.4),
            arm_contrast: 0.0,
        }],
        arms: None,
        dust: Dust::none(),
        reference_kpc: 0.3,
        edge_kpc: 2.0,
    }
}

fn dwarf_irregular() -> Galaxy {
    Galaxy {
        name: "Dwarf irregular".into(),
        morphology: Morphology::DwarfIrregular,
        imf: Imf::galactic(),
        components: vec![Component {
            kind: ComponentKind::ThinDisk,
            mass_msun: 5e8,
            shape: Shape::Disk { scale_length_kpc: 1.0, height: Height::Fixed { pc: 300.0 } },
            history: StarFormation::between(12.0, 0.0),
            metals: Metallicity { at_reference_dex: -0.7, gradient_dex_per_kpc: 0.0, per_gyr_dex: -0.05, scatter_dex: 0.25 },
            arm_contrast: 0.0,
        }],
        arms: None,
        dust: Dust { mag_per_kpc: 0.3, height_pc: 150.0, scale_length_kpc: 1.2, cloud_gas: 1.0 },
        reference_kpc: 1.0,
        edge_kpc: 4.0,
    }
}

fn starburst() -> Galaxy {
    let mut g = Galaxy::milky_way();
    g.name = "Starburst".into();
    g.morphology = Morphology::Starburst;
    let thin = g.components.iter_mut().find(|c| c.kind == ComponentKind::ThinDisk).expect("a thin disk");
    thin.history = StarFormation { burst_share: 0.25, burst_gyr: 0.1, ..thin.history.clone() };
    thin.metals.at_reference_dex = 0.1;
    g.dust = Dust { mag_per_kpc: 1.5, cloud_gas: 4.0, ..g.dust };
    g
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn presets_by_name() {
        for p in GalaxyPreset::ALL {
            assert_eq!(GalaxyPreset::from_name(p.name()), Some(p));
            let g = p.galaxy();
            assert!(!g.components.is_empty() && g.mass_msun() > 1e6, "{}", p.name());
            assert_eq!(by_name(&g.name).name, g.name);
        }
        assert_eq!(Galaxy::milky_way().morphology.label(), "SBbc");
        // An elliptical forms no stars now; a spiral and a dwarf irregular do.
        let forming = |g: &Galaxy| g.components.iter().any(|c| c.history.youngest_gyr == 0.0 && c.history.lookback_density(0.05) > 0.01);
        assert!(!forming(&GalaxyPreset::Elliptical.galaxy()) && forming(&GalaxyPreset::Spiral(SpiralStage::C, false).galaxy()) && forming(&GalaxyPreset::Dwarf.galaxy()));
    }
}
