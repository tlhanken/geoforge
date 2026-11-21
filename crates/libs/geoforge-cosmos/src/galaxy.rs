//! Galaxy models: the stellar population a system is drawn from.
//!
//! A galaxy here is the handful of numbers that decide what a random star of it is like: how
//! stars are born by mass (the initial mass function), when (the star formation history), how
//! metal-rich (by birth time and place), and where (the radial profile). [`Galaxy::milky_way`]
//! is the default; [`GalaxyPreset`] adds an old elliptical, a metal-poor dwarf and a starburst
//! to show what these numbers do.

use serde::{Deserialize, Serialize};

/// The stellar population of a galaxy.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Galaxy {
    pub name: String,
    pub imf: Imf,
    pub history: StarFormation,
    pub metals: Metallicity,
    pub layout: Layout,
}

/// Initial mass function: a broken power law dN/dm ∝ m^−α, continuous at the breaks.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Imf {
    /// Segment edges, M☉, rising: one more than `slopes`.
    pub breaks: Vec<f64>,
    /// α of each segment.
    pub slopes: Vec<f64>,
}

impl Imf {
    /// Kroupa (2001) over hydrogen-burning stars: α = 1.3 from 0.08 to 0.5 M☉, 2.3 above, to
    /// 150 M☉. Brown dwarfs (α = 0.3 below 0.08) are not drawn.
    #[must_use]
    pub fn kroupa() -> Self {
        Self { breaks: vec![0.08, 0.5, 150.0], slopes: vec![1.3, 2.3] }
    }
}

/// When stars were born.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StarFormation {
    /// Age of the oldest stars, Gyr.
    pub age_gyr: f64,
    /// E-folding time of the star formation rate's decline, Gyr: small means most stars are
    /// old; negative means the rate has been rising. `None`: a constant rate.
    pub decline_gyr: Option<f64>,
    /// Share of all stars born in a recent burst, and the burst's length, Gyr.
    pub burst_share: f64,
    pub burst_gyr: f64,
}

/// How metal-rich a star is born: `[Fe/H] = at_reference + gradient·(R − R_ref) +
/// per_gyr·age`, scattered by `scatter_dex`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Metallicity {
    /// [Fe/H] of stars born now at the reference radius, dex.
    pub at_reference_dex: f64,
    /// Radial gradient, dex per kpc.
    pub gradient_dex_per_kpc: f64,
    /// Change per Gyr of age (negative: older stars are poorer), dex.
    pub per_gyr_dex: f64,
    /// Scatter at a given age and radius, dex.
    pub scatter_dex: f64,
}

/// Where stars are: an exponential profile out to an edge.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Layout {
    pub scale_length_kpc: f64,
    pub edge_kpc: f64,
    /// Radius where metallicity takes its reference value (the Sun's orbit, for the Milky Way).
    pub reference_kpc: f64,
}

impl Galaxy {
    /// The Milky Way's thin disk, where nearly all of its stars near the Sun belong:
    ///
    /// - Kroupa (2001) initial mass function;
    /// - stars forming at a roughly constant rate for 11 Gyr (solar-neighbourhood star
    ///   formation began 11–10.5 Gyr ago, ChronoGal 2024; Binney et al. 2000 date it 11.2);
    /// - young stars near the Sun at solar metallicity, older ones poorer by 0.02 dex per Gyr,
    ///   0.2 dex of scatter, and a radial gradient of −0.06 dex/kpc (open clusters, Cepheids,
    ///   red clump giants);
    /// - an exponential disk of scale length 2.6 kpc out to 15 kpc, the Sun at 8.2 kpc.
    #[must_use]
    pub fn milky_way() -> Self {
        Self {
            name: "Milky Way".into(),
            imf: Imf::kroupa(),
            history: StarFormation { age_gyr: 11.0, decline_gyr: None, burst_share: 0.0, burst_gyr: 0.1 },
            metals: Metallicity { at_reference_dex: 0.0, gradient_dex_per_kpc: -0.06, per_gyr_dex: -0.02, scatter_dex: 0.2 },
            layout: Layout { scale_length_kpc: 2.6, edge_kpc: 15.0, reference_kpc: 8.2 },
        }
    }
}

impl Default for Galaxy {
    fn default() -> Self {
        Self::milky_way()
    }
}

/// Galaxies to choose from by name.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GalaxyPreset {
    /// The Milky Way's disk ([`Galaxy::milky_way`]).
    MilkyWay,
    /// A giant elliptical: stars formed in the first 2 Gyr, slightly metal-rich, none forming
    /// now. No blue stars; many white dwarfs and red giants.
    Elliptical,
    /// A dwarf irregular like the Small Magellanic Cloud: a fifth of solar metals, so few
    /// giant planets.
    Dwarf,
    /// A disk galaxy in a starburst: a quarter of its stars born in the last 100 Myr, so
    /// blue stars and fresh remnants are common.
    Starburst,
}

impl GalaxyPreset {
    pub const ALL: [Self; 4] = [Self::MilkyWay, Self::Elliptical, Self::Dwarf, Self::Starburst];

    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Self::MilkyWay => "milky-way",
            Self::Elliptical => "elliptical",
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
        let mw = Galaxy::milky_way();
        match self {
            Self::MilkyWay => mw,
            Self::Elliptical => Galaxy {
                name: "Elliptical".into(),
                history: StarFormation { age_gyr: 12.0, decline_gyr: Some(0.7), burst_share: 0.0, burst_gyr: 0.1 },
                metals: Metallicity { at_reference_dex: 0.15, gradient_dex_per_kpc: -0.02, per_gyr_dex: 0.0, scatter_dex: 0.25 },
                layout: Layout { scale_length_kpc: 4.0, edge_kpc: 30.0, reference_kpc: 4.0 },
                ..mw
            },
            Self::Dwarf => Galaxy {
                name: "Dwarf irregular".into(),
                history: StarFormation { age_gyr: 12.0, decline_gyr: None, burst_share: 0.0, burst_gyr: 0.1 },
                metals: Metallicity { at_reference_dex: -0.7, gradient_dex_per_kpc: 0.0, per_gyr_dex: -0.05, scatter_dex: 0.25 },
                layout: Layout { scale_length_kpc: 1.0, edge_kpc: 4.0, reference_kpc: 1.0 },
                ..mw
            },
            Self::Starburst => Galaxy {
                name: "Starburst".into(),
                history: StarFormation { age_gyr: 11.0, decline_gyr: None, burst_share: 0.25, burst_gyr: 0.1 },
                metals: Metallicity { at_reference_dex: 0.1, ..mw.metals.clone() },
                ..mw
            },
        }
    }
}
