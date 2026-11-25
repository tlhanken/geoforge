//! What sits at a galaxy's centre and round it: its black hole, how active its nucleus is,
//! and its globular clusters — drawn from the galaxy's model, its halo and its seed.
//!
//! - **The black hole** follows the classical bulge's mass where there is one (ellipticals,
//!   lenticulars, early spirals: M = 0.49×10⁹ M☉ (M_bulge/10¹¹ M☉)^1.16, 0.29 dex scatter;
//!   Kormendy & Ho 2013) and otherwise the galaxy's stellar mass (log M = 7.45 + 1.05
//!   log(M★/10¹¹ M☉), 0.55 dex; Reines & Volonteri 2015), the pseudobulges of later spirals
//!   and dwarfs lying under the bulge relation (the Milky Way's 4.3×10⁶ M☉ among them). Not
//!   every dwarf has one: the share with a central black hole falls from all at 10¹⁰ M☉ to
//!   six in ten at 10⁹ M☉ and a fifth at 10⁸ (Greene et al. 2020; M33 has none over 1,500
//!   M☉).
//! - **Its activity**: three in ten nuclei are feeding, at an Eddington ratio with more at
//!   the low end (p(λ) ∝ λ^−0.6 per decade, Aird et al. 2012) from 10⁻⁴ to 0.3 in galaxies
//!   still forming stars and 10⁻⁶ to 10⁻² in old ones (Kauffmann & Heckman 2009); the rest
//!   are quiet, like Sgr A* (10⁻⁹ to 10⁻⁵). The light, λ × 1.26×10³⁸ erg/s per M☉, names
//!   it: a quasar from 10⁴⁵ erg/s, a Seyfert nucleus from 10⁴³, a LINER from 10⁴¹ (Ho 2008
//!   finds a tenth of nearby bright galaxies Seyferts, a fifth to a third LINERs).
//!   Three in ten active ellipticals and lenticulars launch radio jets, one in twelve others.
//! - **Globular clusters** weigh 2.9×10⁻⁵ of the halo (Harris et al. 2017), each of a
//!   log-normal mass (10⁵ M☉ at its peak, 0.52 dex wide: mean 2×10⁵ M☉), so the Milky Way's
//!   halo holds about 180 (it has some 160). Six in ten are metal-poor ([Fe/H] −1.5, 12.5 Gyr
//!   old) and spread through the halo, four in ten metal-rich (−0.55, 11 Gyr) and nearer the
//!   bulge; dwarfs' are nearly all metal-poor. Their number falls as r^−3.5 beyond a core
//!   (3 kpc and 1 kpc for the Milky Way's, scaling as the halo's radius), out to a third of
//!   the halo's radius (the metal-poor) or a twentieth (the metal-rich), so half the Milky
//!   Way's lie within ~5 kpc of its centre (Harris 1996); half their light within ~3 pc;
//!   M/L_V 1.9.

use crate::draw::Draw;
use crate::galaxy::{ComponentKind, Galaxy, Morphology, SpiralStage};
use geoforge_core::seed::Seed;

/// The Milky Way's halo, M☉.
pub const MILKY_WAY_HALO_MSUN: f64 = 1.3e12;
/// Sgr A*'s mass, M☉ (GRAVITY Collaboration 2019).
pub const SGR_A_MSUN: f64 = 4.3e6;
/// Eddington luminosity per M☉, erg/s.
const EDDINGTON_ERG_S: f64 = 1.26e38;

/// How bright a nucleus is.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Activity {
    Quiet,
    Liner,
    Seyfert,
    Quasar,
}

impl Activity {
    #[must_use]
    pub fn label(self) -> &'static str {
        match self {
            Self::Quiet => "quiet",
            Self::Liner => "low-luminosity active nucleus (LINER)",
            Self::Seyfert => "Seyfert nucleus",
            Self::Quasar => "quasar",
        }
    }
}

/// A galaxy's centre.
#[derive(Debug, Clone, PartialEq)]
pub struct Nucleus {
    /// Its black hole's mass, M☉; `None` for a dwarf without one.
    pub black_hole_msun: Option<f64>,
    /// Eddington ratio.
    pub eddington_ratio: f64,
    /// Light from its accretion, erg/s.
    pub luminosity_erg_s: f64,
    pub activity: Activity,
    pub jets: bool,
}

/// A globular cluster.
#[derive(Debug, Clone, PartialEq)]
pub struct Globular {
    pub index: u32,
    pub name: String,
    /// pc from the galaxy's centre.
    pub position_pc: [f64; 3],
    pub mass_msun: f64,
    pub metallicity_dex: f64,
    pub age_gyr: f64,
    /// Half its light within this, pc.
    pub half_light_pc: f64,
    pub absolute_v: f64,
}

impl Globular {
    #[must_use]
    pub fn distance_pc(&self, p: [f64; 3]) -> f64 {
        (0..3).map(|k| (self.position_pc[k] - p[k]).powi(2)).sum::<f64>().sqrt()
    }
}

/// Whether `galaxy`'s bulge is classical (its black hole follows its mass).
fn classical_bulge(galaxy: &Galaxy) -> bool {
    matches!(galaxy.morphology, Morphology::Elliptical | Morphology::Lenticular | Morphology::Spiral { stage: SpiralStage::A, .. })
}

/// The centre of `galaxy`, drawn from `seed`.
#[must_use]
pub fn nucleus(galaxy: &Galaxy, seed: Seed) -> Nucleus {
    let mut d = Draw::new(seed.child(b"nucleus"));
    let stars = galaxy.mass_msun().max(1.0);
    let bulge: f64 = galaxy.components.iter().filter(|c| matches!(c.kind, ComponentKind::Bulge | ComponentKind::Spheroid)).map(|c| c.mass_msun).sum();
    let occupied = (0.6 + 0.4 * (stars / 1e9).log10()).clamp(0.1, 1.0);
    let has = d.u() < occupied;
    let log_m =
        if classical_bulge(galaxy) && bulge > 0.0 { (0.49e9f64).log10() + 1.16 * (bulge / 1e11).log10() + d.normal(0.0, 0.29) } else { 7.45 + 1.05 * (stars / 1e11).log10() + d.normal(0.0, 0.55) };
    let black_hole_msun = has.then(|| 10f64.powf(log_m).max(1e3));
    activity(black_hole_msun, galaxy.morphology, &mut d)
}

/// The Milky Way's centre: Sgr A*, quiet (an Eddington ratio of ~10⁻⁹).
#[must_use]
pub fn milky_way() -> Nucleus {
    let eddington_ratio = 1e-9;
    Nucleus { black_hole_msun: Some(SGR_A_MSUN), eddington_ratio, luminosity_erg_s: eddington_ratio * EDDINGTON_ERG_S * SGR_A_MSUN, activity: Activity::Quiet, jets: false }
}

fn activity(black_hole_msun: Option<f64>, morphology: Morphology, d: &mut Draw) -> Nucleus {
    let Some(m) = black_hole_msun else { return Nucleus { black_hole_msun: None, eddington_ratio: 0.0, luminosity_erg_s: 0.0, activity: Activity::Quiet, jets: false } };
    let early = matches!(morphology, Morphology::Elliptical | Morphology::Lenticular);
    let feeding = d.u() < 0.3;
    // log λ with density ∝ 10^(−0.6 x) per decade if feeding, from −4 to −0.5 in a galaxy
    // still forming stars and −6 to −2 in an old one; −9 to −5 if not feeding.
    let (lo, hi) = if early { (-6.0, -2.0) } else { (-4.0, -0.5) };
    let u = d.u();
    let log_ratio = if feeding {
        let (a, b) = (10f64.powf(-0.6 * lo), 10f64.powf(-0.6 * hi));
        -(a + u * (b - a)).log10() / 0.6
    } else {
        -9.0 + 4.0 * u
    };
    let eddington_ratio = 10f64.powf(log_ratio);
    let luminosity_erg_s = eddington_ratio * EDDINGTON_ERG_S * m;
    let activity = match luminosity_erg_s {
        l if l >= 1e45 => Activity::Quasar,
        l if l >= 1e43 => Activity::Seyfert,
        l if l >= 1e41 => Activity::Liner,
        _ => Activity::Quiet,
    };
    let jets = activity != Activity::Quiet && d.u() < if early { 0.3 } else { 0.08 };
    Nucleus { black_hole_msun, eddington_ratio, luminosity_erg_s, activity, jets }
}

/// The globular clusters of `galaxy` in a halo of `halo_msun`, drawn from `seed`.
#[must_use]
pub fn globulars(galaxy: &Galaxy, halo_msun: f64, seed: Seed) -> Vec<Globular> {
    const PEAK_LOG_M: f64 = 5.0;
    const WIDTH_DEX: f64 = 0.52;
    let mean_mass = 10f64.powf(PEAK_LOG_M + WIDTH_DEX * WIDTH_DEX * std::f64::consts::LN_10 / 2.0);
    let mut d = Draw::new(seed.child(b"globulars"));
    let n = d.poisson(2.9e-5 * halo_msun / mean_mass);
    // Sizes scale with the halo's radius (∝ M^⅓); the Milky Way's virial radius is ~280 kpc.
    let scale = (halo_msun / MILKY_WAY_HALO_MSUN).cbrt();
    let virial_pc = 280_000.0 * scale;
    let dwarf = matches!(galaxy.morphology, Morphology::DwarfSpheroidal | Morphology::DwarfIrregular) || galaxy.mass_msun() < 1e9;
    let rich_share = if dwarf { 0.1 } else { 0.4 };
    (0..n as u32)
        .map(|index| {
            let mut g = Draw::new(seed.child(b"globular").child_u64(b"n", u64::from(index)));
            let rich = g.u() < rich_share;
            let (core_pc, reach_pc, feh, age_gyr) =
                if rich { (1000.0 * scale, virial_pc / 20.0, g.normal(-0.55, 0.25), g.normal(11.0, 0.8)) } else { (3000.0 * scale, virial_pc / 3.0, g.normal(-1.5, 0.3), g.normal(12.5, 0.6)) };
            // N(<r) ∝ 1 − (1 + r/core)^−½, falling as r^−3.5 far out: r = core·((1−u)^−2 − 1),
            // up to the reach.
            let top = 1.0 - (1.0 + reach_pc / core_pc).powf(-0.5);
            let r = core_pc * ((1.0 - g.u() * top).powf(-2.0) - 1.0);
            let position_pc = g.direction().map(|v| v * r.max(50.0));
            let mass_msun = 10f64.powf(g.normal(PEAK_LOG_M, WIDTH_DEX)).clamp(1e3, 5e6);
            let half_light_pc = 3.0 * 10f64.powf(g.normal(0.0, 0.25));
            let absolute_v = 4.83 - 2.5 * (mass_msun / 1.9).log10();
            let name = format!("{} Cluster", geoforge_names::name(seed.child(b"globular name").child_u64(b"n", u64::from(index)).0));
            Globular { index, name, position_pc, mass_msun, metallicity_dex: feh.min(0.2), age_gyr: age_gyr.clamp(9.0, 13.5), half_light_pc, absolute_v }
        })
        .collect()
}

/// A globular cluster as seen from a place.
#[derive(Debug, Clone, PartialEq)]
pub struct SeenGlobular {
    pub cluster: Globular,
    /// Unit vector toward it in the system's frame.
    pub direction: [f64; 3],
    pub distance_pc: f64,
    /// Its integrated magnitude through the dust.
    pub v: f64,
    /// The radius holding most of its light (three half-light radii), degrees.
    pub radius_deg: f64,
}

/// Of `clusters`, those `observer` sees brighter than `limit` through `field`'s dust.
#[must_use]
pub fn seen(clusters: &[Globular], field: &crate::starfield::Starfield, observer: &crate::sky::Observer, limit: f64) -> Vec<SeenGlobular> {
    let here = observer.position_pc;
    // Clouds within a few kpc; the smooth layer the whole way.
    let dust = crate::sectors::dust::DustView::new(field.sectors(), here, 2500.0);
    let mut out: Vec<SeenGlobular> = clusters
        .iter()
        .filter_map(|c| {
            let d = c.distance_pc(here).max(1.0);
            let v = c.absolute_v + 5.0 * (d / 10.0).log10() + dust.to(c.position_pc);
            let g = [0, 1, 2].map(|k| (c.position_pc[k] - here[k]) / d);
            (v <= limit).then(|| SeenGlobular { cluster: c.clone(), direction: observer.to_system_frame(g), distance_pc: d, v, radius_deg: (3.0 * c.half_light_pc / d).atan().to_degrees() })
        })
        .collect();
    out.sort_by(|a, b| a.v.total_cmp(&b.v));
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::galaxy::GalaxyPreset;

    /// Black holes by bulge and by stars: a big elliptical's billions, a Milky Way-like
    /// spiral's millions (Sgr A* inside the scatter), many dwarfs without; activity mostly
    /// quiet, a few per cent Seyferts among big galaxies, quasars rarer.
    #[test]
    fn black_holes_and_activity() {
        let (mut big, mut mw) = (Vec::new(), Vec::new());
        let elliptical = GalaxyPreset::Elliptical.galaxy();
        let spiral = Galaxy::milky_way();
        for s in 0..2000u64 {
            big.push(nucleus(&elliptical, Seed(s)));
            mw.push(nucleus(&spiral, Seed(s)));
        }
        let median = |v: &[Nucleus]| {
            let mut m: Vec<f64> = v.iter().filter_map(|n| n.black_hole_msun).collect();
            m.sort_by(f64::total_cmp);
            m[m.len() / 2]
        };
        assert!(median(&big) > 1e8 && (3e6..3e7).contains(&median(&mw)), "{:e} {:e}", median(&big), median(&mw));
        let share = |v: &[Nucleus], a: Activity| v.iter().filter(|n| n.activity >= a).count() as f64 / v.len() as f64;
        assert!(share(&mw, Activity::Seyfert) > 0.005 && share(&mw, Activity::Seyfert) < 0.1, "{}", share(&mw, Activity::Seyfert));
        assert!(share(&big, Activity::Quasar) < 0.005 && share(&big, Activity::Liner) > 0.15 && share(&big, Activity::Seyfert) < 0.1);
        let dwarf = GalaxyPreset::DwarfSpheroidal.galaxy();
        let without = (0..400u64).filter(|&s| nucleus(&dwarf, Seed(s)).black_hole_msun.is_none()).count();
        assert!(without > 100, "{without} of 400 dwarfs without");
        assert_eq!(nucleus(&spiral, Seed(9)), nucleus(&spiral, Seed(9)));
    }

    /// The Milky Way's halo holds about 180 globulars (it has ~160), three-fifths metal-poor,
    /// half within ~5 kpc and most within 30 kpc, the metal-rich nearer; a dwarf's halo a
    /// handful or none.
    #[test]
    fn globular_systems() {
        let mw = globulars(&Galaxy::milky_way(), MILKY_WAY_HALO_MSUN, Seed(3));
        assert!((120..260).contains(&mw.len()), "{}", mw.len());
        let poor = mw.iter().filter(|g| g.metallicity_dex < -1.0).count() as f64 / mw.len() as f64;
        assert!((0.45..0.75).contains(&poor), "{poor}");
        let r = |g: &Globular| g.distance_pc([0.0; 3]);
        let mut by: Vec<f64> = mw.iter().map(r).collect();
        by.sort_by(f64::total_cmp);
        assert!((2500.0..8000.0).contains(&by[by.len() / 2]) && by[by.len() * 4 / 5] < 30_000.0, "median {}", by[by.len() / 2]);
        let mean = |rich: bool| {
            let v: Vec<f64> = mw.iter().filter(|g| (g.metallicity_dex > -1.0) == rich).map(r).collect();
            v.iter().sum::<f64>() / v.len() as f64
        };
        assert!(mean(true) < mean(false));
        assert!(mw.iter().all(|g| (-11.0..-1.5).contains(&g.absolute_v) && g.half_light_pc > 0.5 && g.half_light_pc < 15.0));
        assert!(globulars(&GalaxyPreset::DwarfSpheroidal.galaxy(), 1e9, Seed(3)).len() < 5);
    }
}
