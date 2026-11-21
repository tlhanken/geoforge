//! A system's phenotype drawn from a galaxy: where it is, when it was born, how metal-rich,
//! and its stars at birth. This is one entry of a galaxy's roster (the solar-systems layer);
//! the orbital tree is built from it.
//!
//! Stars come in singles, binaries and triples, more often the heavier the primary
//! (Duchêne & Kraus 2013): about 25% of red dwarfs have a companion (Winters et al. 2019;
//! Reylé et al. 2021), 46% of Sun-like stars (Raghavan et al. 2010) and over 90% of O stars,
//! and about a quarter of multiples are triples (12% of Sun-like stars, Raghavan et al. 2010).
//! Periods are log-normal: for Sun-like stars a peak at log P = 5.03 days with a spread of
//! 2.28 dex (Raghavan et al. 2010), closer for red dwarfs and much closer for massive stars.
//! Mass ratios are flat; orbits under 12 days are circular, wider ones evenly eccentric.

use crate::draw::Draw;
use crate::galaxy::{Galaxy, Imf, StarFormation};
use crate::stars::evolve;
use geoforge_core::body::{GalacticSite, Star};
use geoforge_core::seed::Seed;

/// One system of a galaxy's roster: its site, birth and stars as they were born.
#[derive(Debug, Clone, PartialEq)]
pub struct Phenotype {
    pub site: GalacticSite,
    pub age_gyr: f64,
    pub metallicity_dex: f64,
    /// Initial mass of the primary, M☉.
    pub primary_mass: f64,
    /// The primary's companion, then (in a triple) an outer star orbiting the inner pair.
    pub companions: Vec<Companion>,
}

/// A companion star at birth.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Companion {
    pub initial_mass: f64,
    /// Semi-major axis at birth, AU: around the primary, or (for the outer star of a triple)
    /// around the inner pair.
    pub a_au: f64,
    pub e: f64,
}

impl Phenotype {
    /// Number of stars.
    #[must_use]
    pub fn multiplicity(&self) -> usize {
        1 + self.companions.len()
    }

    /// The primary as it is now.
    #[must_use]
    pub fn primary(&self) -> Star {
        evolve(self.primary_mass, self.age_gyr, self.metallicity_dex)
    }
}

/// Draw a system from `galaxy`.
#[must_use]
pub fn draw(galaxy: &Galaxy, seed: Seed) -> Phenotype {
    let mut d = Draw::new(seed.child(b"phenotype"));
    let radius_kpc = galactic_radius(&mut d, galaxy);
    let age_gyr = birth_lookback(&mut d, &galaxy.history);
    let m = &galaxy.metals;
    let metallicity_dex = (m.at_reference_dex + m.gradient_dex_per_kpc * (radius_kpc - galaxy.layout.reference_kpc) + m.per_gyr_dex * age_gyr + d.normal(0.0, m.scatter_dex)).clamp(-2.5, 0.6);
    let primary_mass = draw_mass(&mut d, &galaxy.imf);

    let mut companions = Vec::new();
    if d.u() < multiple_fraction(primary_mass) {
        let inner = companion(&mut d, primary_mass, primary_mass, None);
        companions.push(inner);
        if d.u() < TRIPLE_SHARE {
            let inner_log_p = log_period_days(inner.a_au, primary_mass + inner.initial_mass);
            companions.push(companion(&mut d, primary_mass, primary_mass + inner.initial_mass, Some(inner_log_p)));
        }
    }
    Phenotype { site: GalacticSite { galaxy: galaxy.name.clone(), radius_kpc }, age_gyr, metallicity_dex, primary_mass, companions }
}

/// Share of multiples that are triples (Raghavan et al. 2010: 12% of 46%).
const TRIPLE_SHARE: f64 = 0.27;

/// Distance from the centre, kpc: an exponential disk (surface density ∝ e^(−R/h), so radii
/// follow R·e^(−R/h)) cut at its edge.
fn galactic_radius(d: &mut Draw, g: &Galaxy) -> f64 {
    let h = g.layout.scale_length_kpc;
    loop {
        let r = -h * (d.u().max(1e-12).ln() + d.u().max(1e-12).ln());
        if r <= g.layout.edge_kpc {
            return r;
        }
    }
}

/// How long ago the system was born, Gyr.
fn birth_lookback(d: &mut Draw, h: &StarFormation) -> f64 {
    let u = d.u();
    if d.u() < h.burst_share {
        return u * h.burst_gyr;
    }
    match h.decline_gyr {
        // Rate ∝ e^(−t/τ) since the first stars: lookback density ∝ e^(lookback/τ).
        Some(tau) if tau.abs() > 1e-6 => {
            let k = (h.age_gyr / tau).clamp(-700.0, 700.0);
            (tau * (1.0 + u * k.exp_m1()).ln()).clamp(0.0, h.age_gyr)
        }
        _ => u * h.age_gyr,
    }
}

/// An initial mass from the IMF, by inverting its cumulative distribution segment by segment.
fn draw_mass(d: &mut Draw, imf: &Imf) -> f64 {
    // Continuity: k₀ = 1, k_{i+1} = k_i · b_{i+1}^(α_{i+1} − α_i).
    let mut k = 1.0;
    let mut segments = Vec::with_capacity(imf.slopes.len());
    for (i, &alpha) in imf.slopes.iter().enumerate() {
        if i > 0 {
            k *= imf.breaks[i].powf(alpha - imf.slopes[i - 1]);
        }
        let (lo, hi) = (imf.breaks[i], imf.breaks[i + 1]);
        segments.push((lo, hi, alpha, k * power_integral(lo, hi, alpha)));
    }
    let total: f64 = segments.iter().map(|s| s.3).sum();
    let mut u = d.u() * total;
    let v = d.u();
    for &(lo, hi, alpha, w) in &segments {
        if u <= w {
            return invert_power(lo, hi, alpha, v);
        }
        u -= w;
    }
    imf.breaks[imf.breaks.len() - 1]
}

/// ∫ m^(−α) dm from `lo` to `hi`.
fn power_integral(lo: f64, hi: f64, alpha: f64) -> f64 {
    if (alpha - 1.0).abs() < 1e-9 {
        (hi / lo).ln()
    } else {
        (hi.powf(1.0 - alpha) - lo.powf(1.0 - alpha)) / (1.0 - alpha)
    }
}

/// The mass at which ∫ m^(−α) from `lo` reaches fraction `v` of its total up to `hi`.
fn invert_power(lo: f64, hi: f64, alpha: f64, v: f64) -> f64 {
    if (alpha - 1.0).abs() < 1e-9 {
        lo * (hi / lo).powf(v)
    } else {
        let p = 1.0 - alpha;
        (lo.powf(p) + v * (hi.powf(p) - lo.powf(p))).powf(1.0 / p)
    }
}

/// Interpolated in log mass between `(mass, value)` points.
fn by_log_mass(m: f64, table: &[(f64, f64)]) -> f64 {
    let x = m.ln();
    let first = table[0];
    if x <= first.0.ln() {
        return first.1;
    }
    for w in table.windows(2) {
        let (a, b) = (w[0], w[1]);
        if x <= b.0.ln() {
            let t = (x - a.0.ln()) / (b.0.ln() - a.0.ln());
            return a.1 + t * (b.1 - a.1);
        }
    }
    table[table.len() - 1].1
}

/// Share of stars of initial mass `m` with at least one companion.
#[must_use]
pub fn multiple_fraction(m: f64) -> f64 {
    by_log_mass(m, &[(0.08, 0.20), (0.3, 0.25), (1.0, 0.46), (2.0, 0.55), (5.0, 0.65), (10.0, 0.80), (20.0, 0.90), (60.0, 0.95)])
}

/// Log-normal period distribution (log₁₀ days): mean and spread by primary mass.
fn period_law(m: f64) -> (f64, f64) {
    let mean = by_log_mass(m, &[(0.1, 4.3), (0.5, 4.7), (1.0, 5.03), (3.0, 4.5), (10.0, 3.0), (30.0, 1.8)]);
    let spread = by_log_mass(m, &[(0.1, 1.3), (0.5, 1.8), (1.0, 2.28), (3.0, 2.4), (10.0, 2.2), (30.0, 1.8)]);
    (mean, spread)
}

fn log_period_days(a_au: f64, total_mass: f64) -> f64 {
    ((a_au.powi(3) / total_mass).sqrt() * 365.25).log10()
}

/// A companion to a primary of mass `primary` (M☉) orbiting `central` M☉; an outer one
/// (`inner_log_p` given) orbits at least eight times slower than the inner pair, for stability.
fn companion(d: &mut Draw, primary: f64, central: f64, inner_log_p: Option<f64>) -> Companion {
    let (mean, spread) = period_law(primary);
    let mut log_p = d.normal(mean, spread).clamp(0.2, 8.0);
    if let Some(inner) = inner_log_p {
        log_p = log_p.max(inner + 0.9).min(8.5);
    }
    let q_min = (0.08 / primary).clamp(0.1, 1.0);
    let initial_mass = primary * d.range(q_min, 1.0);
    let p_yr = 10f64.powf(log_p) / 365.25;
    let a_au = (p_yr * p_yr * (central + initial_mass)).cbrt();
    let e = if log_p < 12f64.log10() { 0.0 } else { d.range(0.0, 0.8) };
    Companion { initial_mass, a_au, e }
}

#[cfg(test)]
mod tests {
    use super::*;
    use geoforge_core::body::StellarStage;

    fn sample(galaxy: &Galaxy, n: u64) -> Vec<Phenotype> {
        (0..n).map(|s| draw(galaxy, Seed(s))).collect()
    }

    #[test]
    fn deterministic() {
        let g = Galaxy::milky_way();
        assert_eq!(draw(&g, Seed(9)), draw(&g, Seed(9)));
    }

    /// The Milky Way near the Sun, star by star: about three quarters are red dwarfs, 5–8%
    /// white dwarfs (4.8×10⁻³ pc⁻³ of ~0.08), giants and blue stars rare, neutron stars
    /// ~0.5% and black holes ~0.2% of stars ever born (10⁸–10⁹ of each in the Galaxy).
    #[test]
    fn milky_way_census() {
        let ph = sample(&Galaxy::milky_way(), 20_000);
        let n = ph.len() as f64;
        let share = |f: &dyn Fn(&Star) -> bool| ph.iter().filter(|p| f(&p.primary())).count() as f64 / n;
        let red_dwarfs = share(&|s| s.stage == StellarStage::MainSequence && s.initial_mass() < 0.6);
        let white_dwarfs = share(&|s| s.stage == StellarStage::WhiteDwarf);
        let giants = share(&|s| matches!(s.stage, StellarStage::Giant | StellarStage::Subgiant));
        let blue = share(&|s| s.stage == StellarStage::MainSequence && s.temperature_k > 10_000.0);
        let neutron = share(&|s| s.stage == StellarStage::NeutronStar);
        let holes = share(&|s| s.stage == StellarStage::BlackHole);
        assert!((0.70..0.82).contains(&red_dwarfs), "red dwarfs {red_dwarfs}");
        assert!((0.04..0.09).contains(&white_dwarfs), "white dwarfs {white_dwarfs}");
        assert!((0.003..0.03).contains(&giants), "giants {giants}");
        assert!(blue < 0.003, "B stars {blue}");
        assert!((0.002..0.008).contains(&neutron), "neutron stars {neutron}");
        assert!((0.0005..0.004).contains(&holes), "black holes {holes}");
        // Multiplicity by primary type.
        let mf = |lo: f64, hi: f64| {
            let set: Vec<_> = ph.iter().filter(|p| (lo..hi).contains(&p.primary_mass)).collect();
            set.iter().filter(|p| p.multiplicity() > 1).count() as f64 / set.len() as f64
        };
        assert!((0.2..0.3).contains(&mf(0.08, 0.5)), "M dwarf multiplicity {}", mf(0.08, 0.5));
        assert!((0.4..0.52).contains(&mf(0.8, 1.2)), "Sun-like multiplicity {}", mf(0.8, 1.2));
        // Metallicity near the Sun: centred a little under solar, ~0.2 dex wide.
        let local: Vec<f64> = ph.iter().filter(|p| (7.2..9.2).contains(&p.site.radius_kpc)).map(|p| p.metallicity_dex).collect();
        let mean = local.iter().sum::<f64>() / local.len() as f64;
        assert!((-0.2..0.0).contains(&mean), "local mean [Fe/H] {mean}");
    }

    #[test]
    fn galaxies_differ() {
        let old = sample(&crate::galaxy::GalaxyPreset::Elliptical.galaxy(), 5_000);
        let young = sample(&crate::galaxy::GalaxyPreset::Starburst.galaxy(), 5_000);
        let blue = |set: &[Phenotype]| set.iter().filter(|p| p.primary().temperature_k > 10_000.0 && p.primary().stage == StellarStage::MainSequence).count();
        assert_eq!(blue(&old), 0, "no blue stars in an old elliptical");
        assert!(blue(&young) > 5, "a starburst has blue stars: {}", blue(&young));
        let wd = |set: &[Phenotype]| set.iter().filter(|p| p.primary().stage == StellarStage::WhiteDwarf).count() as f64 / set.len() as f64;
        assert!(wd(&old) > 1.3 * wd(&sample(&Galaxy::milky_way(), 5_000)), "ellipticals are rich in white dwarfs");
        let dwarf = sample(&crate::galaxy::GalaxyPreset::Dwarf.galaxy(), 2_000);
        assert!(dwarf.iter().map(|p| p.metallicity_dex).sum::<f64>() / 2_000.0 < -0.6);
    }

    #[test]
    fn imf_inverts() {
        // Kroupa: 76% of stars from 0.08 to 150 M☉ are under 0.5 M☉.
        let mut d = Draw::new(Seed(1));
        let low = (0..20_000).filter(|_| draw_mass(&mut d, &Imf::kroupa()) < 0.5).count() as f64 / 20_000.0;
        assert!((low - 0.761).abs() < 0.015, "{low}");
    }
}
