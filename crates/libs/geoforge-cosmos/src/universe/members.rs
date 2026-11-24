//! The galaxies of a halo: its central and its satellites, with their types, sizes and light.
//!
//! - **Stars.** A galaxy's stellar mass comes from its halo's (at infall, for a satellite) by
//!   Moster, Naab & White's (2013) relation: 2N[(M/M₁)^−β + (M/M₁)^γ]^−1 of it, M₁ =
//!   10^11.59 M☉, N = 0.0351, β = 1.376, γ = 0.608 — 3% for a Milky Way's halo, ~10⁻⁴ for a
//!   10¹⁰ M☉ dwarf's — with 0.15 dex of scatter.
//! - **Satellites.** Subhalos by infall mass: N(>m) = 0.05 (m/M)^−0.9 down to 10⁹ M☉, those
//!   over a tenth of the host mostly sunk into it by dynamical friction (kept with chance
//!   e^−(m/M)/0.15). A Milky Way's halo then holds about a dozen satellites with 10⁵ M☉ of
//!   stars or more, and one in seven an LMC (Liu et al. 2011: 11% one, 3.5% two). They lie
//!   inside the virial radius (200 times the critical density) as an NFW profile of
//!   concentration 5, less concentrated than the dark matter, as satellites are.
//! - **Types.** Dwarf satellites are spheroidals; field dwarfs irregular. Above ~3×10⁹ M☉ the
//!   early-type share rises with mass (half at 6×10¹⁰ M☉; Bamford et al. 2009), more so for
//!   satellites in groups and clusters (the morphology–density relation, Dressler 1980);
//!   early types are ellipticals (more of them among the massive) or lenticulars. Spirals
//!   run from Sd for the light to Sa for the heavy; three in five are barred; 2% of the
//!   star-forming are bursting.
//! - **Sizes and light.** Half-light radii by Shen et al. (2003): late types 0.1 kpc
//!   M^0.14 (1 + M/3.98×10¹⁰)^0.25, early types 2.88×10⁻⁶ kpc M^0.56, dwarf spheroidals 0.3 kpc
//!   (M/10⁶)^0.3 (Fornax 0.7 kpc, Sculptor 0.3), 0.15 dex of scatter. Light from V-band
//!   mass-to-light ratios typical of each type, 4 for an elliptical to 0.8 for a dwarf
//!   irregular and 0.5 in a starburst. Disks are flat (15–40% as thick as wide), ellipticals
//!   rounder; each turned at random.

use super::halos::{Halo, LIGHTEST};
use crate::draw::Draw;
use crate::galaxy::presets::{Morphology, SpiralStage};
use geoforge_core::galactic::GalaxyId;
use geoforge_core::seed::Seed;

/// Critical density for h = 0.7, M☉ per Mpc³.
const RHO_CRIT: f64 = 1.36e11;
/// Subhalos above m per host: this times (m/M)^−0.9.
const SUBHALO_NORM: f64 = 0.05;
/// Satellites heavier than this share of their host mostly sink.
const FRICTION_SHARE: f64 = 0.15;
/// Concentration of the satellites' NFW profile.
const SATELLITE_CONCENTRATION: f64 = 5.0;
/// Inside about this share of the virial radius the central's tides (a disk's above all)
/// tear satellites apart: one survives at radius r with chance 1 − exp(−(r/this)²), so
/// a Milky Way keeps few within 20 kpc (Garrison-Kimmel et al. 2017).
const TIDAL_SHARE: f64 = 0.1;

/// A galaxy in the universe.
#[derive(Debug, Clone, PartialEq)]
pub struct UniverseGalaxy {
    pub id: GalaxyId,
    pub seed: Seed,
    pub position_mpc: [f64; 3],
    /// Its halo's mass (at infall for a satellite), M☉.
    pub halo_msun: f64,
    /// Its host halo's mass, M☉ (its own for a central).
    pub host_msun: f64,
    pub stellar_msun: f64,
    pub morphology: Morphology,
    /// Half-light radius, kpc.
    pub half_light_kpc: f64,
    /// V-band mass-to-light ratio, M☉/L☉.
    pub mass_to_light: f64,
    /// Its disk's (or short) axis, a unit vector in the universe's frame.
    pub axis: [f64; 3],
    /// How thick it is against its width.
    pub thickness: f64,
}

impl UniverseGalaxy {
    #[must_use]
    pub fn is_central(&self) -> bool {
        self.id.member == 0
    }

    /// V-band luminosity, L☉.
    #[must_use]
    pub fn luminosity_v(&self) -> f64 {
        self.stellar_msun / self.mass_to_light
    }

    /// Absolute V magnitude.
    #[must_use]
    pub fn absolute_v(&self) -> f64 {
        4.83 - 2.5 * self.luminosity_v().log10()
    }

    #[must_use]
    pub fn name(&self) -> String {
        geoforge_names::name(self.seed.0)
    }

    /// Hubble type.
    #[must_use]
    pub fn kind(&self) -> String {
        self.morphology.label()
    }
}

/// Stellar mass (M☉) of a halo of `m` M☉, without scatter.
#[must_use]
pub fn stellar_mass(m: f64) -> f64 {
    let x = m / 10f64.powf(11.59);
    m * 2.0 * 0.0351 / (x.powf(-1.376) + x.powf(0.608))
}

/// Virial radius (Mpc) of a halo of `m` M☉.
#[must_use]
pub fn virial_radius_mpc(m: f64) -> f64 {
    (3.0 * m / (4.0 * std::f64::consts::PI * 200.0 * RHO_CRIT)).cbrt()
}

/// The central galaxy of `halo`.
#[must_use]
pub fn central(universe: Seed, halo: &Halo) -> UniverseGalaxy {
    let id = GalaxyId { bin: halo.bin, cube: halo.cube, n: halo.n, member: 0 };
    galaxy(universe, id, halo.mass_msun, halo.mass_msun, halo.position_mpc)
}

/// The galaxies of `halo` in a universe seeded `universe`: its central first.
#[must_use]
pub fn of_halo(universe: Seed, halo: &Halo) -> Vec<UniverseGalaxy> {
    let id = |member: u32| GalaxyId { bin: halo.bin, cube: halo.cube, n: halo.n, member };
    let mut out = vec![galaxy(universe, id(0), halo.mass_msun, halo.mass_msun, halo.position_mpc)];
    let mut d = Draw::new(id(0).seed(universe).child(b"satellites"));
    let x_min = LIGHTEST / halo.mass_msun;
    if x_min >= 0.5 {
        return out;
    }
    let candidates = d.poisson(SUBHALO_NORM * x_min.powf(-0.9));
    let r_vir = virial_radius_mpc(halo.mass_msun);
    let mut member = 1;
    for _ in 0..candidates {
        let x = x_min * d.u().max(1e-12).powf(-1.0 / 0.9);
        let (sinks, r, dir) = (d.u() >= (-x / FRICTION_SHARE).exp(), nfw_radius(&mut d, SATELLITE_CONCENTRATION), d.direction());
        let torn = d.u() < (-(r / TIDAL_SHARE).powi(2)).exp();
        if x > 0.5 || sinks || torn {
            continue;
        }
        let p = [0, 1, 2].map(|k| halo.position_mpc[k] + r * r_vir * dir[k]);
        out.push(galaxy(universe, id(member), x * halo.mass_msun, halo.mass_msun, p));
        member += 1;
    }
    out
}

/// A radius (in virial radii) drawn from an NFW profile of concentration `c`.
fn nfw_radius(d: &mut Draw, c: f64) -> f64 {
    let m = |x: f64| (1.0 + c * x).ln() - c * x / (1.0 + c * x);
    let target = d.u() * m(1.0);
    let (mut lo, mut hi) = (0.0, 1.0);
    for _ in 0..40 {
        let mid = 0.5 * (lo + hi);
        if m(mid) < target {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

fn galaxy(universe: Seed, id: GalaxyId, halo: f64, host: f64, position_mpc: [f64; 3]) -> UniverseGalaxy {
    let seed = id.seed(universe);
    let mut d = Draw::new(seed.child(b"galaxy"));
    let stellar = stellar_mass(halo) * 10f64.powf(0.15 * d.gauss());
    let central = id.member == 0;
    let morphology = morphology(&mut d, stellar, central, host);
    let scatter = 10f64.powf(0.15 * d.gauss());
    let half_light_kpc = scatter
        * match morphology {
            Morphology::Elliptical | Morphology::Lenticular => 2.88e-6 * stellar.powf(0.56),
            Morphology::DwarfSpheroidal => 0.3 * (stellar / 1e6).powf(0.3),
            _ => 0.1 * stellar.powf(0.14) * (1.0 + stellar / 3.98e10).powf(0.25),
        };
    let (mass_to_light, thickness) = match morphology {
        Morphology::Elliptical => (4.0, d.range(0.5, 0.9)),
        Morphology::Lenticular => (3.0, 0.3),
        Morphology::Spiral { stage, .. } => (
            match stage {
                SpiralStage::A => 2.5,
                SpiralStage::B => 2.2,
                SpiralStage::Bc => 1.9,
                SpiralStage::C => 1.6,
                SpiralStage::D => 1.2,
            },
            0.15,
        ),
        Morphology::DwarfIrregular => (0.8, 0.4),
        Morphology::DwarfSpheroidal => (2.0, 0.6),
        Morphology::Starburst => (0.5, 0.25),
    };
    UniverseGalaxy { id, seed, position_mpc, halo_msun: halo, host_msun: host, stellar_msun: stellar, morphology, half_light_kpc, mass_to_light, axis: d.direction(), thickness }
}

fn morphology(d: &mut Draw, stellar: f64, central: bool, host: f64) -> Morphology {
    let log_m = stellar.log10();
    let in_group = !central && host > 1e13;
    if log_m < 8.0 {
        return if central && d.u() > 0.15 { Morphology::DwarfIrregular } else { Morphology::DwarfSpheroidal };
    }
    if log_m < 9.5 {
        return if in_group && d.u() < 0.5 || !central && d.u() < 0.3 { Morphology::DwarfSpheroidal } else { Morphology::DwarfIrregular };
    }
    let mut early = 1.0 / (1.0 + (-(log_m - 10.8) / 0.25).exp());
    if in_group {
        early = (early + 0.35 * ((host.log10() - 13.0) / 2.0).clamp(0.0, 1.0)).min(0.95);
    }
    if d.u() < early {
        return if d.u() < if log_m > 11.0 { 0.85 } else { 0.6 } { Morphology::Elliptical } else { Morphology::Lenticular };
    }
    if d.u() < 0.02 {
        return Morphology::Starburst;
    }
    let stage = match log_m {
        m if m < 9.8 => SpiralStage::D,
        m if m < 10.2 => SpiralStage::C,
        m if m < 10.6 => SpiralStage::Bc,
        m if m < 10.9 => SpiralStage::B,
        _ => SpiralStage::A,
    };
    Morphology::Spiral { stage, barred: d.u() < 0.6 }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn milky_way_halo(n: u32) -> Halo {
        Halo { bin: 6, cube: [0, 0, 0], n, mass_msun: 1.3e12, position_mpc: [0.0; 3] }
    }

    #[test]
    fn stellar_masses_follow_moster() {
        let mw = stellar_mass(1.3e12);
        assert!((3e10..6e10).contains(&mw), "Milky Way {mw:e}");
        assert!(stellar_mass(1e10) < 1e7 && stellar_mass(1e15) / 1e15 < 2e-3);
        assert!((0.2..0.26).contains(&virial_radius_mpc(1.3e12)), "R200 {}", virial_radius_mpc(1.3e12));
    }

    /// A Milky Way's halo: about a dozen satellites with 10⁵ M☉ of stars, most of them
    /// spheroidals within ~100 kpc and hardly any within 20 kpc (the Milky Way: only the
    /// Sagittarius dwarf, being torn apart); one halo in six to ten with an LMC.
    #[test]
    fn a_milky_way_has_a_dozen_satellites() {
        let (mut bright, mut lmcs, mut distances) = (0.0, 0.0, Vec::new());
        let n = 400;
        for k in 0..n {
            let g = of_halo(Seed(3), &milky_way_halo(k));
            let sats = &g[1..];
            bright += sats.iter().filter(|s| s.stellar_msun > 1e5).count() as f64;
            lmcs += f64::from(u8::from(sats.iter().any(|s| s.stellar_msun > 1e9)));
            distances.extend(sats.iter().filter(|s| s.stellar_msun > 1e5).map(|s| halos_dist(s.position_mpc) * 1000.0));
            assert!(sats.iter().filter(|s| s.stellar_msun < 1e7).all(|s| s.morphology == Morphology::DwarfSpheroidal));
        }
        let (bright, lmcs) = (bright / f64::from(n), lmcs / f64::from(n));
        distances.sort_by(f64::total_cmp);
        let median = distances[distances.len() / 2];
        assert!((7.0..20.0).contains(&bright), "{bright} satellites over 10⁵ M☉");
        assert!((0.08..0.25).contains(&lmcs), "{lmcs} with an LMC");
        assert!((50.0..130.0).contains(&median), "median satellite distance {median} kpc");
        let inner = distances.iter().filter(|&&r| r < 20.0).count() as f64 / distances.len() as f64;
        assert!(inner < 0.03, "{inner} of satellites within 20 kpc");
    }

    fn halos_dist(p: [f64; 3]) -> f64 {
        (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt()
    }

    #[test]
    fn heavy_galaxies_are_early_types() {
        let share = |m: f64| {
            let halo = |n| Halo { bin: 10, cube: [0, 0, 1], n, mass_msun: m, position_mpc: [0.0; 3] };
            let g: Vec<UniverseGalaxy> = (0..600).map(|n| of_halo(Seed(4), &halo(n)).swap_remove(0)).collect();
            g.iter().filter(|g| matches!(g.morphology, Morphology::Elliptical | Morphology::Lenticular)).count() as f64 / 600.0
        };
        let (mw, group) = (share(1.3e12), share(3e13));
        assert!(mw < 0.4 && group > 0.7, "early-type centrals: {mw} for a Milky Way's halo, {group} for a group's");
        let mw = of_halo(Seed(5), &milky_way_halo(1)).swap_remove(0);
        assert!((2.0..9.0).contains(&mw.half_light_kpc) && (-22.5..-19.5).contains(&mw.absolute_v()), "{} kpc, M_V {}", mw.half_light_kpc, mw.absolute_v());
    }
}
