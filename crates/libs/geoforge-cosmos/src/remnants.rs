//! What a dead star's remnant does: a neutron star's spin and field, a black hole's horizon,
//! and either feeding on a close companion.
//!
//! - **Neutron stars** are born spinning every 0.3 ± 0.15 s with fields of 10^(12.65 ± 0.55) G
//!   (Faucher-Giguère & Kaspi 2006) and spin down as magnetic dipoles: P·Ṗ = (B / 3.2×10¹⁹ G)²,
//!   so P² = P₀² + 2(B / 3.2×10¹⁹)²·t (the Crab, 970 years old at 3.8×10¹² G: 33 ms, Ṗ =
//!   4.2×10⁻¹³, 4.5×10³⁸ erg/s). One shines as a radio pulsar until it crosses the death line,
//!   B/P² = 0.17×10¹² G s⁻² (Chen & Ruderman 1993), after some 10–100 Myr; one with a field
//!   over the quantum limit (4.4×10¹³ G) younger than 100,000 years is a magnetar. A neutron
//!   star whose companion has given it mass (now a white dwarf) or that holds planets (PSR
//!   B1257+12's) has been spun up: a millisecond pulsar, 1.5–10 ms at 10⁸–10⁹ G, shining for
//!   billions of years.
//! - **Black holes**: the horizon, 2GM/c² = 2.95 km per M☉.
//! - **Feeding**: a remnant whose companion fills its Roche lobe at periastron draws its gas
//!   through a disk (a low-mass X-ray binary, Sco X-1); one beside an O or B supergiant
//!   within 1 AU catches its wind (Cygnus X-1). Their X-rays are not counted in the light.

use crate::binary::roche_lobe_share;
use crate::draw::Draw;
use crate::stars::end_of_life_gyr;
use geoforge_core::body::{Elements, NodeKind, SolarSystem, Star, StellarStage};
use geoforge_core::seed::Seed;
use geoforge_core::Key;

/// Seconds in a gigayear.
const GYR_S: f64 = 3.156e16;
/// The magnetic-dipole constant: B = 3.2×10¹⁹ G √(P·Ṗ).
const DIPOLE_G: f64 = 3.2e19;
/// A neutron star's moment of inertia, g cm².
const MOMENT_G_CM2: f64 = 1e45;
/// The death line: radio pulsing stops below this B/P², G s⁻².
const DEATH_LINE: f64 = 0.17e12;
/// The quantum critical field, G: magnetars lie above it.
const MAGNETAR_G: f64 = 4.4e13;
/// How long a magnetar stays active, Gyr.
const MAGNETAR_GYR: f64 = 1e-4;
const L_SUN_ERG_S: f64 = 3.828e33;
const AU_PER_R_SUN: f64 = 215.032;

/// What a neutron star shows.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PulsarKind {
    /// A young neutron star with a field past the quantum limit: X-ray bursts and flares.
    Magnetar,
    /// Pulsing in radio, young and slowing.
    Radio,
    /// Spun up by its companion: pulsing every few milliseconds.
    Millisecond,
    /// Past the death line: silent.
    Quiet,
}

impl PulsarKind {
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Self::Magnetar => "magnetar",
            Self::Radio => "radio pulsar",
            Self::Millisecond => "millisecond pulsar",
            Self::Quiet => "quiet neutron star",
        }
    }
}

/// A neutron star's spin and field.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Pulsar {
    pub period_s: f64,
    pub field_gauss: f64,
    /// Rotational power lost to its wind and beams, L☉.
    pub spin_down_l_solar: f64,
    pub kind: PulsarKind,
}

/// The spin and field of `star` if it is a neutron star, drawn from `seed`; `recycled` if a
/// companion has given it mass or it holds planets.
#[must_use]
pub fn pulsar(star: &Star, seed: Seed, recycled: bool) -> Option<Pulsar> {
    if star.stage != StellarStage::NeutronStar {
        return None;
    }
    let mut d = Draw::new(seed.child(b"pulsar"));
    let age_gyr = (star.age_gyr - end_of_life_gyr(star.initial_mass())).max(1e-9);
    let (period_s, field_gauss) = if recycled {
        (10f64.powf(d.range(0.18, 1.0)) * 1e-3, 10f64.powf(d.range(8.0, 9.0)))
    } else {
        let p0 = d.normal(0.3, 0.15).max(0.01);
        let b = 10f64.powf(d.normal(12.65, 0.55));
        ((p0 * p0 + 2.0 * (b / DIPOLE_G).powi(2) * age_gyr * GYR_S).sqrt(), b)
    };
    let p_dot = (field_gauss / DIPOLE_G).powi(2) / period_s;
    let spin_down_l_solar = 4.0 * std::f64::consts::PI.powi(2) * MOMENT_G_CM2 * p_dot / period_s.powi(3) / L_SUN_ERG_S;
    let kind = if recycled {
        PulsarKind::Millisecond
    } else if field_gauss > MAGNETAR_G && age_gyr < MAGNETAR_GYR {
        PulsarKind::Magnetar
    } else if field_gauss / (period_s * period_s) >= DEATH_LINE {
        PulsarKind::Radio
    } else {
        PulsarKind::Quiet
    };
    Some(Pulsar { period_s, field_gauss, spin_down_l_solar, kind })
}

/// A black hole's horizon (Schwarzschild radius), km.
#[must_use]
pub fn horizon_km(star: &Star) -> Option<f64> {
    (star.stage == StellarStage::BlackHole).then_some(2.953 * star.mass_solar)
}

/// How a remnant feeds on its companion.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Feeding {
    /// The companion overfills its Roche lobe: a disk, a low-mass X-ray binary.
    RocheLobe,
    /// A hot supergiant's wind, caught close in: a high-mass X-ray binary.
    Wind,
}

impl Feeding {
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Self::RocheLobe => "X-ray binary, fed through a disk",
            Self::Wind => "X-ray binary, fed by its companion's wind",
        }
    }
}

/// Whether `remnant` (a neutron star or black hole) feeds on `companion`, `a_au` apart with
/// eccentricity `e`.
#[must_use]
pub fn feeding(remnant: &Star, companion: &Star, a_au: f64, e: f64) -> Option<Feeding> {
    if !matches!(remnant.stage, StellarStage::NeutronStar | StellarStage::BlackHole) || companion.stage.is_remnant() {
        return None;
    }
    let lobe_r_sun = roche_lobe_share(companion.mass_solar / remnant.mass_solar) * a_au * (1.0 - e) * AU_PER_R_SUN;
    if companion.radius_solar >= lobe_r_sun {
        Some(Feeding::RocheLobe)
    } else if a_au < 1.0 && companion.temperature_k > 10_000.0 && companion.luminosity_solar > 1e4 {
        Some(Feeding::Wind)
    } else {
        None
    }
}

/// Lines on the remnant at `key` of `system` (none for a living star): a neutron star's
/// spin, field and kind, a black hole's horizon, and whether it feeds on a companion.
#[must_use]
pub fn notes(system: &SolarSystem, key: Key) -> Vec<String> {
    let Some(node) = system.node(key) else { return Vec::new() };
    let NodeKind::Star(star) = &node.node else { return Vec::new() };
    let mut out = Vec::new();
    // Its partners and their separations: one orbiting the other, or both a barycentre.
    let partners: Vec<(&Star, f64, f64)> = system
        .nodes
        .iter()
        .filter(|n| n.key != key)
        .filter_map(|n| {
            let NodeKind::Star(s) = &n.node else { return None };
            let sep = if n.parent == Some(key) {
                n.elements?
            } else if node.parent == Some(n.key) {
                node.elements?
            } else if n.parent.is_some() && n.parent == node.parent && system.node(n.parent?).is_some_and(|p| matches!(p.node, NodeKind::Barycenter)) {
                let (a, b) = (n.elements?, node.elements?);
                Elements { a_au: a.a_au + b.a_au, ..a }
            } else {
                return None;
            };
            Some((s, sep.a_au, sep.e))
        })
        .collect();
    if star.stage == StellarStage::NeutronStar {
        let planets = system.children(key).iter().any(|c| matches!(c.node, NodeKind::Body(_)));
        let recycled = planets || partners.iter().any(|(s, _, _)| s.stage == StellarStage::WhiteDwarf);
        if let Some(p) = pulsar(star, system.root.seed.child_u64(b"remnant", key.0), recycled) {
            out.push(format!("{}: every {} s, field 10^{:.1} G, spinning down at {:.1e} L☉", p.kind.name(), sig(p.period_s), p.field_gauss.log10(), p.spin_down_l_solar));
        }
    }
    if let Some(r) = horizon_km(star) {
        out.push(format!("black hole: horizon {r:.0} km in radius"));
    }
    for (s, a, e) in partners {
        if let Some(f) = feeding(star, s, a, e) {
            out.push(f.name().to_string());
        }
    }
    out
}

/// `x` to three significant figures.
fn sig(x: f64) -> String {
    let digits = (2 - x.abs().log10().floor() as i32).max(0) as usize;
    format!("{x:.digits$}")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stars::evolve;

    fn neutron(age_after_gyr: f64) -> Star {
        evolve(12.0, end_of_life_gyr(12.0) + age_after_gyr, 0.0)
    }

    /// The Crab, by the spin-down law: 33 ms and 4.5×10³⁸ erg/s at 970 years and 3.8×10¹² G
    /// from a 19 ms birth.
    #[test]
    fn the_crab_spins_down() {
        let (p0, b, t) = (0.019f64, 3.8e12f64, 970.0 * 3.156e7);
        let p = (p0 * p0 + 2.0 * (b / DIPOLE_G).powi(2) * t).sqrt();
        let p_dot = (b / DIPOLE_G).powi(2) / p;
        let e_dot = 4.0 * std::f64::consts::PI.powi(2) * MOMENT_G_CM2 * p_dot / p.powi(3);
        assert!((0.030..0.040).contains(&p) && (3e-13..6e-13).contains(&p_dot) && (2e38..8e38).contains(&e_dot), "{p} {p_dot} {e_dot}");
    }

    /// Young neutron stars pulse; almost all those older than a few hundred million years
    /// have crossed the death line; a few young ones are magnetars; recycled ones pulse in
    /// milliseconds; the same seed gives the same pulsar.
    #[test]
    fn pulsars_age_into_silence() {
        let count = |age: f64, k: PulsarKind| (0..400u64).filter(|&s| pulsar(&neutron(age), Seed(s), false).unwrap().kind == k).count();
        assert!(count(1e-3, PulsarKind::Radio) > 300, "young pulsars");
        assert!(count(1.0, PulsarKind::Quiet) > 390, "old ones silent");
        let magnetars = count(1e-5, PulsarKind::Magnetar);
        assert!((3..60).contains(&magnetars), "{magnetars} magnetars of 400");
        let ms = pulsar(&neutron(3.0), Seed(1), true).unwrap();
        assert!(ms.kind == PulsarKind::Millisecond && (1.5e-3..1e-2).contains(&ms.period_s));
        assert_eq!(pulsar(&neutron(0.1), Seed(9), false), pulsar(&neutron(0.1), Seed(9), false));
        assert!(pulsar(&evolve(1.0, 4.6, 0.0), Seed(1), false).is_none());
    }

    #[test]
    fn horizons_and_feeding() {
        let hole = evolve(40.0, 1.0, 0.0);
        assert_eq!(hole.stage, StellarStage::BlackHole);
        assert!((horizon_km(&hole).unwrap() / hole.mass_solar - 2.95).abs() < 0.01);
        // A swelling subgiant (2.8 R☉) 0.02 AU from a neutron star overflows onto it; a
        // dwarf far off does not.
        let sub = evolve(1.2, 6.0, 0.0);
        let ns = neutron(1.0);
        assert_eq!(feeding(&ns, &sub, 0.02, 0.0), Some(Feeding::RocheLobe), "{} R☉", sub.radius_solar);
        assert_eq!(feeding(&ns, &evolve(0.5, 3.0, 0.0), 2.0, 0.0), None);
        // An O supergiant close to a black hole feeds it its wind (Cygnus X-1).
        let o = evolve(30.0, 0.004, 0.0);
        assert_eq!(feeding(&hole, &o, 0.2, 0.0), Some(Feeding::Wind), "{:?} {} L☉", o.stage, o.luminosity_solar);
    }
}
