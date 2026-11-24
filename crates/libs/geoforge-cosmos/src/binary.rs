//! Binary evolution, simply: what a pair of stars close enough to touch becomes.
//!
//! A star that swells past its Roche lobe (Eggleton 1983, at periastron) pours its envelope onto
//! its companion. Which way it goes turns on the mass ratio q = donor/accretor against a
//! critical one (Hurley, Tout & Pols 2002): 4 for a subgiant, and for a giant with a core of
//! `Mc` in a star of `M`, (1.67 − 0.3 + 2(Mc/M)⁵)/2.13, about 0.65.
//!
//! - **Stable transfer** (q under the critical ratio): the donor gives up its whole envelope and
//!   is left as its core, a white dwarf (a helium one under ~0.45 M☉). Half of what it gives
//!   reaches the companion; the rest leaves the pair from beside it, carrying its orbital
//!   momentum (isotropic re-emission), which sets the new separation. The companion takes the
//!   mass and the fresh hydrogen with it, and is rejuvenated: it keeps the share of its main
//!   sequence it had lived (Hurley et al. 2002), now at its new mass (a blue straggler).
//! - **Common envelope** (q over it): the envelope swallows the companion, and the drag that
//!   spirals the pair in spends orbital energy to blow it off, with an efficiency αλ = 0.5
//!   (Webbink 1984): α·[G·Mc·Ma/2a_f − G·M·Ma/2a_i] = G·M·M_env/(λ·R). What is left is a tight
//!   pair of the core and the companion, unless the companion would then overfill its own lobe,
//!   when the two merge.
//! - **Merging**: one star, of both masses less half the envelope, as far through its life as
//!   the donor was. On the main sequence a donor merges only when q is over 3; below that the
//!   pair would be a contact binary, which is not modelled (the stars live on as if apart).
//!
//! The core of a giant comes from its brightness, which a shell-burning star's core alone sets:
//! L ≈ 2.3×10⁵ L☉·(Mc/M☉)⁶ (Refsdal & Weigert 1970; 0.47 M☉ at the Sun's tip). Tides round the
//! orbit while the stars touch. Donors born at 8 M☉ or more, whose cores go on to collapse, and
//! a second transfer (the rejuvenated companion later filling its lobe onto the white dwarf: a
//! cataclysmic variable, a type Ia supernova) are not modelled. Stars a galaxy's sky shows stay
//! as if alone: a pair's members are worked out with this only where its system is built.

use crate::stars::{classify, evolve, main_sequence_lifetime_gyr, white_dwarf_luminosity, white_dwarf_mass, white_dwarf_radius, HYDROGEN_BURNING_LIMIT, WHITE_DWARF_LIMIT};
use geoforge_core::body::{Star, StellarStage};
use geoforge_core::derived::effective_temperature_k;

/// Solar radii in an AU.
const R_SUN_PER_AU: f64 = 215.032;
/// Share of the donor's envelope that reaches the companion in stable transfer.
const ACCRETED_SHARE: f64 = 0.5;
/// Common-envelope efficiency times the envelope's binding-energy factor, αλ.
const ALPHA_LAMBDA: f64 = 0.5;
/// Critical mass ratio of a main-sequence donor, past which it merges.
const MAIN_SEQUENCE_Q_CRIT: f64 = 3.0;
/// Critical mass ratio of a subgiant donor.
const SUBGIANT_Q_CRIT: f64 = 4.0;
/// Share of the donor's envelope that a merger throws off.
const MERGER_LOSS_SHARE: f64 = 0.5;

/// What a pair has been through.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum History {
    /// Never touched: each star lived as if alone.
    Detached,
    /// The donor's envelope went over to the companion, `at_gyr` after birth.
    Stable { at_gyr: f64 },
    /// The donor's envelope swallowed the companion and was thrown off.
    CommonEnvelope { at_gyr: f64 },
    /// The two became one.
    Merged { at_gyr: f64 },
}

/// A pair as it is now: the first-named star, the second (gone if they merged), and their orbit.
#[derive(Debug, Clone, PartialEq)]
pub struct PairNow {
    pub primary: Star,
    pub companion: Option<Star>,
    /// Today's separation, AU.
    pub a_au: f64,
    pub e: f64,
    pub history: History,
}

/// Roche lobe radius over the separation of a star of mass ratio `q` = its mass / the other's
/// (Eggleton 1983: within 1% for all q).
#[must_use]
pub fn roche_lobe_share(q: f64) -> f64 {
    let q23 = q.powf(2.0 / 3.0);
    0.49 * q23 / (0.6 * q23 + (1.0 + q.cbrt()).ln())
}

/// Core mass (M☉) of a star past its main sequence shining with `luminosity` L☉, at most the
/// white dwarf a star born with `initial_mass` leaves.
#[must_use]
pub fn core_mass(luminosity: f64, initial_mass: f64) -> f64 {
    (luminosity / 2.3e5).powf(1.0 / 6.0).clamp(0.1, white_dwarf_mass(initial_mass))
}

/// The pair of stars born with `m1` and `m2` M☉, `a_au` apart with eccentricity `e`, after
/// `age_gyr`, of metallicity `feh`.
#[must_use]
pub fn evolve_pair(m1: f64, m2: f64, a_au: f64, e: f64, age_gyr: f64, feh: f64) -> PairNow {
    let detached = || {
        let (p, c) = (evolve(m1, age_gyr, feh), evolve(m2, age_gyr, feh));
        // Mass lost in winds widens the orbit in proportion.
        let a = a_au * (m1 + m2) / (p.mass_solar + c.mass_solar);
        PairNow { primary: p, companion: Some(c), a_au: a, e, history: History::Detached }
    };
    let first_donates = m1 >= m2;
    let (md0, ma0) = if first_donates { (m1, m2) } else { (m2, m1) };
    if md0 >= WHITE_DWARF_LIMIT || ma0 < HYDROGEN_BURNING_LIMIT {
        return detached();
    }
    let Some((t, donor, acc, a_i)) = overflow(md0, ma0, a_au, e, age_gyr, feh) else { return detached() };
    let order = |d: Star, a: Option<Star>| if first_donates { (d, a) } else { (a.expect("the pair kept its donor's partner"), Some(d)) };
    let q = donor.mass_solar / acc.mass_solar;
    let later = age_gyr - t;
    if donor.stage == StellarStage::MainSequence {
        if q <= MAIN_SEQUENCE_Q_CRIT {
            return detached();
        }
        let m = donor.mass_solar + acc.mass_solar;
        return merged(m, t / main_sequence_lifetime_gyr(md0), later, feh, t);
    }
    let mc = core_mass(donor.luminosity_solar, md0).min(donor.mass_solar);
    let env = donor.mass_solar - mc;
    let q_crit = if donor.stage == StellarStage::Subgiant { SUBGIANT_Q_CRIT } else { (1.67 - 0.3 + 2.0 * (mc / donor.mass_solar).powi(5)) / 2.13 };
    if q <= q_crit {
        // Stable: the envelope goes, half of it to the companion.
        let ma = acc.mass_solar + ACCRETED_SHARE * env;
        let a_f = a_i * stable_widening(donor.mass_solar, acc.mass_solar, mc);
        let companion = rejuvenated(acc.mass_solar, ma, t, ma0, later, feh);
        let wd = white_dwarf(mc, md0, later, age_gyr, feh);
        let a_now = a_f * (mc + ma) / (wd.mass_solar + companion.mass_solar);
        let (p, c) = order(wd, Some(companion));
        return PairNow { primary: p, companion: c, a_au: a_now, e: 0.0, history: History::Stable { at_gyr: t } };
    }
    // A common envelope: the orbit's energy pays for the envelope's binding.
    let r_au = donor.radius_solar / R_SUN_PER_AU;
    let a_f = mc * acc.mass_solar / (2.0 * donor.mass_solar * env / (ALPHA_LAMBDA * r_au) + donor.mass_solar * acc.mass_solar / a_i);
    if acc.radius_solar > roche_lobe_share(acc.mass_solar / mc) * a_f * R_SUN_PER_AU {
        let m = donor.mass_solar + acc.mass_solar - MERGER_LOSS_SHARE * env;
        return merged(m, t / main_sequence_lifetime_gyr(md0), later, feh, t);
    }
    let companion = rejuvenated(acc.mass_solar, acc.mass_solar, t, ma0, later, feh);
    let wd = white_dwarf(mc, md0, later, age_gyr, feh);
    let a_now = a_f * (mc + acc.mass_solar) / (wd.mass_solar + companion.mass_solar);
    let (p, c) = order(wd, Some(companion));
    PairNow { primary: p, companion: c, a_au: a_now, e: 0.0, history: History::CommonEnvelope { at_gyr: t } }
}

/// When (Gyr) the donor born with `md0` M☉ first fills its Roche lobe at periastron, before
/// `age_gyr`, with both stars then and their separation (AU, widened by winds so far).
fn overflow(md0: f64, ma0: f64, a_au: f64, e: f64, age_gyr: f64, feh: f64) -> Option<(f64, Star, Star, f64)> {
    let t_ms = main_sequence_lifetime_gyr(md0);
    // Coarse through the main sequence, fine through the giant branch, where the star swells.
    let shares = (0..=50).map(|k| f64::from(k) / 50.0).chain((1..=400).map(|k| 1.0 + 0.14 * f64::from(k) / 400.0));
    for f in shares {
        let t = f * t_ms;
        if t > age_gyr {
            return None;
        }
        let (donor, acc) = (evolve(md0, t, feh), evolve(ma0, t, feh));
        let a = a_au * (md0 + ma0) / (donor.mass_solar + acc.mass_solar);
        let lobe = roche_lobe_share(donor.mass_solar / acc.mass_solar) * a * (1.0 - e) * R_SUN_PER_AU;
        if donor.radius_solar > lobe {
            return Some((t, donor, acc, a));
        }
    }
    None
}

/// Factor by which stable transfer widens the orbit: the donor of `md` M☉ shrinks to its core
/// `mc`, the accretor of `ma` takes in [`ACCRETED_SHARE`] of it, and the rest leaves from beside
/// the accretor with its specific orbital momentum. With J = μ√(GMa), a ∝ J²·M/(Md·Ma)², and
/// each bit lost takes dJ/J = dm·Md/(Ma·M).
fn stable_widening(md: f64, ma: f64, mc: f64) -> f64 {
    const STEPS: usize = 200;
    let dm = (md - mc) / STEPS as f64;
    let (mut d, mut a, mut ln_j) = (md, ma, 0.0);
    for _ in 0..STEPS {
        ln_j -= (1.0 - ACCRETED_SHARE) * dm * d / (a * (d + a));
        d -= dm;
        a += ACCRETED_SHARE * dm;
    }
    (2.0 * ln_j).exp() * ((d + a) / (md + ma)) * (md * ma / (d * a)).powi(2)
}

/// A companion of `before` M☉ grown to `after` at `t_gyr`, born with `born` M☉, after `later`
/// Gyr more: as far through its main sequence as it was, at its new mass.
fn rejuvenated(before: f64, after: f64, t_gyr: f64, born: f64, later: f64, feh: f64) -> Star {
    let lived = (t_gyr / main_sequence_lifetime_gyr(born)).min(1.0) * (before / after);
    let mut s = evolve(after, lived * main_sequence_lifetime_gyr(after) + later, feh);
    s.age_gyr = t_gyr + later;
    s
}

/// The core of `mass` M☉ left by a donor born with `born` M☉, cooled for `cooling_gyr`, in a
/// system `age_gyr` old.
fn white_dwarf(mass: f64, born: f64, cooling_gyr: f64, age_gyr: f64, feh: f64) -> Star {
    let (l, r) = (white_dwarf_luminosity(mass, cooling_gyr.max(1e-3)), white_dwarf_radius(mass));
    let t = effective_temperature_k(l, r);
    let stage = StellarStage::WhiteDwarf;
    Star { mass_solar: mass, luminosity_solar: l, radius_solar: r, temperature_k: t, age_gyr, metallicity_dex: feh, spectral_class: classify(stage, t), stage, initial_mass_solar: Some(born) }
}

/// The one star two became at `t_gyr`: `mass` M☉, `share` of the way through its main
/// sequence (the donor's), `later` Gyr ago.
fn merged(mass: f64, share: f64, later: f64, feh: f64, t_gyr: f64) -> PairNow {
    let mut s = evolve(mass, share * main_sequence_lifetime_gyr(mass) + later, feh);
    s.age_gyr = t_gyr + later;
    PairNow { primary: s, companion: None, a_au: 0.0, e: 0.0, history: History::Merged { at_gyr: t_gyr } }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roche_lobes() {
        // Equal masses: 0.38 of the separation; the Earth–Sun (q = 3×10⁻⁶): 0.49·q^⅓.
        assert!((roche_lobe_share(1.0) - 0.379).abs() < 0.002);
        assert!((roche_lobe_share(3e-6) - 0.0071).abs() < 0.0002, "{}", roche_lobe_share(3e-6));
        assert!(roche_lobe_share(2.0) > roche_lobe_share(1.0));
    }

    #[test]
    fn the_suns_core_at_its_tip() {
        assert!((core_mass(2500.0, 1.0) - 0.47).abs() < 0.02);
        assert!(core_mass(3.0, 1.0) < 0.16);
    }

    /// Pairs too far apart to touch live as if alone; their orbit widens as they lose mass.
    #[test]
    fn wide_pairs_stay_detached() {
        let p = evolve_pair(1.0, 0.8, 50.0, 0.3, 12.0, 0.0);
        assert_eq!(p.history, History::Detached);
        assert_eq!(p.primary, evolve(1.0, 12.0, 0.0));
        assert!(p.a_au > 50.0 && p.e == 0.3);
        // Young: nothing has swelled yet.
        assert_eq!(evolve_pair(1.0, 0.8, 0.5, 0.0, 3.0, 0.0).history, History::Detached);
    }

    /// A Sun-like star 0.7 AU from a lighter one fills its lobe as a red giant and throws off its
    /// envelope round its companion: a white dwarf and a red dwarf days apart (post-common-
    /// envelope binaries such as V471 Tau and the Hyades' pairs: 0.01–0.1 AU).
    #[test]
    fn a_giant_swallows_its_companion() {
        let p = evolve_pair(1.0, 0.4, 0.7, 0.0, 12.0, 0.0);
        assert!(matches!(p.history, History::CommonEnvelope { .. }), "{:?}", p.history);
        assert_eq!(p.primary.stage, StellarStage::WhiteDwarf);
        assert!(p.primary.mass_solar < 0.47, "a helium white dwarf: {}", p.primary.mass_solar);
        let c = p.companion.expect("survives");
        assert!(c.stage.is_main_sequence() && (c.mass_solar - 0.4).abs() < 1e-9);
        assert!((0.005..0.15).contains(&p.a_au) && p.e == 0.0, "{} AU", p.a_au);
    }

    /// A subgiant donor gives its envelope over: the companion grows, younger than its age,
    /// and the orbit widens (a blue straggler with a helium white dwarf, as in the old cluster
    /// NGC 188).
    #[test]
    fn stable_transfer_makes_a_blue_straggler() {
        let p = evolve_pair(1.3, 1.25, 0.025, 0.0, 4.7, 0.0);
        assert!(matches!(p.history, History::Stable { .. }), "{:?}", p.history);
        assert!(p.primary.stage == StellarStage::WhiteDwarf && p.primary.mass_solar < 0.45, "{:?}", p.primary);
        let c = p.companion.expect("survives");
        assert!(c.mass_solar > 1.6 && c.stage.is_main_sequence(), "{c:?}");
        assert!(p.a_au > 0.1 && p.e == 0.0, "{} AU", p.a_au);
    }

    /// A star far heavier than its main-sequence partner that touches it while both burn
    /// hydrogen merges; a giant that would leave its partner inside its own lobe merges too.
    #[test]
    fn close_unequal_pairs_merge() {
        let p = evolve_pair(4.0, 0.9, 0.025, 0.0, 0.2, 0.0);
        assert!(matches!(p.history, History::Merged { .. }), "{:?}", p.history);
        assert!(p.companion.is_none() && p.primary.mass_solar > 4.0);
        let g = evolve_pair(2.0, 0.1, 0.3, 0.0, 2.0, 0.0);
        assert!(matches!(g.history, History::Merged { .. }), "{:?}", g.history);
    }

    /// The donor is the heavier star whichever is named first; the names keep their stars.
    #[test]
    fn the_heavier_star_gives() {
        let (a, b) = (evolve_pair(1.0, 0.4, 0.7, 0.0, 12.0, 0.0), evolve_pair(0.4, 1.0, 0.7, 0.0, 12.0, 0.0));
        assert_eq!(a.primary, b.companion.clone().unwrap());
        assert_eq!(b.primary, a.companion.clone().unwrap());
        assert_eq!(a, evolve_pair(1.0, 0.4, 0.7, 0.0, 12.0, 0.0));
    }
}
