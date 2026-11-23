//! Planetary systems for any star: what forms at birth, and what is left when the star has
//! swollen into a giant or died.
//!
//! Planets form on a ladder of orbits from a few days out. Inside the frost line the small
//! planets come as Kepler's compact systems do: a chain of neighbours whose period ratios are
//! log-normal about 1.8 (0.3 in ln; a peak near 1.6 and a long tail, as Kepler's multi-planet
//! systems, Fabrycky et al. 2014), with masses alike within a system ("peas in a pod", Weiss et
//! al. 2018: 0.25 dex about the system's own). Beyond it each rung is 1.45–2.1 times farther
//! than the last. Then they settle ([`crate::dynamics`]): neighbours at least 8 mutual Hill
//! radii apart, too-close or too-eccentric pairs scattered, some scattered giants made hot
//! Jupiters. Belts sit at the giants' resonances. Rates, by the host's initial mass `M` (M☉)
//! and metallicity `[Fe/H]`:
//!
//! | Planets | Rate | After |
//! |---|---|---|
//! | an inner system at all | 70% of red dwarfs, 50% of Sun-like stars, 25% at 2 M☉, 2% at 10; around brown dwarfs falling as √M, with masses in proportion to the host's and no ice giants (assumed) | Kepler: 2.5 small planets per red dwarf within 200 days (Dressing & Charbonneau 2015); 30% of Sun-like stars with Kepler-like systems of ~3 (Zhu et al. 2018) |
//! | planets in the chain | 1 + Poisson(2), 1 + Poisson(2.6) for red dwarfs | |
//! | a planet's mass | the system's, 10^U(−1.2, 1.1) M⊕·√M, ± 0.25 dex; a sub-Neptune above 3 M⊕ | about a third sub-Neptunes around Sun-like stars, fewer around red dwarfs |
//! | giant planets | 0.14·M·10^(1.2[Fe/H]) up to 2.5 M☉, falling as M⁻² above | ∝ M and 10^(1.2[Fe/H]) (Johnson et al. 2010); ~14% of Sun-like stars |
//! | giants, given any | 1 + Poisson(0.7), packed at period ratios 1.5–2 (between the 3:2 and 2:1 resonances disk migration leaves them in, as the Nice model starts the Solar System's) | many scatter: survivors' eccentricities as radial velocities find (Kipping 2013) |
//! | a hot Jupiter by disk migration, given giants | 2% | with the high-eccentricity channel, ~1% of Sun-like stars (Wright et al. 2012: 1.2%; Kepler 0.4–0.6%) |
//! | each outer rung without a giant | 20% an ice giant, 15% a dwarf planet | |
//! | an outer belt, at the outermost giant's 3:2–2:1 (else beyond the planets and 10× the frost line) | 90% | debris disks are seen around ~20% of Sun-like stars, but a Kuiper belt like ours would not be |
//!
//! Giants form on near-circular orbits (e under 0.2); the larger eccentricities come from
//! scattering.
//! A brown dwarf's disk is lit by its glow at 10 Myr.
//!
//! A close companion stunts planets around one star: under 1 AU none form (they may orbit the
//! pair instead), at 10 AU a fifth as many, beyond 200 AU no effect (Moe & Kratter 2021);
//! orbits must also be stable (Holman & Wiegert 1999).
//!
//! After the main sequence: a giant swallows what orbits within three of its radii; a white
//! dwarf keeps only what formed beyond ~2 AU·M^0.8 (Mustill & Villaver 2012), on orbits
//! widened by the mass the star lost; a supernova unbinds everything, and a neutron star has
//! planets only if a second generation formed from fallback (under 0.5% of pulsars, Nițu et
//! al. 2022, like PSR B1257+12's); a black hole has none.

use crate::draw::Draw;
use crate::dynamics;
use geoforge_core::body::{BodyClass, Star, StellarStage};
use geoforge_core::derived::frost_line_au;

/// Solar radius in AU.
const SOLAR_RADIUS_AU: f64 = 0.004_65;

/// What planets form around: a single star, or a close pair together.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Disk {
    /// Initial mass, M☉.
    pub mass: f64,
    /// Main-sequence luminosity, L☉: sets the frost line and the outer extent.
    pub luminosity: f64,
    pub metallicity: f64,
    /// Where orbits are stable, AU (birth orbits).
    pub inner_au: f64,
    pub outer_au: f64,
    /// Chance a companion lets planets form at all.
    pub survival: f64,
}

/// A planet at birth.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Formed {
    pub a_au: f64,
    pub e: f64,
    pub class: BodyClass,
    pub mass_earth: f64,
}

/// A belt: inner and outer edge, AU.
pub(crate) type Belt = (f64, f64);

/// Everything orbiting one host.
#[derive(Debug, Clone, Default)]
pub(crate) struct Architecture {
    pub planets: Vec<Formed>,
    pub asteroid_belt: Option<Belt>,
    pub outer_belt: Option<Belt>,
}

/// Share of stars of mass `m` with an inner (inside the frost line) planetary system.
fn inner_system_rate(m: f64) -> f64 {
    match m {
        m if m < 0.08 => 0.7 * (m / 0.08).sqrt(),
        m if m <= 0.6 => 0.7,
        m if m <= 1.2 => 0.7 - 0.25 * (m / 0.6).log2(),
        m if m <= 2.0 => 0.45 - 0.2 * ((m - 1.2) / 0.8),
        m => (0.25 * (2.0 / m).powi(2)).max(0.02),
    }
}

/// Mean number of planets in an inner chain beyond the first, by host mass.
fn chain_extra(m: f64) -> f64 {
    if m < 0.6 {
        2.6
    } else {
        2.0
    }
}

/// Period ratio of neighbours in an inner chain: log-normal about 1.8, 0.3 wide in ln.
pub(crate) fn period_ratio(d: &mut Draw) -> f64 {
    (1.8f64.ln() + 0.3 * d.gauss()).exp().max(1.15)
}

/// Share of giant systems whose first giant migrates in the disk to a few days' orbit.
const DISK_HOT_JUPITERS: f64 = 0.02;
/// Giants are born on orbits up to this eccentric: the disk damps them, the neighbours' secular
/// pull stirs them a little.
const GIANT_BIRTH_E: f64 = 0.2;
/// Mean number of giants beyond the first in a system that makes giants.
const GIANTS_EXTRA: f64 = 0.7;
/// Period ratios of giants born in sequence: between the 3:2 and 2:1 resonances that disk
/// migration packs them into (the Solar System's giants as the Nice model starts them).
const GIANT_PACKING: (f64, f64) = (1.5, 2.0);

/// A planet in a compact chain heavier than this (Earth masses) keeps a thick envelope.
const SUB_NEPTUNE_EARTH: f64 = 3.0;

/// Share of stars of mass `m` and metallicity `feh` with giant planets.
fn giant_rate(m: f64, feh: f64) -> f64 {
    let base = 0.14 * m.min(2.5) * 10f64.powf(1.2 * feh);
    (base * if m > 2.5 { (2.5 / m).powi(2) } else { 1.0 }).min(0.9)
}

/// The planets and belts that form in `disk`.
pub(crate) fn form(d: &mut Draw, disk: &Disk) -> Architecture {
    if d.u() >= disk.survival {
        return Architecture::default();
    }
    let m = disk.mass;
    let frost = frost_line_au(disk.luminosity);
    // The ladder, from an innermost period of 1.5–12 days (outside dust sublimation).
    let p0_yr = 10f64.powf(d.range(1.5f64.log10(), 12f64.log10())) / 365.25;
    let first = ((p0_yr * p0_yr * m).cbrt()).max(0.05 * disk.luminosity.sqrt()).max(disk.inner_au * 1.05);
    let last = (40.0 * disk.luminosity.sqrt()).min(150.0 * m.sqrt()).min(disk.outer_au * 0.95);
    let mut planets = Vec::new();
    let chain = if d.u() < inner_system_rate(m) { 1 + d.poisson(chain_extra(m)) } else { 0 };
    // Around a brown dwarf, fewer and lighter planets, in step with its small disk: no ice
    // giants, and the inner chain's masses in proportion to the host's (assumed).
    let scale = (m / 0.08).min(1.0);
    // The chain's own mass, which its planets share.
    let chain_log_mass = d.range(-1.2, 1.1) + 0.5 * m.log10() + scale.log10();
    let giants = if d.u() < giant_rate(m, disk.metallicity) { (1 + d.poisson(GIANTS_EXTRA)).min(4) } else { 0 };
    let hot = giants > 0 && d.u() < DISK_HOT_JUPITERS;
    let warm = giants > 0 && !hot && d.u() < 0.2;
    let (mut giants_placed, mut small_placed) = (0, 0);
    let mut a = first;
    let mut i = 0;
    while a <= last {
        let inside = a < frost;
        let giant_here = giants_placed < giants
            && match giants_placed {
                0 if hot => i == 0,
                0 => (warm && inside && a > 0.1 && d.u() < 0.5) || !inside,
                // After a hot Jupiter the others start beyond the frost line; otherwise each
                // follows the last on the next rung, as giants migrating in the disk pack.
                _ => !inside || !hot,
            };
        i += 1;
        if giant_here {
            // Giant masses 0.1–10 M_J, more Saturns than super-Jupiters: dN/dlog m ∝ m^−0.3.
            let k = 0.3 * std::f64::consts::LN_10;
            let mass = 10f64.powf(1.5 - (1.0 - d.u() * (1.0 - (-2.0 * k).exp())).ln() / k);
            planets.push(Formed { a_au: a, e: d.range(0.0, GIANT_BIRTH_E), class: BodyClass::GasGiant, mass_earth: mass });
            giants_placed += 1;
            let packed = giants_placed < giants && !(hot && giants_placed == 1);
            a *= if packed {
                d.range(GIANT_PACKING.0, GIANT_PACKING.1).powf(2.0 / 3.0)
            } else if inside {
                period_ratio(d).powf(2.0 / 3.0)
            } else {
                d.range(1.45, 2.1)
            };
            continue;
        }
        if inside {
            if !hot && small_placed < chain {
                let mass = 10f64.powf(chain_log_mass + 0.25 * d.gauss()).clamp(0.02, 20.0);
                let class = if mass > SUB_NEPTUNE_EARTH { BodyClass::IceGiant } else { BodyClass::Terrestrial };
                planets.push(Formed { a_au: a, e: d.range(0.0, 0.08), class, mass_earth: mass });
                small_placed += 1;
            }
        } else {
            let u = d.u();
            if u < 0.2 && m >= 0.08 {
                planets.push(Formed { a_au: a, e: d.range(0.0, 0.1), class: BodyClass::IceGiant, mass_earth: d.range(8.0, 25.0) });
            } else if u < 0.35 {
                planets.push(Formed { a_au: a, e: d.range(0.0, 0.2), class: BodyClass::SmallBody, mass_earth: d.range(0.001, 0.02) });
            }
        }
        a *= if inside { period_ratio(d).powf(2.0 / 3.0) } else { d.range(1.45, 2.1) };
    }

    dynamics::settle(d, &mut planets, m, None);
    planets.retain(|p| p.a_au < disk.outer_au);

    // An asteroid belt inside the first giant, between its 4:1 and 2:1 resonances.
    let first_giant = planets.iter().find(|p| p.class == BodyClass::GasGiant && p.a_au > frost * 0.5).map(|p| p.a_au);
    let last_rock = planets.iter().filter(|p| p.class == BodyClass::Terrestrial && first_giant.is_none_or(|g| p.a_au < g)).map(|p| p.a_au).fold(0.0, f64::max);
    let asteroid_belt = first_giant.and_then(|g| dynamics::asteroid_belt(g, last_rock));
    // A Kuiper belt: outside the outermost giant, between its 3:2 and 2:1 resonances; with no
    // giant, beyond the planets and well beyond the frost line (ours is 11–18 times it).
    let outermost_giant = planets.iter().rev().find(|p| matches!(p.class, BodyClass::GasGiant | BodyClass::IceGiant) && p.a_au > frost);
    let belt = outermost_giant.map_or_else(
        || {
            let inner = planets.last().map_or(0.0, |p| p.a_au * 1.3).max(10.0 * frost);
            (inner, inner * 1.5)
        },
        |g| dynamics::kuiper_belt(g.a_au),
    );
    let outer_belt = (d.u() < 0.9 && belt.1 < disk.outer_au).then_some(belt);
    Architecture { planets, asteroid_belt, outer_belt }
}

/// What is left of `arch` around a host that is now `star`, with the host's mass now and at
/// birth (a pair's together) setting how far orbits have widened.
pub(crate) fn fate(d: &mut Draw, arch: Architecture, star: &Star, mass_now: f64, mass_birth: f64) -> Architecture {
    let widen = mass_birth / mass_now.max(1e-9);
    // Orbits widen as the star loses mass, but Hill spheres widen more: neighbours that were
    // stable may not be, and scatter (as white dwarfs polluted by infalling rock suggest).
    let mut keep_beyond = |cut: f64, a: Architecture| -> Architecture {
        let mut planets: Vec<Formed> = a.planets.into_iter().filter(|p| p.a_au > cut).map(|p| Formed { a_au: p.a_au * widen, ..p }).collect();
        dynamics::settle(d, &mut planets, mass_now, None);
        Architecture {
            planets,
            asteroid_belt: a.asteroid_belt.filter(|b| b.1 > cut).map(|b| (b.0.max(cut) * widen, b.1 * widen)),
            outer_belt: a.outer_belt.filter(|b| b.1 > cut).map(|b| (b.0.max(cut) * widen, b.1 * widen)),
        }
    };
    match star.stage {
        StellarStage::MainSequence | StellarStage::BrownDwarf => arch,
        StellarStage::Subgiant | StellarStage::Giant | StellarStage::Supergiant => keep_beyond(3.0 * star.radius_solar * SOLAR_RADIUS_AU, arch),
        StellarStage::WhiteDwarf => keep_beyond(2.0 * star.initial_mass().powf(0.8), arch),
        StellarStage::NeutronStar if d.u() < 0.005 => {
            // A second generation from the supernova's fallback disk, like PSR B1257+12's.
            let count = 1 + (d.u() * 3.0) as usize;
            let mut a = d.range(0.15, 0.25);
            let planets = (0..count)
                .map(|_| {
                    let p = Formed { a_au: a, e: d.range(0.0, 0.03), class: BodyClass::Terrestrial, mass_earth: 10f64.powf(d.range(-1.7, 0.7)) };
                    a *= d.range(1.3, 1.9);
                    p
                })
                .collect();
            Architecture { planets, asteroid_belt: None, outer_belt: None }
        }
        StellarStage::NeutronStar | StellarStage::BlackHole => Architecture::default(),
    }
}

/// Holman & Wiegert (1999): the widest stable orbit around one star of a binary (S-type), as
/// a multiple of the binary's semi-major axis; `mu` is the companion's share of the mass.
pub(crate) fn s_type_limit(mu: f64, e: f64) -> f64 {
    (0.464 - 0.380 * mu - 0.631 * e + 0.586 * mu * e + 0.150 * e * e - 0.198 * mu * e * e).max(0.02)
}

/// Holman & Wiegert (1999): the closest stable orbit around both stars (P-type), as a
/// multiple of the binary's semi-major axis.
pub(crate) fn p_type_limit(mu: f64, e: f64) -> f64 {
    1.60 + 5.10 * e - 2.22 * e * e + 4.12 * mu - 4.27 * e * mu - 5.09 * mu * mu + 4.61 * e * e * mu * mu
}

/// Chance planets form around one star of a binary `a_au` wide: none under 1 AU, a fifth at
/// 10 AU, half at 50, all beyond 200 (Moe & Kratter 2021, after Kraus et al. 2016).
pub(crate) fn binary_survival(a_au: f64) -> f64 {
    if a_au <= 1.0 {
        0.0
    } else {
        (a_au.ln() / 200f64.ln()).min(1.0).powi(2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stars::evolve;
    use geoforge_core::derived::main_sequence_luminosity;
    use geoforge_core::seed::Seed;

    fn disk(m: f64, feh: f64) -> Disk {
        Disk { mass: m, luminosity: crate::stars::birth_luminosity(m), metallicity: feh, inner_au: 0.0, outer_au: 1e4, survival: 1.0 }
    }

    fn sample(m: f64, feh: f64, n: u64) -> Vec<Architecture> {
        (0..n).map(|s| form(&mut Draw::new(Seed(s)), &disk(m, feh))).collect()
    }

    #[test]
    fn red_dwarfs_have_many_small_planets_and_few_giants() {
        let set = sample(0.35, 0.0, 4000);
        let n = set.len() as f64;
        // Small planets (rocky or sub-Neptune) within 200 days.
        let p200 = (200.0f64 / 365.25).powf(2.0 / 3.0) * 0.35f64.cbrt();
        let small = set.iter().flat_map(|a| &a.planets).filter(|p| p.a_au < p200 && p.class != BodyClass::GasGiant).count() as f64 / n;
        assert!((1.5..3.5).contains(&small), "small planets per red dwarf: {small}");
        let giants = set.iter().filter(|a| a.planets.iter().any(|p| p.class == BodyClass::GasGiant)).count() as f64 / n;
        assert!((0.02..0.08).contains(&giants), "red dwarfs with giants: {giants}");
    }

    #[test]
    fn giants_follow_metals_and_mass() {
        let rate = |m: f64, feh: f64| sample(m, feh, 4000).iter().filter(|a| a.planets.iter().any(|p| p.class == BodyClass::GasGiant)).count() as f64 / 4000.0;
        let sun = rate(1.0, 0.0);
        assert!((0.1..0.18).contains(&sun), "Sun-like stars with giants: {sun}");
        assert!(rate(1.0, 0.3) > 1.8 * sun && rate(1.0, -0.5) < 0.5 * sun, "metal-rich stars have more giants");
        assert!(rate(2.0, 0.0) > 1.6 * sun, "A stars have more giants");
        let hot = sample(1.0, 0.0, 8000).iter().filter(|a| a.planets.first().is_some_and(|p| p.class == BodyClass::GasGiant && p.a_au < 0.1)).count() as f64 / 8000.0;
        assert!((0.003..0.02).contains(&hot), "hot Jupiters {hot}");
    }

    /// Small planets inside the frost line around Sun-like stars, system by system.
    fn inner_chains(set: &[Architecture], m: f64) -> Vec<Vec<Formed>> {
        let frost = frost_line_au(main_sequence_luminosity(m));
        set.iter().map(|a| a.planets.iter().filter(|p| p.a_au < frost && p.class != BodyClass::GasGiant && p.mass_earth < 25.0).copied().collect()).collect()
    }

    /// Kepler's compact systems: neighbours' period ratios peak near 1.6 with a tail (Fabrycky
    /// et al. 2014), sizes alike within a system (Weiss et al. 2018), ~1 small planet per
    /// Sun-like star within 200 days, a third or so of them sub-Neptunes.
    #[test]
    fn inner_systems_are_kepler_like() {
        let set = sample(1.0, 0.0, 6000);
        let chains = inner_chains(&set, 1.0);
        let mut ratios: Vec<f64> = chains.iter().flat_map(|c| c.windows(2).map(|w| (w[1].a_au / w[0].a_au).powf(1.5)).collect::<Vec<_>>()).collect();
        ratios.sort_by(f64::total_cmp);
        let median = ratios[ratios.len() / 2];
        let share = |lo: f64, hi: f64| ratios.iter().filter(|r| (lo..hi).contains(*r)).count() as f64 / ratios.len() as f64;
        assert!((1.6..2.2).contains(&median), "median period ratio {median}");
        assert!(share(1.5, 2.0) > 0.3, "1.5–2: {}", share(1.5, 2.0));
        assert!(share(0.0, 1.2) < 0.03, "under 1.2: {}", share(0.0, 1.2));
        // Peas in a pod: neighbours' masses closer than two planets from different systems'.
        let spread = |pairs: &[(f64, f64)]| (pairs.iter().map(|(a, b)| (a / b).log10().powi(2)).sum::<f64>() / pairs.len() as f64).sqrt();
        let neighbours: Vec<(f64, f64)> = chains.iter().flat_map(|c| c.windows(2).map(|w| (w[0].mass_earth, w[1].mass_earth)).collect::<Vec<_>>()).collect();
        let firsts: Vec<f64> = chains.iter().filter_map(|c| c.first().map(|p| p.mass_earth)).collect();
        let strangers: Vec<(f64, f64)> = firsts.windows(2).map(|w| (w[0], w[1])).collect();
        assert!(spread(&neighbours) < 0.6 * spread(&strangers), "neighbours {} strangers {}", spread(&neighbours), spread(&strangers));
        let p200 = (200.0f64 / 365.25).powf(2.0 / 3.0);
        let per_star = chains.iter().flatten().filter(|p| p.a_au < p200).count() as f64 / 6000.0;
        assert!((0.7..2.2).contains(&per_star), "small planets within 200 days per Sun-like star: {per_star}");
        let all: Vec<&Formed> = chains.iter().flatten().collect();
        let sub = all.iter().filter(|p| p.class == BodyClass::IceGiant).count() as f64 / all.len() as f64;
        assert!((0.2..0.45).contains(&sub), "sub-Neptunes {sub}");
    }

    /// Giants from radial velocities: eccentric by Kipping's (2013) Beta law (median 0.17, 18%
    /// over 0.4); hot Jupiters lonely, warm ones often with neighbours (Huang et al. 2016).
    #[test]
    fn giants_are_eccentric_and_hot_jupiters_lonely() {
        let set = sample(1.0, 0.2, 12_000);
        let mut e: Vec<f64> = set.iter().flat_map(|a| &a.planets).filter(|p| p.class == BodyClass::GasGiant && (0.1..5.0).contains(&p.a_au)).map(|p| p.e).collect();
        e.sort_by(f64::total_cmp);
        let median = e[e.len() / 2];
        let over = e.iter().filter(|x| **x > 0.4).count() as f64 / e.len() as f64;
        assert!((0.08..0.3).contains(&median), "median giant e {median} of {}", e.len());
        assert!((0.05..0.3).contains(&over), "giants over e 0.4: {over}");
        let neighbour_share = |lo: f64, hi: f64| {
            let with: Vec<&Architecture> = set.iter().filter(|a| a.planets.iter().any(|p| p.class == BodyClass::GasGiant && (lo..hi).contains(&p.a_au))).collect();
            let crowded = with.iter().filter(|a| a.planets.iter().filter(|p| p.a_au < 1.0).count() > 1).count() as f64;
            (crowded / with.len() as f64, with.len())
        };
        let (hot, n_hot) = neighbour_share(0.0, 0.1);
        let (warm, _) = neighbour_share(0.1, 0.5);
        assert!(n_hot > 50 && hot < 0.25 && warm > 1.5 * hot, "hot Jupiters with neighbours inside 1 AU: {hot} of {n_hot}; warm: {warm}");
    }

    #[test]
    fn brown_dwarfs_have_few_small_planets() {
        let set = sample(0.04, 0.0, 3000);
        let n = set.iter().flat_map(|a| &a.planets).filter(|p| p.class != BodyClass::SmallBody).count() as f64 / 3000.0;
        let heavy = set.iter().flat_map(|a| &a.planets).filter(|p| p.mass_earth > 5.0).count() as f64 / (n * 3000.0);
        assert!((0.3..2.5).contains(&n), "planets per brown dwarf {n}");
        assert!(heavy < 0.02, "{heavy} of planets over 5 M⊕");
    }

    #[test]
    fn dead_stars_keep_little() {
        let mut d = Draw::new(Seed(3));
        let arch = Architecture { planets: [0.5, 1.0, 1.5, 3.0, 6.0].iter().map(|&a| Formed { a_au: a, e: 0.0, class: BodyClass::Terrestrial, mass_earth: 1.0 }).collect(), ..Default::default() };
        let wd = evolve(1.0, 12.0, 0.0);
        let left = fate(&mut d, arch.clone(), &wd, wd.mass_solar, 1.0);
        assert_eq!(left.planets.len(), 2, "only planets beyond ~2 AU survive the giant phases");
        assert!((left.planets[0].a_au - 3.0 / wd.mass_solar).abs() < 1e-9, "orbits widen as the star loses mass");
        // Near the tip of the red giant branch, ~100 R☉: what orbits within 1.5 AU is gone.
        let giant = evolve(1.0, 11.398, 0.0);
        assert!(giant.radius_solar > 80.0, "{giant:?}");
        assert!(fate(&mut d, arch.clone(), &giant, giant.mass_solar, 1.0).planets.len() < arch.planets.len());
        let hole = evolve(40.0, 1.0, 0.0);
        assert!(fate(&mut d, arch, &hole, hole.mass_solar, 40.0).planets.is_empty());
    }

    #[test]
    fn binaries() {
        assert!(s_type_limit(0.5, 0.0) < 0.3 && s_type_limit(0.1, 0.0) > 0.4);
        assert!((2.0..3.5).contains(&p_type_limit(0.3, 0.1)));
        assert_eq!(binary_survival(0.5), 0.0);
        assert!((0.15..0.25).contains(&binary_survival(10.0)));
        assert!(binary_survival(300.0) >= 1.0);
    }
}
