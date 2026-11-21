//! Planetary systems for any star: what forms at birth, and what is left when the star has
//! swollen into a giant or died.
//!
//! Planets form on a geometric ladder of orbits (period ratios 1.45–2.1, as in Kepler's
//! multi-planet systems) from a few days out. Rates, by the host's initial mass `M` (M☉) and
//! metallicity `[Fe/H]`:
//!
//! | Planets | Rate | After |
//! |---|---|---|
//! | an inner system at all | 70% of red dwarfs, 50% of Sun-like stars, 25% at 2 M☉, 2% at 10 | Kepler: 2.5 small planets per red dwarf within 200 days (Dressing & Charbonneau 2015) |
//! | each inner rung, given one | 70%; a third sub-Neptunes (a quarter for red dwarfs), else rock | |
//! | giant planets | 0.14·M·10^(1.2[Fe/H]) up to 2.5 M☉, falling as M⁻² above | ∝ M and 10^(1.2[Fe/H]) (Johnson et al. 2010); ~14% of Sun-like stars |
//! | hot Jupiter, given giants | 7% (about 1% of Sun-like stars) | lonely: no inner planets |
//! | each outer rung without a giant | 20% an ice giant, 15% a dwarf planet | |
//! | an outer belt, beyond the planets and 10× the frost line | 90% | debris disks are seen around ~20% of Sun-like stars, but a Kuiper belt like ours would not be |
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
        m if m <= 0.6 => 0.7,
        m if m <= 1.2 => 0.7 - 0.25 * (m / 0.6).log2(),
        m if m <= 2.0 => 0.45 - 0.2 * ((m - 1.2) / 0.8),
        m => (0.25 * (2.0 / m).powi(2)).max(0.02),
    }
}

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
    let mut rungs = Vec::new();
    let mut a = first;
    while a <= last {
        rungs.push(a);
        a *= d.range(1.45, 2.1);
    }

    let mut planets = Vec::new();
    let inner = d.u() < inner_system_rate(m);
    let has_giants = d.u() < giant_rate(m, disk.metallicity);
    let hot = has_giants && d.u() < 0.07;
    let warm = has_giants && !hot && d.u() < 0.2;
    let sub_neptunes = if m < 0.6 { 0.25 } else { 0.33 };
    let mut giants_placed = 0;
    for (i, &a) in rungs.iter().enumerate() {
        let inside = a < frost;
        let giant_here = has_giants
            && giants_placed < 4
            && if hot {
                i == 0 || (!inside && d.u() < 0.3)
            } else if giants_placed == 0 {
                (warm && inside && a > 0.1 && d.u() < 0.5) || (!inside && a < frost * 2.5) || (!inside && i + 1 == rungs.len())
            } else {
                !inside && d.u() < 0.35
            };
        if giant_here {
            // Giant masses 0.1–10 M_J, more Saturns than super-Jupiters: dN/dlog m ∝ m^−0.3.
            let k = 0.3 * std::f64::consts::LN_10;
            let mass = 10f64.powf(1.5 - (1.0 - d.u() * (1.0 - (-2.0 * k).exp())).ln() / k);
            planets.push(Formed { a_au: a, e: d.range(0.0, 0.3), class: BodyClass::GasGiant, mass_earth: mass });
            giants_placed += 1;
            continue;
        }
        if inside {
            if inner && !hot && d.u() < 0.7 {
                let (class, mass) = if d.u() < sub_neptunes { (BodyClass::IceGiant, d.range(3.0, 15.0)) } else { (BodyClass::Terrestrial, 10f64.powf(d.range(-1.5, 0.9))) };
                planets.push(Formed { a_au: a, e: d.range(0.0, 0.08), class, mass_earth: mass });
            }
        } else {
            let u = d.u();
            if u < 0.2 {
                planets.push(Formed { a_au: a, e: d.range(0.0, 0.1), class: BodyClass::IceGiant, mass_earth: d.range(8.0, 25.0) });
            } else if u < 0.35 {
                planets.push(Formed { a_au: a, e: d.range(0.0, 0.2), class: BodyClass::SmallBody, mass_earth: d.range(0.001, 0.02) });
            }
        }
    }

    // Asteroid belt between the last rocky planet and the first giant, if the gap is wide.
    let first_giant = planets.iter().find(|p| p.class == BodyClass::GasGiant && p.a_au > frost * 0.5).map(|p| p.a_au);
    let last_rock = planets.iter().filter(|p| p.class == BodyClass::Terrestrial && first_giant.is_none_or(|g| p.a_au < g)).map(|p| p.a_au).fold(0.0, f64::max);
    let asteroid_belt = first_giant.and_then(|g| {
        let inner = (last_rock * 1.35).max(frost * 0.6);
        (g * 0.72 > inner && g / last_rock.max(1e-9) > 2.2).then_some((inner, g * 0.72))
    });
    // A Kuiper belt: beyond the planets and well beyond the frost line (ours: 11–18 times it).
    let inner_edge = planets.last().map_or(0.0, |p| p.a_au * 1.3).max(10.0 * frost);
    let outer_belt = (d.u() < 0.9 && inner_edge * 1.5 < disk.outer_au).then_some((inner_edge, inner_edge * 1.5));
    Architecture { planets, asteroid_belt, outer_belt }
}

/// What is left of `arch` around a host that is now `star`, with the host's mass now and at
/// birth (a pair's together) setting how far orbits have widened.
pub(crate) fn fate(d: &mut Draw, arch: Architecture, star: &Star, mass_now: f64, mass_birth: f64) -> Architecture {
    let widen = mass_birth / mass_now.max(1e-9);
    let keep_beyond = |cut: f64, a: Architecture| -> Architecture {
        Architecture {
            planets: a.planets.into_iter().filter(|p| p.a_au > cut).map(|p| Formed { a_au: p.a_au * widen, ..p }).collect(),
            asteroid_belt: a.asteroid_belt.filter(|b| b.1 > cut).map(|b| (b.0.max(cut) * widen, b.1 * widen)),
            outer_belt: a.outer_belt.filter(|b| b.1 > cut).map(|b| (b.0.max(cut) * widen, b.1 * widen)),
        }
    };
    match star.stage {
        StellarStage::MainSequence => arch,
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
        Disk { mass: m, luminosity: main_sequence_luminosity(m), metallicity: feh, inner_au: 0.0, outer_au: 1e4, survival: 1.0 }
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

    #[test]
    fn dead_stars_keep_little() {
        let mut d = Draw::new(Seed(3));
        let arch = Architecture {
            planets: [0.5, 1.0, 1.5, 3.0, 6.0].iter().map(|&a| Formed { a_au: a, e: 0.0, class: BodyClass::Terrestrial, mass_earth: 1.0 }).collect(),
            ..Default::default()
        };
        let wd = evolve(1.0, 12.0, 0.0);
        let left = fate(&mut d, arch.clone(), &wd, wd.mass_solar, 1.0);
        assert_eq!(left.planets.len(), 2, "only planets beyond ~2 AU survive the giant phases");
        assert!((left.planets[0].a_au - 3.0 / wd.mass_solar).abs() < 1e-9, "orbits widen as the star loses mass");
        // Near the tip of the red giant branch, ~100 R☉: what orbits within 1.5 AU is gone.
        let giant = evolve(1.0, 11.48, 0.0);
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
