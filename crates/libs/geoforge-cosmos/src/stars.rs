//! Stellar evolution, simply: what a star of a given initial mass and age is now.
//!
//! | Initial mass | Life | End |
//! |---|---|---|
//! | 0.08–8 M☉ | main sequence → subgiant → red giant | white dwarf |
//! | 8–20 M☉ | main sequence → supergiant | neutron star |
//! | over 20 M☉ | main sequence → supergiant | black hole |
//!
//! The main sequence lasts `10 Gyr · M/L` (the Sun: 10 Gyr; a 2 M☉ star: 1.25 Gyr; a red
//! dwarf: longer than the universe is old), and the star brightens through it, from 0.72 of its
//! mid-life luminosity to twice it (the Sun: 0.72 at birth, 1 now, 2 at the end). The giant
//! phases follow: the red giant branch takes 7% as long again, quickening toward its tip, then
//! helium burns in the core for a fifth of the main-sequence life, at most 120 Myr (the red
//! clump); massive stars spend 10% as long again as supergiants. Remnants then cool: white dwarfs by
//! Mestel's law, with masses from the initial–final mass relation of Cummings et al. (2018).
//!
//! Brown dwarfs (under 0.08 M☉) never burn hydrogen: they shine by their birth heat and fade,
//! by Burrows & Liebert's (1993) fits, L ≈ 4×10⁻⁵ L☉·(t/Gyr)^−1.3·(M/0.05 M☉)^2.64 and
//! T ≈ 1550 K·(t/Gyr)^−0.32·(M/0.05 M☉)^0.83, which put them about a Jupiter radius across,
//! through the spectral types M, L (from 2400 K), T (1300 K) and Y (550 K). This is a star
//! alone; a pair close enough to touch lives through [`crate::binary`].

use geoforge_core::body::{Star, StellarStage};
use geoforge_core::derived::{effective_temperature_k, main_sequence_luminosity, main_sequence_radius};

/// Below this mass (M☉) an object never burns hydrogen: a brown dwarf.
pub const HYDROGEN_BURNING_LIMIT: f64 = 0.08;
/// Stars born below this mass (M☉) end as white dwarfs.
pub const WHITE_DWARF_LIMIT: f64 = 8.0;
/// Stars born above this mass (M☉) end as black holes; between the limits, neutron stars.
pub const BLACK_HOLE_LIMIT: f64 = 20.0;
/// End of the subgiant phase and of the red giant branch, in main-sequence lifetimes.
const SUBGIANT_END: f64 = 1.07;
const RED_GIANT_END: f64 = 1.14;
/// How the red giant branch quickens: time per factor of luminosity ∝ L^−CLIMB.
const CLIMB: f64 = 0.85;
/// End of the supergiant phase of a massive star, in main-sequence lifetimes.
const SUPERGIANT_END: f64 = 1.10;
/// Neutron star mass (M☉) and radius (km).
const NEUTRON_STAR_MASS: f64 = 1.4;
const NEUTRON_STAR_RADIUS_KM: f64 = 12.0;
const SOLAR_RADIUS_KM: f64 = 695_700.0;

/// Main-sequence lifetime, Gyr.
#[must_use]
pub fn main_sequence_lifetime_gyr(initial_mass: f64) -> f64 {
    let m = initial_mass;
    if m <= 10.0 {
        10.0 * m / main_sequence_luminosity(m)
    } else {
        // Massive stars burn a larger share of their hydrogen: lifetimes level off near 3 Myr.
        10.0 * 10.0 / main_sequence_luminosity(10.0) * (m / 10.0).powf(-1.1)
    }
}

/// When the star leaves the stage it is in now and enters the next, Gyr after its birth.
#[must_use]
pub fn end_of_life_gyr(initial_mass: f64) -> f64 {
    let t_ms = main_sequence_lifetime_gyr(initial_mass);
    if initial_mass < WHITE_DWARF_LIMIT {
        RED_GIANT_END * t_ms + helium_burning_gyr(t_ms)
    } else {
        SUPERGIANT_END * t_ms
    }
}

/// How long a giant burns helium in its core, Gyr: a fifth of its main-sequence life, but no
/// more than 120 Myr for stars whose cores ignite in a flash at a fixed mass (Schaller et
/// al. 1992: 0.1 Gyr for the Sun, 0.25 Gyr for 2 M☉, 0.08 Gyr for 3 M☉).
fn helium_burning_gyr(t_ms: f64) -> f64 {
    (0.2 * t_ms).min(0.12)
}

/// Luminosity through the main sequence relative to the mass–luminosity relation's value,
/// at `f` of the way through it (0 at birth, 1 at the end).
fn brightening(f: f64) -> f64 {
    0.72 + 0.037 * f + 1.243 * f * f
}

/// The star born with `initial_mass` (M☉) and metallicity `metallicity_dex`, `age_gyr` later.
#[must_use]
pub fn evolve(initial_mass: f64, age_gyr: f64, metallicity_dex: f64) -> Star {
    let m0 = initial_mass;
    if m0 < HYDROGEN_BURNING_LIMIT {
        let (l, t) = (brown_dwarf_luminosity(m0, age_gyr), brown_dwarf_temperature_k(m0, age_gyr));
        let r = l.sqrt() * (5772.0 / t).powi(2);
        let stage = StellarStage::BrownDwarf;
        return Star { mass_solar: m0, luminosity_solar: l, radius_solar: r, temperature_k: t, age_gyr, metallicity_dex, spectral_class: classify(stage, t), stage, initial_mass_solar: Some(m0) };
    }
    let t_ms = main_sequence_lifetime_gyr(m0);
    let f = age_gyr / t_ms;
    let l_ms = main_sequence_luminosity(m0);
    let r_ms = main_sequence_radius(m0);
    let make = |stage, mass: f64, l: f64, r: f64| {
        let t = if l > 0.0 { effective_temperature_k(l, r) } else { 0.0 };
        Star { mass_solar: mass, luminosity_solar: l, radius_solar: r, temperature_k: t, age_gyr, metallicity_dex, spectral_class: classify(stage, t), stage, initial_mass_solar: Some(m0) }
    };
    if f < 1.0 {
        let l = l_ms * brightening(f);
        return make(StellarStage::MainSequence, m0, l, r_ms * (l / l_ms).powf(0.45));
    }
    // At the end of the main sequence.
    let l_end = l_ms * brightening(1.0);
    let t_end = effective_temperature_k(l_end, r_ms * brightening(1.0).powf(0.45));
    let with_temperature = |l: f64, t: f64| l.sqrt() * (5772.0 / t).powi(2);
    if m0 < WHITE_DWARF_LIMIT {
        if f < SUBGIANT_END {
            // Across the Hertzsprung gap: cooling toward 5100 K and brightening.
            let x = (f - 1.0) / (SUBGIANT_END - 1.0);
            let l = l_end * (1.0 + 2.0 * x);
            let t = t_end + (t_end.min(5100.0) - t_end) * x;
            return make(StellarStage::Subgiant, m0, l, with_temperature(l, t));
        }
        if f < RED_GIANT_END {
            // Up the red giant branch, quickening as it brightens: the time spent per factor
            // of luminosity falls as L^−0.85 (the core grows at a rate ∝ L, and L ∝ M_core^7),
            // so most red giants are near its foot. Heavier stars start brighter and climb
            // less far; the tip is 3600 K.
            let y = (f - SUBGIANT_END) / (RED_GIANT_END - SUBGIANT_END);
            let base = (3.0 * l_end).max(10.0);
            let span = (2.2 - 0.4 * m0.log2().max(0.0)).clamp(1.0, 2.2);
            let top = 10f64.powf(span);
            let l = base * (1.0 - y * (1.0 - top.powf(-CLIMB))).powf(-1.0 / CLIMB);
            let t = 5000.0 - 1400.0 * (l / base).log10() / span;
            return make(StellarStage::Giant, m0 * (1.0 - 0.2 * y), l, with_temperature(l, t));
        }
        let helium_end = RED_GIANT_END * t_ms + helium_burning_gyr(t_ms);
        if age_gyr < helium_end {
            // Burning helium in the core: the red clump, 50 L☉ and 4800 K for Sun-like stars;
            // heavier stars brighter and warmer (Capella's pair: 79 and 73 L☉, 4970 and 5730 K).
            let l = (1.5 * l_end).max(50.0);
            let t = (4800.0 + 400.0 * m0.log2().max(0.0)).min(5600.0);
            return make(StellarStage::Giant, 0.8 * m0, l, with_temperature(l, t));
        }
        let mass = white_dwarf_mass(m0);
        let cooling = (age_gyr - helium_end).max(1e-3);
        return make(StellarStage::WhiteDwarf, mass, white_dwarf_luminosity(mass, cooling), white_dwarf_radius(mass));
    }
    if f < SUPERGIANT_END {
        // Across the top of the diagram at nearly constant luminosity, blue to red.
        let x = (f - 1.0) / (SUPERGIANT_END - 1.0);
        let l = 1.5 * l_ms;
        let t = t_end * (3500.0 / t_end).powf(x);
        return make(StellarStage::Supergiant, m0 * (1.0 - 0.3 * x), l, with_temperature(l, t));
    }
    let cooling = (age_gyr - SUPERGIANT_END * t_ms).max(1e-6);
    if m0 < BLACK_HOLE_LIMIT {
        // Thermal glow only; a pulsar's spin-down power goes into its wind and beams.
        let t = (5e5 * (cooling / 1e-3).powf(-0.25)).clamp(3e4, 2e6);
        let r = NEUTRON_STAR_RADIUS_KM / SOLAR_RADIUS_KM;
        let l = r * r * (t / 5772.0).powi(4);
        return make(StellarStage::NeutronStar, NEUTRON_STAR_MASS, l, r);
    }
    let mass = (0.4 * m0 - 3.0).clamp(5.0, 20.0);
    // Schwarzschild radius: 2.95 km per solar mass.
    make(StellarStage::BlackHole, mass, 0.0, 2.95 * mass / SOLAR_RADIUS_KM)
}

/// White dwarf mass (M☉) from the initial mass: Cummings et al. (2018), whose slopes are
/// 0.08, 0.19 and 0.11 per M☉ below 2.85 M☉, between 2.85 and 3.6, and above.
#[must_use]
pub fn white_dwarf_mass(initial_mass: f64) -> f64 {
    let m = initial_mass;
    let mf = if m < 2.85 {
        0.080 * m + 0.489
    } else if m < 3.60 {
        0.187 * m + 0.184
    } else {
        0.107 * m + 0.471
    };
    mf.clamp(0.5, 1.38)
}

/// Mestel cooling: luminosity (L☉) of a carbon–oxygen white dwarf of `mass` (M☉) that has
/// cooled for `cooling_gyr` (0.6 M☉: 8×10⁻⁴ L☉ after 1 Gyr, 3×10⁻⁵ after 10).
#[must_use]
pub fn white_dwarf_luminosity(mass: f64, cooling_gyr: f64) -> f64 {
    (8.8e-3 * mass.powf(5.0 / 7.0) / cooling_gyr.max(1e-4)).powf(1.4).min(100.0)
}

/// Brown dwarfs younger than this (Gyr) are taken at this age: the cooling fits start once
/// deuterium burning is over.
const BROWN_DWARF_YOUNGEST_GYR: f64 = 0.01;

/// Luminosity (L☉) of a brown dwarf of `mass` (M☉) at `age_gyr` (Burrows & Liebert 1993).
#[must_use]
pub fn brown_dwarf_luminosity(mass: f64, age_gyr: f64) -> f64 {
    4e-5 * age_gyr.max(BROWN_DWARF_YOUNGEST_GYR).powf(-1.3) * (mass / 0.05).powf(2.64)
}

/// Effective temperature (K) of a brown dwarf of `mass` (M☉) at `age_gyr` (Burrows & Liebert
/// 1993).
#[must_use]
pub fn brown_dwarf_temperature_k(mass: f64, age_gyr: f64) -> f64 {
    1550.0 * age_gyr.max(BROWN_DWARF_YOUNGEST_GYR).powf(-0.32) * (mass / 0.05).powf(0.83)
}

/// The light a disk around a star born with `mass` (M☉) feels while its planets form: the
/// main-sequence luminosity, or for a brown dwarf its glow at 10 Myr, when disks clear.
#[must_use]
pub fn birth_luminosity(mass: f64) -> f64 {
    if mass < HYDROGEN_BURNING_LIMIT {
        brown_dwarf_luminosity(mass, BROWN_DWARF_YOUNGEST_GYR)
    } else {
        main_sequence_luminosity(mass)
    }
}

/// White dwarf radius (R☉) from Nauenberg's mass–radius relation (0.6 M☉: 0.0125 R☉).
#[must_use]
pub fn white_dwarf_radius(mass: f64) -> f64 {
    0.0127 * mass.powf(-1.0 / 3.0) * (1.0 - (mass / 1.44).powf(4.0 / 3.0)).max(0.01).sqrt()
}

/// Spectral class: Morgan–Keenan letter, subclass and luminosity class for living stars
/// ("G2V", "K1III"); "DA7" (50 400 K / T) or "DC" for white dwarfs; "NS" and "BH".
#[must_use]
pub fn classify(stage: StellarStage, teff: f64) -> String {
    let suffix = match stage {
        StellarStage::MainSequence => "V",
        StellarStage::Subgiant => "IV",
        StellarStage::Giant => "III",
        StellarStage::Supergiant => "I",
        StellarStage::WhiteDwarf if teff >= 5000.0 => return format!("DA{:.0}", (50_400.0 / teff).clamp(1.0, 12.0)),
        StellarStage::WhiteDwarf => return "DC".into(),
        StellarStage::NeutronStar => return "NS".into(),
        StellarStage::BlackHole => return "BH".into(),
        StellarStage::BrownDwarf => {
            let (letter, sub) = match teff {
                t if t >= 2_400.0 => mk_letter(t),
                t if t >= 1_300.0 => ("L", ((2_400.0 - t) / 110.0) as u32),
                t if t >= 550.0 => ("T", ((1_300.0 - t) / 75.0) as u32),
                t => ("Y", (((550.0 - t) / 100.0) as u32).min(4)),
            };
            return format!("{letter}{}", sub.min(9));
        }
    };
    let (letter, sub) = mk_letter(teff);
    format!("{letter}{sub}{suffix}")
}

/// Main-sequence spectral class from temperature ("G2V").
#[must_use]
pub fn spectral_class(teff: f64) -> String {
    classify(StellarStage::MainSequence, teff)
}

fn mk_letter(teff: f64) -> (&'static str, u32) {
    let (letter, lo, hi) = match teff {
        t if t >= 30_000.0 => ("O", 30_000.0, 50_000.0),
        t if t >= 10_000.0 => ("B", 10_000.0, 30_000.0),
        t if t >= 7_500.0 => ("A", 7_500.0, 10_000.0),
        t if t >= 6_000.0 => ("F", 6_000.0, 7_500.0),
        t if t >= 5_200.0 => ("G", 5_200.0, 6_000.0),
        t if t >= 3_700.0 => ("K", 3_700.0, 5_200.0),
        _ => ("M", 2_400.0, 3_700.0),
    };
    // O2 is the earliest type in use.
    (letter, ((hi - teff) / (hi - lo) * 10.0).clamp(if letter == "O" { 2.0 } else { 0.0 }, 9.0) as u32)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_sun_through_its_life() {
        let now = evolve(1.0, 4.6, 0.0);
        assert_eq!(now.stage, StellarStage::MainSequence);
        assert!((now.luminosity_solar - 1.0).abs() < 0.02 && (now.radius_solar - 1.0).abs() < 0.02, "{now:?}");
        assert!(now.spectral_class.starts_with('G') && now.spectral_class.ends_with('V'));
        let young = evolve(1.0, 0.0, 0.0);
        assert!((young.luminosity_solar - 0.72).abs() < 0.01, "faint young Sun");
        assert_eq!(evolve(1.0, 10.3, 0.0).stage, StellarStage::Subgiant);
        let giant = evolve(1.0, 11.2, 0.0);
        assert_eq!(giant.stage, StellarStage::Giant);
        assert!(giant.radius_solar > 10.0 && giant.spectral_class.ends_with("III"), "{giant:?}");
        let wd = evolve(1.0, 12.0, 0.0);
        assert_eq!(wd.stage, StellarStage::WhiteDwarf);
        assert!((0.55..0.6).contains(&wd.mass_solar) && wd.radius_solar < 0.015 && wd.luminosity_solar < 0.01, "{wd:?}");
    }

    #[test]
    fn lifetimes_and_remnants() {
        // Main-sequence lifetimes against stellar models: 2 M☉ ~1.2 Gyr, 5 ~0.1, 20 ~0.01.
        for (m, t) in [(2.0, 1.2), (5.0, 0.1), (20.0, 0.0095), (60.0, 0.0037)] {
            let got = main_sequence_lifetime_gyr(m);
            assert!((got / t - 1.0).abs() < 0.35, "{m} M☉: {got} Gyr vs {t}");
        }
        assert!(main_sequence_lifetime_gyr(0.3) > 100.0, "red dwarfs outlive the universe");
        assert_eq!(evolve(12.0, 1.0, 0.0).stage, StellarStage::NeutronStar);
        assert_eq!(evolve(40.0, 1.0, 0.0).stage, StellarStage::BlackHole);
        assert_eq!(evolve(40.0, 0.001, 0.0).spectral_class.chars().next(), Some('O'));
        assert_eq!(evolve(5.0, 0.01, 0.0).spectral_class.chars().next(), Some('B'));
        // Sirius B: born ~5 M☉, now ~1 M☉.
        assert!((white_dwarf_mass(5.0) - 1.0).abs() < 0.05);
        // Supergiants run blue to red at high luminosity.
        let blue = evolve(15.0, main_sequence_lifetime_gyr(15.0) * 1.005, 0.0);
        let red = evolve(15.0, main_sequence_lifetime_gyr(15.0) * 1.095, 0.0);
        assert_eq!(blue.stage, StellarStage::Supergiant);
        assert!(blue.temperature_k > 10_000.0 && red.temperature_k < 4_000.0 && red.radius_solar > 300.0, "{blue:?} {red:?}");
    }

    #[test]
    fn habitable_classes_are_unchanged() {
        assert_eq!(spectral_class(5772.0), "G2V");
        assert_eq!(spectral_class(4500.0), "K4V");
    }
}
