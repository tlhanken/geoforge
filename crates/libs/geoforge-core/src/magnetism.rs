//! A body's magnetic field: whether a dynamo runs, how strong it is, how far it holds off the
//! star's wind, and where its auroras fall. Derived from the body, never stored.
//!
//! - **Dynamo.** A rocky world's iron core convects only while the mantle draws its heat away.
//!   Moving plates cool the mantle fast and keep it going for billions of years (Earth's has
//!   run at least 3.5 Gyr); under one lid the core soon stalls (Mars's stopped ~0.5 Gyr after
//!   it formed, Venus has none). Lifetime: 10 Gyr·√M with plates, 1.5 Gyr·√M under a
//!   stagnant lid, 1 Gyr·√M under an episodic one (M in Earth masses: bigger worlds stay hot
//!   longer). Giants always have one: metallic hydrogen in gas giants, ionic water in ice
//!   giants.
//! - **Strength.** Earth's equatorial surface field, 30 µT, scaled as M^(1/3) for rocky worlds
//!   (a bigger core, a surface farther out); Jupiter's 420 µT scaled by mass for gas giants;
//!   ice giants ~20 µT. A slow spinner (a day over 200 h, as tidally locked planets of red
//!   dwarfs) makes a weaker multipolar field (×0.3 at the surface), as ice giants do.
//! - **Tilt.** Up to 12° from the spin axis (Earth ~10°, Jupiter ~10°); ice giants 40–60°
//!   (Uranus 59°, Neptune 47°).
//! - **Wind.** The star's wind ram pressure, 2 nPa at Earth, falls as 1/d² (d from the flux
//!   and a main-sequence luminosity L ≈ (T/5772 K)^7) and is stronger around young stars
//!   (× 4.6 Gyr / age, within 0.5–20).
//! - **Magnetopause** (Chapman–Ferraro): where the field's pressure balances the wind,
//!   R = 1.26·(B²/(2µ₀P))^(1/6) planet radii (the 1.26 for the compressed field doubling at
//!   the boundary): Earth ~10.
//! - **Auroras**: on the field lines that reach the magnetopause, sin²θ = 1/R from the
//!   magnetic pole: Earth's oval at ~71° magnetic latitude (observed 65–75°).
//!
//! Exceptions the rules miss: Mercury keeps a weak field (0.3 µT) under its lid; Ganymede
//! has a dynamo inside its ice shell; Saturn's field is weak (21 µT) for its mass; and moons
//! inside a giant's magnetosphere feel the giant's plasma, not the star's wind.

use crate::body::{BodyClass, BodyContext, ShellKind, TectonicRegime, EARTH_MASS_KG};
use crate::seed::{unit, Seed};

const EARTH_FIELD_UT: f64 = 30.0;
const JUPITER_FIELD_UT: f64 = 420.0;
const JUPITER_MASS_KG: f64 = 1.898e27;
const ICE_GIANT_FIELD_UT: f64 = 20.0;
/// Wind ram pressure at Earth, nPa.
const EARTH_WIND_NPA: f64 = 2.0;
/// A day longer than this makes a multipolar field.
const MULTIPOLAR_DAY_H: f64 = 200.0;
const MU0: f64 = 4e-7 * std::f64::consts::PI;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MagneticField {
    /// A dynamo runs today.
    pub dynamo: bool,
    /// How long a rocky world's dynamo lasts, Gyr after it formed; `None` for giants (always)
    /// and bodies with no core to run one.
    pub lifetime_gyr: Option<f64>,
    /// Equatorial surface field, µT; 0 without a dynamo.
    pub surface_ut: f64,
    /// Dipole-dominated (Earth, Jupiter) rather than multipolar (Uranus, slow spinners).
    pub dipolar: bool,
    /// Angle between the magnetic and spin axes, degrees.
    pub tilt_deg: f64,
    /// The star's wind ram pressure at the body, nPa.
    pub wind_npa: f64,
    /// Magnetopause distance on the star's side, planet radii; `None` without a field.
    pub magnetopause_radii: Option<f64>,
    /// Magnetic latitude of the auroral ovals, degrees; `None` where the wind reaches the air
    /// (no field, or a magnetopause at the surface).
    pub aurora_lat_deg: Option<f64>,
}

impl MagneticField {
    #[must_use]
    pub fn of(body: &BodyContext) -> Self {
        let ext = &body.exterior;
        let age = body.epoch.system_age_gyr;
        let mass = ext.mass_kg / EARTH_MASS_KG;
        let slow = ext.rotation_h.abs() > MULTIPOLAR_DAY_H;
        // Keys repeat between systems; the mass tells bodies apart.
        let tilt = |most: f64| most * unit(Seed(body.key.0 ^ ext.mass_kg.to_bits()).child(b"magnetic-tilt"));
        // A slow spinner's field is multipolar, weaker at the surface.
        let spin = if slow { 0.3 } else { 1.0 };
        let (lifetime_gyr, dynamo, surface_ut, dipolar, tilt_deg) = match ext.class {
            BodyClass::GasGiant => (None, true, JUPITER_FIELD_UT * ext.mass_kg / JUPITER_MASS_KG * spin, !slow, tilt(12.0)),
            BodyClass::IceGiant => (None, true, ICE_GIANT_FIELD_UT, false, 40.0 + tilt(20.0)),
            _ => {
                let has_core = body.interior.shells.iter().any(|s| s.kind == ShellKind::Core);
                let per_root_mass = match body.interior.regime {
                    TectonicRegime::Plates => Some(10.0),
                    TectonicRegime::StagnantLid => Some(1.5),
                    TectonicRegime::Episodic => Some(1.0),
                    TectonicRegime::IceShell | TectonicRegime::Inactive => None,
                };
                let life = per_root_mass.filter(|_| has_core).map(|b| b * mass.sqrt());
                let on = life.is_some_and(|l| l > age);
                (life, on, if on { EARTH_FIELD_UT * mass.cbrt() * spin } else { 0.0 }, !slow, tilt(12.0))
            }
        };
        let wind_npa = wind_pressure_npa(body.irradiance.flux_w_m2, body.irradiance.star_temperature_k, age);
        let magnetopause_radii = dynamo.then(|| magnetopause_radii(surface_ut, wind_npa));
        let aurora_lat_deg = magnetopause_radii.filter(|&r| r > 1.0).map(|r| 90.0 - (1.0 / r).sqrt().asin().to_degrees());
        Self { dynamo, lifetime_gyr, surface_ut, dipolar, tilt_deg, wind_npa, magnetopause_radii, aurora_lat_deg }
    }
}

/// The star's wind ram pressure, nPa, at a body getting `flux` W/m² from a star of
/// temperature `star_k` and age `age_gyr`.
#[must_use]
pub fn wind_pressure_npa(flux: f64, star_k: f64, age_gyr: f64) -> f64 {
    let luminosity = (star_k / 5772.0).powi(7);
    let inverse_d2 = (flux / 1361.0) / luminosity;
    let youth = (4.6 / age_gyr.max(1e-3)).clamp(0.5, 20.0);
    EARTH_WIND_NPA * inverse_d2 * youth
}

/// Magnetopause distance, planet radii, for an equatorial surface field `b_ut` against a wind
/// of `p_npa`.
#[must_use]
pub fn magnetopause_radii(b_ut: f64, p_npa: f64) -> f64 {
    let b = b_ut * 1e-6;
    1.26 * (b * b / (2.0 * MU0 * p_npa * 1e-9)).powf(1.0 / 6.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Archetype, Key};

    fn field(a: Archetype) -> MagneticField {
        MagneticField::of(&a.context(Key(1), "x"))
    }

    #[test]
    fn earth() {
        let f = field(Archetype::Earth);
        assert!(f.dynamo && f.dipolar && (f.surface_ut - 30.0).abs() < 1.0 && f.tilt_deg <= 12.0, "{f:?}");
        let r = f.magnetopause_radii.unwrap();
        assert!((9.0..11.0).contains(&r), "magnetopause {r}");
        let oval = f.aurora_lat_deg.unwrap();
        assert!((65.0..75.0).contains(&oval), "auroras at {oval}");
        assert!((f.wind_npa - 2.0).abs() < 0.3, "{}", f.wind_npa);
    }

    #[test]
    fn dead_and_live_dynamos() {
        for a in [Archetype::Mars, Archetype::Venus, Archetype::Moon, Archetype::Mercury] {
            let f = field(a);
            assert!(!f.dynamo && f.surface_ut == 0.0 && f.aurora_lat_deg.is_none(), "{a:?}: {f:?}");
        }
        assert!(field(Archetype::Mars).lifetime_gyr.unwrap() < 1.0, "Mars's dynamo stopped early");
        assert!(field(Archetype::Europa).lifetime_gyr.is_none());
        let jupiter = field(Archetype::GasGiant);
        assert!(jupiter.dynamo && (400.0..440.0).contains(&jupiter.surface_ut), "{jupiter:?}");
        // A pure dipole gives ~40 Jupiter radii; Io's plasma inflates the real one to 50–100.
        assert!((30.0..60.0).contains(&jupiter.magnetopause_radii.unwrap()), "{jupiter:?}");
        let ice = field(Archetype::IceGiant);
        assert!(ice.dynamo && !ice.dipolar && ice.tilt_deg >= 40.0 && (5.0..30.0).contains(&ice.surface_ut), "{ice:?}");
    }

    #[test]
    fn bodies_and_stars_shape_the_field() {
        // A super-Earth keeps its dynamo longer and stronger; a tidally locked world of a red
        // dwarf, slow and close, has a weak multipolar field squeezed by a strong wind.
        let mut big = Archetype::Earth.context(Key(1), "big");
        big.exterior.mass_kg *= 5.0;
        let big = MagneticField::of(&big);
        assert!(big.lifetime_gyr.unwrap() > 20.0 && big.surface_ut > 45.0);
        let mut locked = Archetype::Earth.context(Key(1), "locked");
        locked.irradiance.star_temperature_k = 3300.0;
        locked.irradiance.year_days = 20.0;
        locked.exterior.rotation_h = 480.0;
        let locked = MagneticField::of(&locked);
        let earth = field(Archetype::Earth);
        assert!(!locked.dipolar && locked.surface_ut < 10.0 && locked.wind_npa > 30.0 * earth.wind_npa, "{locked:?}");
        assert!(locked.magnetopause_radii.unwrap() < 0.5 * earth.magnetopause_radii.unwrap());
        // A young star blows harder.
        let mut young = Archetype::Earth.context(Key(1), "young");
        young.epoch.system_age_gyr = 0.5;
        assert!(MagneticField::of(&young).wind_npa > 5.0 * earth.wind_npa);
    }
}
