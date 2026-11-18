//! Derived values: pure functions of primary data, computed on read and never stored as truth.

use crate::body::{AU_KM, SOLAR_CONSTANT_W_M2, SOLAR_MASS_KG};
use std::f64::consts::PI;

pub const G: f64 = 6.674_30e-11;
pub const STEFAN_BOLTZMANN: f64 = 5.670_374e-8;

/// Orbital period in days from Kepler's third law.
#[must_use]
pub fn orbital_period_days(a_au: f64, central_mass_kg: f64) -> f64 {
    let a_m = a_au * AU_KM * 1000.0;
    2.0 * PI * (a_m.powi(3) / (G * central_mass_kg)).sqrt() / 86_400.0
}

#[must_use]
pub fn surface_gravity_m_s2(mass_kg: f64, radius_km: f64) -> f64 {
    G * mass_kg / (radius_km * 1000.0).powi(2)
}

/// Stellar flux at distance `a_au` from a star of `luminosity_solar`, W/m².
#[must_use]
pub fn flux_w_m2(luminosity_solar: f64, a_au: f64) -> f64 {
    SOLAR_CONSTANT_W_M2 * luminosity_solar / (a_au * a_au)
}

/// Equilibrium temperature, K, for a rapidly rotating body with Bond albedo `albedo`.
#[must_use]
pub fn equilibrium_temperature_k(flux_w_m2: f64, albedo: f64) -> f64 {
    (flux_w_m2 * (1.0 - albedo) / (4.0 * STEFAN_BOLTZMANN)).powf(0.25)
}

/// Conservative habitable zone (inner, outer) in AU, after Kopparapu et al. (2013), Sun-like
/// effective flux limits 1.107 and 0.356 S☉ scaled by luminosity.
#[must_use]
pub fn habitable_zone_au(luminosity_solar: f64) -> (f64, f64) {
    ((luminosity_solar / 1.107).sqrt(), (luminosity_solar / 0.356).sqrt())
}

/// Main-sequence luminosity (L☉) from mass (M☉), piecewise mass–luminosity relation.
#[must_use]
pub fn main_sequence_luminosity(mass_solar: f64) -> f64 {
    let m = mass_solar;
    if m < 0.43 {
        0.23 * m.powf(2.3)
    } else if m < 2.0 {
        m.powi(4)
    } else if m < 55.0 {
        1.4 * m.powf(3.5)
    } else {
        32_000.0 * m
    }
}

/// Main-sequence radius (R☉) from mass (M☉).
#[must_use]
pub fn main_sequence_radius(mass_solar: f64) -> f64 {
    if mass_solar < 1.0 { mass_solar.powf(0.8) } else { mass_solar.powf(0.57) }
}

/// Effective temperature (K) from luminosity and radius in solar units.
#[must_use]
pub fn effective_temperature_k(luminosity_solar: f64, radius_solar: f64) -> f64 {
    5772.0 * (luminosity_solar / (radius_solar * radius_solar)).powf(0.25)
}

/// Frost line, AU (water ice condenses beyond it).
#[must_use]
pub fn frost_line_au(luminosity_solar: f64) -> f64 {
    2.7 * luminosity_solar.sqrt()
}

#[must_use]
pub fn solar_masses(kg: f64) -> f64 {
    kg / SOLAR_MASS_KG
}

/// Daily-mean top-of-atmosphere insolation (W/m²) at latitude `lat` for solar declination
/// `decl` (radians) and instantaneous flux `s`.
#[must_use]
pub fn daily_insolation(s: f64, lat: f64, decl: f64) -> f64 {
    let x = -lat.tan() * decl.tan();
    let h0 = if x >= 1.0 { 0.0 } else if x <= -1.0 { PI } else { x.acos() };
    (s / PI) * (h0 * lat.sin() * decl.sin() + lat.cos() * decl.cos() * h0.sin())
}

/// Annual-mean, January-like and July-like insolation at latitude `lat_deg`.
///
/// Integrates daily insolation over solar longitude for a body with obliquity
/// `obliquity_deg`, mean flux `mean_flux` and eccentricity `e` (to first order).
#[must_use]
pub fn insolation_by_latitude(mean_flux: f64, lat_deg: f64, obliquity_deg: f64, e: f64) -> Insolation {
    let lat = lat_deg.to_radians().clamp(-PI / 2.0 + 1e-9, PI / 2.0 - 1e-9);
    let eps = obliquity_deg.to_radians();
    let steps = 360;
    let mut sum = 0.0;
    let s = mean_flux / (1.0 - e * e).sqrt();
    for i in 0..steps {
        let ls = 2.0 * PI * (i as f64 + 0.5) / steps as f64;
        let decl = (eps.sin() * ls.sin()).asin();
        sum += daily_insolation(s, lat, decl);
    }
    let june = daily_insolation(mean_flux, lat, eps);
    let december = daily_insolation(mean_flux, lat, -eps);
    Insolation { annual_w_m2: sum / steps as f64, june_w_m2: june, december_w_m2: december }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Insolation {
    pub annual_w_m2: f64,
    pub june_w_m2: f64,
    pub december_w_m2: f64,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::body::EARTH_MASS_KG;

    #[test]
    fn earth_numbers() {
        let p = orbital_period_days(1.0, SOLAR_MASS_KG);
        assert!((p - 365.25).abs() < 0.5, "{p}");
        let g = surface_gravity_m_s2(EARTH_MASS_KG, 6371.0);
        assert!((g - 9.82).abs() < 0.02, "{g}");
        let t = equilibrium_temperature_k(1361.0, 0.3);
        assert!((t - 255.0).abs() < 1.5, "{t}");
        // Annual-mean insolation: ~417 W/m² at the equator, ~173 at the poles for Earth.
        let eq = insolation_by_latitude(1361.0, 0.0, 23.44, 0.0167).annual_w_m2;
        let pole = insolation_by_latitude(1361.0, 90.0, 23.44, 0.0167).annual_w_m2;
        assert!((eq - 417.0).abs() < 6.0, "{eq}");
        assert!((pole - 173.0).abs() < 6.0, "{pole}");
        // Global mean is S/4.
        let mut acc = 0.0;
        let mut w = 0.0;
        for i in 0..180 {
            let lat = -89.5 + i as f64;
            let c = lat.to_radians().cos();
            acc += insolation_by_latitude(1361.0, lat, 23.44, 0.0).annual_w_m2 * c;
            w += c;
        }
        assert!((acc / w - 340.25).abs() < 1.5, "{}", acc / w);
    }

    #[test]
    fn habitable_zone_sun() {
        let (i, o) = habitable_zone_au(1.0);
        assert!(i < 1.0 && o > 1.5);
    }
}
