//! Energy balance by latitude (Budyko–Sellers) on the body's own radiative properties.
//!
//! Absorbed starlight by latitude balances outgoing longwave from the gray atmosphere
//! ([`Greenhouse`]) and poleward heat transport `C·(T − T̄)`, carried by the ocean where
//! there is one and by the air in proportion to its mass. Frost brightens cold ground:
//! water ice and snow on wet worlds, seasonal CO₂ frost under thin CO₂ air (Mars). Heat transport and the damping of seasons grow with the mass of
//! the air column, so thick atmospheres (Venus, Titan) are nearly uniform and thin ones (Mars)
//! have strong gradients and large seasons.

use geoforge_core::derived::{daily_insolation, insolation_by_latitude};
use geoforge_core::radiation::{balance_temperature_c, open_albedo, Greenhouse};
use geoforge_core::body::Species;
use geoforge_core::BodyContext;
use std::f64::consts::PI;

/// Earth's radiative damping, W/m²/K: converts the user's greenhouse offset to W/m².
const EARTH_OLR_SLOPE: f64 = 2.09;
/// How much an ocean damps the seasons of the whole air mass, interiors included; fitted so
/// Earth's continental interiors get 0.10 °C per W/m² of peak-to-peak insolation.
const OCEAN_MODERATION: f64 = 0.42;
/// Share of Earth's poleward heat transport carried by the ocean, which does not depend on
/// the air's pressure; the rest scales with the air column.
const OCEAN_TRANSPORT_SHARE: f64 = 0.4;
/// Extra albedo of frost-covered ground.
const FROST_ALBEDO: f64 = 0.32;
/// Annual mean temperature at which half the ground is frosted, °C.
const WATER_FROST_C: f64 = -8.0;
const CO2_FROST_C: f64 = -75.0;

#[derive(Debug, Clone, Copy)]
pub struct Radiation {
    greenhouse: Greenhouse,
    open_albedo: f64,
    /// Annual mean temperature where frost covers half the ground, if any can form.
    frost_c: Option<f64>,
    /// Poleward heat transport coefficient C, W/m²/K.
    transport: f64,
    /// Outgoing longwave removed by the user's greenhouse offset, W/m².
    olr_offset: f64,
    flux: f64,
    obliquity: f64,
    eccentricity: f64,
    /// Seasonal range, °C per W/m² of peak-to-peak insolation.
    seasonal_gain: f64,
    /// Radiative balance without ice or transport: the solver's starting point.
    start_c: f64,
}

impl Radiation {
    /// `heat_transport` is Earth's coefficient (scaled here with air mass);
    /// `greenhouse_offset_c` adds (or removes) warming on top of the atmosphere's own.
    #[must_use]
    pub fn for_body(body: &BodyContext, heat_transport: f64, greenhouse_offset_c: f64) -> Self {
        let air = &body.interior;
        let greenhouse = Greenhouse::of(body);
        let pressure = air.surface_pressure_bar();
        let start_c = balance_temperature_c(body);
        let wet = air.has_surface_liquid_water();
        let open_albedo = open_albedo(air);
        let damping = greenhouse.olr_slope(start_c).max(0.05);
        let impedance = seasonal_impedance(damping, column_heat_j_m2_k(pressure, body.gravity_m_s2()), body.irradiance.year_days);
        let moderation = if wet { OCEAN_MODERATION } else { 1.0 };
        Self {
            greenhouse,
            open_albedo,
            frost_c: if wet {
                Some(WATER_FROST_C)
            } else if air.partial_pressure_bar(Species::CarbonDioxide) > 0.5 * pressure && pressure < 0.1 {
                Some(CO2_FROST_C)
            } else {
                None
            },
            transport: heat_transport * (if wet { OCEAN_TRANSPORT_SHARE } else { 0.0 } + (1.0 - OCEAN_TRANSPORT_SHARE) * (pressure / 1.013).sqrt()),
            olr_offset: EARTH_OLR_SLOPE * greenhouse_offset_c,
            flux: body.irradiance.flux_w_m2,
            obliquity: body.irradiance.obliquity_deg,
            eccentricity: body.irradiance.eccentricity,
            seasonal_gain: (1.0 - open_albedo) * moderation / impedance,
            start_c,
        }
    }

    /// Bond albedo at surface temperature `t` (°C).
    #[must_use]
    pub fn albedo(&self, t: f64) -> f64 {
        let frost = self.frost_c.map_or(0.0, |f| 0.5 * FROST_ALBEDO * (1.0 - ((t - f) / 10.0).tanh()));
        self.open_albedo + frost
    }

    fn olr(&self, t: f64) -> f64 {
        self.greenhouse.olr_w_m2(t) - self.olr_offset
    }

    /// Zonal sea-level mean temperatures (°C) at the given latitudes, and their global mean.
    #[must_use]
    pub fn zonal_temperature(&self, lats_deg: &[f64]) -> (Vec<f64>, f64) {
        let s: Vec<f64> = lats_deg.iter().map(|&l| insolation_by_latitude(self.flux, l, self.obliquity, self.eccentricity).annual_w_m2).collect();
        self.balance(lats_deg, &s, &vec![0.0; lats_deg.len()])
    }

    /// Temperatures (°C) of a tidally locked world by "substellar latitude": 90° under the
    /// star, 0 on the terminator, −90° opposite (Koll & Abbot's tidally locked coordinates),
    /// and the global mean. Sunlight is the star's flux times the sine of that latitude by
    /// day, nothing by night; the same heat transport spreads it. On a wet world the rising
    /// air under the star builds a thick cloud deck that reflects much of the light and keeps
    /// the day side from running away (Yang, Cowan & Abbot 2013): up to 0.35 more albedo
    /// straight under the star, over water warmer than 30 °C (none below freezing).
    #[must_use]
    pub fn substellar_temperature(&self, lats_deg: &[f64]) -> (Vec<f64>, f64) {
        let s: Vec<f64> = lats_deg.iter().map(|&l| self.flux * l.to_radians().sin().max(0.0)).collect();
        let deck = if self.frost_c == Some(WATER_FROST_C) { SUBSTELLAR_CLOUDS } else { 0.0 };
        let clouds: Vec<f64> = lats_deg.iter().map(|&l| deck * l.to_radians().sin().max(0.0)).collect();
        // Day–night contrasts are what the whole circulation works on, with little rotation
        // to hold them apart: heat moves more readily than from equator to pole.
        let locked = Self { transport: self.transport * LOCKED_TRANSPORT, ..*self };
        locked.balance(lats_deg, &s, &clouds)
    }

    /// Energy balance by band: sunlight `s`, extra albedo `clouds`.
    fn balance(&self, lats_deg: &[f64], s: &[f64], clouds: &[f64]) -> (Vec<f64>, f64) {
        let w: Vec<f64> = lats_deg.iter().map(|l| l.to_radians().cos()).collect();
        let wsum: f64 = w.iter().sum();
        let c = self.transport;
        let mut t = vec![self.start_c; lats_deg.len()];
        let mut mean = self.start_c;
        for _ in 0..400 {
            mean = t.iter().zip(&w).map(|(a, b)| a * b).sum::<f64>() / wsum;
            // Newton step on absorbed − emitted − exported, damped by half.
            let mut worst: f64 = 0.0;
            for ((ti, &si), &ci) in t.iter_mut().zip(s).zip(clouds) {
                // Deep convective cloud needs warm water under it: none below 0 °C, full by 30.
                let deck = ci * (*ti / 30.0).clamp(0.0, 1.0);
                let net = si * (1.0 - self.albedo(*ti) - deck) - self.olr(*ti) - c * (*ti - mean);
                let step = 0.5 * net / (self.greenhouse.olr_slope(*ti).max(0.05) + c);
                *ti += step;
                worst = worst.max(step.abs());
            }
            if worst < 1e-4 {
                break;
            }
        }
        (t, mean)
    }

    /// Signed seasonal range at a latitude, °C: positive when the warm season is northern
    /// summer. Land in the continental interior gets all of it, the sea about an eighth.
    #[must_use]
    pub fn seasonal_range(&self, lat_deg: f64) -> f64 {
        seasonal_swing(self.flux, lat_deg, self.obliquity, self.eccentricity) * self.seasonal_gain
    }
}

/// Heat transport of a locked world relative to a rotating one.
const LOCKED_TRANSPORT: f64 = 1.6;

/// Extra albedo of the cloud deck straight under the star of a locked wet world.
const SUBSTELLAR_CLOUDS: f64 = 0.35;

/// Heat capacity of the air column, J/m²/K (cₚ ≈ 1000 J/kg/K).
fn column_heat_j_m2_k(pressure_bar: f64, gravity_m_s2: f64) -> f64 {
    pressure_bar * 1e5 / gravity_m_s2 * 1000.0
}

/// |B + iωC|: how strongly a column resists an annual forcing, W/m²/K.
fn seasonal_impedance(damping: f64, column_heat: f64, year_days: f64) -> f64 {
    let omega = 2.0 * PI / (year_days.max(1.0) * 86_400.0);
    damping.hypot(omega * column_heat)
}

/// Peak-to-peak annual harmonic of daily insolation at a latitude (W/m²), positive when the
/// peak is in northern summer.
#[must_use]
pub fn seasonal_swing(flux: f64, lat_deg: f64, obliquity: f64, e: f64) -> f64 {
    let lat = lat_deg.to_radians().clamp(-PI / 2.0 + 1e-6, PI / 2.0 - 1e-6);
    let eps = obliquity.to_radians();
    let n = 48;
    let (mut re, mut im) = (0.0, 0.0);
    for i in 0..n {
        let ls = 2.0 * PI * (f64::from(i) + 0.5) / f64::from(n);
        let decl = (eps.sin() * ls.sin()).asin();
        // Distance factor: perihelion near northern winter, as on Earth.
        let s = flux * (1.0 + 2.0 * e * (ls + PI / 2.0 + 0.2).cos());
        let q = daily_insolation(s, lat, decl);
        re += q * ls.sin();
        im += q * ls.cos();
    }
    let amp = 2.0 * (re * re + im * im).sqrt() / f64::from(n);
    2.0 * amp * re.signum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use geoforge_core::{Archetype, Key};

    fn zonal(a: Archetype) -> (Vec<f64>, f64, Radiation) {
        let body = a.context(Key(1), a.slug());
        let r = Radiation::for_body(&body, 3.2, 0.0);
        let lats: Vec<f64> = (0..180).map(|i| 89.5 - f64::from(i)).collect();
        let (t, mean) = r.zonal_temperature(&lats);
        (t, mean, r)
    }

    #[test]
    /// Earth locked to the Sun: ~30 °C under the star and ~−40 °C opposite (locked aquaplanet
    /// GCMs: 300–310 K and 230–250 K), cooler overall than turning (Yang et al. 2013).
    fn earth_locked() {
        let (_, mean_z, r) = zonal(Archetype::Earth);
        let lats: Vec<f64> = (0..180).map(|i| -89.5 + f64::from(i)).collect();
        let (t, mean) = r.substellar_temperature(&lats);
        assert!((20.0..45.0).contains(&t[179]) && t[0] < -30.0, "substellar {:.0} °C, antistellar {:.0} °C", t[179], t[0]);
        assert!((mean_z - 45.0..mean_z).contains(&mean), "locked worlds run cooler (clouds by day, ice by night): {mean} vs {mean_z}");
    }

    #[test]
    fn earth_zonal_means() {
        let (t, mean, _) = zonal(Archetype::Earth);
        let at = |l: f64| t[(89.5 - l) as usize];
        // Earth (sea-level, zonal, annual): ~14–15 °C global, ~26 °C equator, ~0 °C at 60°,
        // ~-20 °C near the poles.
        assert!((13.5..17.0).contains(&mean), "global {mean}");
        assert!((25.0..28.5).contains(&at(0.5)), "equator {}", at(0.5));
        assert!((8.0..13.0).contains(&at(45.5)), "45N {}", at(45.5));
        assert!((-5.0..3.0).contains(&at(60.5)), "60N {}", at(60.5));
        assert!((-22.0..-10.0).contains(&at(89.5)), "pole {}", at(89.5));
    }

    #[test]
    fn other_worlds() {
        let span = |t: &[f64]| t.iter().copied().fold(f64::NEG_INFINITY, f64::max) - t.iter().copied().fold(f64::INFINITY, f64::min);
        // Venus: hot and nearly uniform, with no seasons.
        let (t, mean, r) = zonal(Archetype::Venus);
        assert!((430.0..490.0).contains(&mean), "Venus {mean}");
        assert!(span(&t) < 15.0, "Venus span {}", span(&t));
        assert!(r.seasonal_range(45.0).abs() < 1.0);
        // Mars: cold, steep gradient, big seasons.
        let (t, mean, r) = zonal(Archetype::Mars);
        assert!((-75.0..-50.0).contains(&mean), "Mars {mean}");
        assert!(span(&t) > 50.0, "Mars span {}", span(&t));
        assert!(r.seasonal_range(50.0) > 30.0, "Mars range {}", r.seasonal_range(50.0));
        // Titan: cold, thick air, small gradient.
        let (t, mean, _) = zonal(Archetype::Titan);
        assert!((-185.0..-173.0).contains(&mean), "Titan {mean}");
        assert!(span(&t) < 8.0, "Titan span {}", span(&t));
    }

    #[test]
    fn seasons_follow_latitude() {
        let s0 = seasonal_swing(1361.0, 0.0, 23.44, 0.0).abs();
        let s60 = seasonal_swing(1361.0, 60.0, 23.44, 0.0);
        let s60s = seasonal_swing(1361.0, -60.0, 23.44, 0.0);
        assert!(s0 < 60.0, "{s0}");
        assert!(s60 > 380.0 && s60s < -380.0, "{s60} {s60s}");
        // No obliquity, no seasons.
        assert!(seasonal_swing(1361.0, 45.0, 0.0, 0.0).abs() < 1.0);
        // Earth's own gain is the reference.
        let (_, _, r) = zonal(Archetype::Earth);
        assert!((r.seasonal_range(60.0) / seasonal_swing(1361.0, 60.0, 23.44, 0.0167) - 0.10).abs() < 0.005);
    }
}
