//! Bands layer: the cloud tops of a world without a surface, a gas or ice giant.
//!
//! **Jets.** Zonal jets are spaced by the Rhines scale, L = π·√(2U/β) with β = 2Ω·cos 45°/a,
//! so a hemisphere has about (πa/2)/L of them. With U ≈ 60 m/s for gas giants this gives
//! Jupiter six and Saturn five; ice giants, with winds of ~250 m/s, get one (Uranus and
//! Neptune). Gas giants superrotate at the equator (Jupiter +100 m/s, Saturn +400), and
//! between the jets the wind runs weakly backwards; ice giants' equators run retrograde
//! (Neptune −400 m/s) with a prograde jet near 70°.
//!
//! **Bands.** Zones (bright, rising ammonia cloud) lie where the shear is anticyclonic,
//! belts (dark, where deeper cloud shows) where it is cyclonic; the equator is a bright zone.
//! Ice giants show the same pattern faintly, under clear methane-blue air.
//!
//! **Weather.** Band edges are wavy (festoons and chevrons where the jets are unstable), and
//! long-lived anticyclones sit in the zones: white ovals, a few per jet. Half of giants also
//! have one great storm as wide as its zone (Jupiter's Great Red Spot; Saturn has none that
//! lasts; on ice giants, a dark spot like Neptune's, which come and go).
//!
//! Products: `zonal_wind` (m/s, eastward) and `cloud_bands` (cloud-top brightness, 0 for the
//! darkest belt to 1 for the brightest zone or oval).

use geoforge_core::body::BodyClass;
use geoforge_core::grid::Cell;
use geoforge_core::product::{Generator, LayerId, Precondition, ProductData, ProductId, Products, RunContext};
use geoforge_core::seed::{unit, Seed};
use geoforge_core::shapes::Field;
use geoforge_core::{BodyContext, Error};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct BandsParams {
    /// Typical jet speed, m/s. `None`: 60 for gas giants, 250 for ice giants.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub jet_speed_m_s: Option<f64>,
    /// Long-lived ovals per jet, on average.
    pub storms_per_jet: f64,
    /// How wavy the band edges are (0 = smooth).
    pub turbulence: f64,
    /// Chance of one great storm, as wide as its zone. Jupiter has kept its Great Red Spot
    /// since at least 1830; Saturn has none that lasts (a Great White Spot erupts about once a
    /// Saturn year and fades); Neptune's dark spots come and go within years. So: a coin toss.
    pub great_storm_chance: f64,
}

impl Default for BandsParams {
    fn default() -> Self {
        Self { jet_speed_m_s: None, storms_per_jet: 1.5, turbulence: 1.0, great_storm_chance: 0.5 }
    }
}

/// Jets per hemisphere of a giant of radius `radius_km` turning in `rotation_h` with jets of
/// `speed` m/s, from the Rhines scale.
#[must_use]
pub fn giant_jets(rotation_h: f64, radius_km: f64, speed: f64) -> u32 {
    let a = radius_km * 1000.0;
    let omega = 2.0 * PI / (rotation_h.abs().max(0.5) * 3600.0);
    let beta = 2.0 * omega * 45f64.to_radians().cos() / a;
    let rhines = PI * (2.0 * speed / beta).sqrt();
    ((PI * a / 2.0) / rhines).round().clamp(1.0, 20.0) as u32
}

/// The zonal wind profile of a giant: eastward wind (m/s) by latitude.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Profile {
    pub jets: u32,
    pub speed: f64,
    /// Ice giants: retrograde equator.
    pub ice: bool,
}

impl Profile {
    #[must_use]
    pub fn for_body(body: &BodyContext, p: &BandsParams) -> Self {
        let ice = body.exterior.class == BodyClass::IceGiant;
        let speed = p.jet_speed_m_s.unwrap_or(if ice { 250.0 } else { 60.0 });
        Self { jets: giant_jets(body.exterior.rotation_h, body.exterior.radius_km, speed), speed, ice }
    }

    /// Half-width of the equatorial jet, degrees.
    fn equatorial(&self) -> f64 {
        180.0 / f64::from(2 * self.jets + 1)
    }

    /// Eastward wind at `lat` (degrees), m/s.
    #[must_use]
    pub fn wind(&self, lat: f64) -> f64 {
        let l = lat.abs().min(90.0);
        let u = self.speed;
        if self.ice && self.jets <= 1 {
            // Neptune: retrograde to 50°, a prograde jet near 70°.
            return if l < 50.0 { -1.6 * u * (0.5 * PI * l / 50.0).cos() } else { u * (PI * (l - 50.0) / 40.0).sin() };
        }
        let weq = self.equatorial();
        let sign = if self.ice { -1.0 } else { 1.0 };
        if l < weq {
            // A flat-topped equatorial jet.
            let x = ((l / weq - 0.5) / 0.5).max(0.0);
            return sign * 1.6 * u * (0.5 * PI * x).cos().powi(2);
        }
        // Alternating jets out to the pole: weak retrograde troughs, strong prograde peaks.
        let pair = (90.0 - weq) / f64::from(self.jets);
        let s = (2.0 * PI * (l - weq) / pair).sin();
        u * 0.7 * l.to_radians().cos().sqrt() * (-s * (0.75 - 0.25 * s))
    }

    /// Cloud-top brightness of the bands at `lat`: bright where the shear is anticyclonic.
    #[must_use]
    pub fn brightness(&self, lat: f64) -> f64 {
        let h = 0.25;
        let shear = (self.wind(lat + h) - self.wind(lat - h)) / (2.0 * h) * lat.signum();
        let scale = self.speed * f64::from(self.jets.max(1)) / 40.0;
        // Bands are sharpest near the equator and fade toward the mottled polar regions.
        let contrast = if self.ice { 0.12 } else { 0.42 } * lat.to_radians().cos().powf(1.5);
        let eq = 0.25 * (-(lat / (0.4 * self.equatorial())).powi(2)).exp();
        (0.5 + contrast * (shear / scale).tanh() + eq).clamp(0.0, 1.0)
    }
}

/// A long-lived oval storm.
#[derive(Debug, Clone, Copy)]
struct Oval {
    lat: f64,
    lon: f64,
    half_lat: f64,
    half_lon: f64,
    /// Brightness inside: white ovals on gas giants, dark spots on ice giants.
    value: f64,
}

/// Deterministic draws from one seed.
struct Draws {
    seed: Seed,
    n: u64,
}

impl Draws {
    fn u(&mut self) -> f64 {
        self.n += 1;
        unit(self.seed.child_u64(b"d", self.n))
    }
    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.u()
    }
}

/// Ovals in the zones: centred where the brightness peaks (the anticyclonic zones).
fn storms(d: &mut Draws, profile: &Profile, p: &BandsParams) -> Vec<Oval> {
    let zones: Vec<(f64, f64)> = (-89..=89)
        .map(f64::from)
        .filter(|&l| l.abs() > 3.0)
        .filter(|&l| profile.brightness(l) > profile.brightness(l - 1.0) && profile.brightness(l) >= profile.brightness(l + 1.0) && profile.brightness(l) > 0.6)
        .map(|l| {
            let width = (1..30).map(f64::from).take_while(|k| profile.brightness(l + k) > 0.5 && profile.brightness(l - k) > 0.5).count() as f64;
            (l, width.max(1.0))
        })
        .collect();
    let mut out = Vec::new();
    // Perhaps one great storm in the widest subtropical zone: a dark red oval on a gas giant
    // (the Great Red Spot), a dark spot on an ice giant.
    let has_great = d.u() < p.great_storm_chance;
    let great = zones.iter().filter(|z| (12.0..45.0).contains(&z.0.abs())).max_by(|a, b| a.1.total_cmp(&b.1)).copied();
    // Ice giants' faint bands may have no subtropical zone: the spot sits anywhere at 15–35°.
    let great = great.or_else(|| Some((d.range(15.0, 35.0) * if d.u() < 0.5 { -1.0 } else { 1.0 }, 5.0)));
    if let (true, Some((lat, width))) = (has_great, great) {
        let half_lat = (0.9 * width).max(4.0);
        out.push(Oval { lat, lon: d.range(-180.0, 180.0), half_lat, half_lon: half_lat * 1.9 / lat.to_radians().cos(), value: if profile.ice { 0.08 } else { 0.22 } });
    }
    // White ovals spread over the zones; on ice giants, bright methane streaks instead.
    let expected = p.storms_per_jet * f64::from(2 * profile.jets);
    let share = expected / zones.len().max(1) as f64;
    for &(lat, width) in zones.iter().filter(|z| !profile.ice || z.0.abs() > 15.0) {
        for _ in 0..(share + d.u()).floor() as usize {
            let half_lat = (width * d.range(0.3, 0.55)).max(1.5);
            let stretch = if profile.ice { d.range(4.0, 8.0) } else { d.range(1.4, 2.2) };
            out.push(Oval {
                lat: lat + d.range(-0.2, 0.2) * width,
                lon: d.range(-180.0, 180.0),
                half_lat: if profile.ice { half_lat * 0.4 } else { half_lat },
                half_lon: (half_lat * stretch / lat.to_radians().cos().max(0.2)).min(30.0),
                value: 0.97,
            });
        }
    }
    out
}

pub struct GiantBands;

impl Generator for GiantBands {
    type Params = BandsParams;
    const NAME: &'static str = "giant_bands";
    const VERSION: u32 = 1;
    const LAYER: LayerId = LayerId::Bands;

    fn needs() -> &'static [ProductId] {
        &[]
    }
    fn applies(body: &BodyContext) -> Result<(), Precondition> {
        crate::preconditions::deep_atmosphere(body)
    }

    fn run(ctx: &RunContext<'_>, inputs: &Products, p: &BandsParams) -> Result<Products, Error> {
        let mut out = Products::new();
        if inputs.contains(ProductId::ZonalWind) && inputs.contains(ProductId::CloudBands) {
            return Ok(out);
        }
        let profile = Profile::for_body(ctx.body, p);
        let mut d = Draws { seed: ctx.seed(ProductId::CloudBands), n: 0 };
        let ovals = storms(&mut d, &profile, p);
        // Waves along the band edges: a few zonal wavenumbers, stronger at mid-latitudes.
        let waves: Vec<(f64, f64, f64)> = (0..6).map(|k| (f64::from(3 + 4 * k) + d.range(0.0, 3.0), d.range(0.0, 2.0 * PI), d.range(0.0, 2.0 * PI))).collect();
        // Band edges wander by up to a fifth of a band (a couple of degrees on Jupiter).
        let amp = p.turbulence * (0.1 * (90.0 - profile.equatorial()) / f64::from(profile.jets)).min(2.5);
        let grid = ctx.grid;
        let cells: Vec<(f32, f32)> = (0..grid.len())
            .into_par_iter()
            .map(|i| {
                let ll = grid.center(i as Cell).to_lat_lon();
                let (lat, lon) = (ll.lat_deg, ll.lon_deg);
                let lam = lon.to_radians();
                let wobble: f64 = waves.iter().enumerate().map(|(k, &(m, ph, q))| (m * lam + ph + q * lat.to_radians()).sin() / (1.0 + k as f64)).sum();
                let lat_w = (lat + amp * wobble * (2.0 * lat.to_radians()).sin().abs()).clamp(-90.0, 90.0);
                let u = profile.wind(lat_w);
                let mut b = profile.brightness(lat_w);
                for o in &ovals {
                    let dlon = ((lon - o.lon + 540.0).rem_euclid(360.0)) - 180.0;
                    let r2 = ((lat - o.lat) / o.half_lat).powi(2) + (dlon / o.half_lon).powi(2);
                    if r2 < 1.0 {
                        let core = (1.0 - r2).sqrt();
                        b += (o.value - b) * (0.4 + 0.6 * core);
                    }
                }
                (u as f32, b.clamp(0.0, 1.0) as f32)
            })
            .collect();
        let spec = grid.spec();
        if !inputs.contains(ProductId::ZonalWind) {
            out.insert(ProductData::ZonalWind(Field::from_values(spec, cells.iter().map(|c| c.0).collect())));
        }
        if !inputs.contains(ProductId::CloudBands) {
            out.insert(ProductData::CloudBands(Field::from_values(spec, cells.iter().map(|c| c.1).collect())));
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use geoforge_core::{Archetype, Key};

    #[test]
    fn jets_of_the_solar_system_giants() {
        assert_eq!(giant_jets(9.93, 69_911.0, 60.0), 6, "Jupiter");
        assert_eq!(giant_jets(10.7, 58_232.0, 60.0), 5, "Saturn");
        assert_eq!(giant_jets(17.2, 25_559.0, 250.0), 1, "Uranus");
        assert_eq!(giant_jets(16.1, 24_764.0, 250.0), 1, "Neptune");
    }

    #[test]
    fn profiles() {
        let jupiter = Profile { jets: 6, speed: 60.0, ice: false };
        assert!(jupiter.wind(0.0) > 80.0, "equatorial superrotation");
        let signs: Vec<bool> = (15..80).map(|l| jupiter.wind(f64::from(l)) > 0.0).collect();
        assert!(signs.windows(2).filter(|w| w[0] != w[1]).count() >= 9, "alternating jets");
        assert!(jupiter.brightness(0.0) > 0.7, "bright equatorial zone");
        let bands: Vec<bool> = (5..80).map(|l| jupiter.brightness(f64::from(l)) > 0.5).collect();
        assert!(bands.windows(2).filter(|w| w[0] != w[1]).count() >= 8, "zones and belts");
        let neptune = Profile { jets: 1, speed: 250.0, ice: true };
        assert!(neptune.wind(0.0) < -350.0 && neptune.wind(70.0) > 200.0 && neptune.wind(90.0).abs() < 1.0);
    }

    #[test]
    fn great_storms_are_a_coin_toss() {
        let profile = Profile { jets: 6, speed: 60.0, ice: false };
        let p = BandsParams::default();
        let great = (0..200).filter(|&k| storms(&mut Draws { seed: Seed(k), n: 0 }, &profile, &p).iter().any(|o| o.value < 0.5)).count();
        assert!((60..140).contains(&great), "{great} of 200 with a great storm");
    }

    #[test]
    fn applies_to_giants_only() {
        assert!(GiantBands::applies(&Archetype::GasGiant.context(Key(1), "j")).is_ok());
        assert!(GiantBands::applies(&Archetype::IceGiant.context(Key(1), "n")).is_ok());
        assert!(GiantBands::applies(&Archetype::Earth.context(Key(1), "e")).is_err());
    }
}
