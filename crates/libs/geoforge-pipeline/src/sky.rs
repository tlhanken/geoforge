//! A world's sky on a date: the naked-eye stars of its neighbourhood, the galaxy's glow, the
//! other stars and planets of its system, its moons and their phases, and the eclipses to
//! come. All of it is derived from the system file and its seed; nothing is stored.
//!
//! Coordinates are the world's own right ascension and declination: its equator, and its
//! northern spring equinox (the present) at RA 0.
//!
//! Planets shine by reflected starlight with geometric albedos by class (rocky planets 0.3 —
//! Earth 0.37, Venus 0.69, Mars 0.17, Mercury 0.14; moons 0.15 — the Moon 0.12, Callisto
//! 0.20, Europa 0.67; gas giants 0.5 — Jupiter 0.54, Saturn 0.50; ice giants 0.45 — Uranus
//! 0.49, Neptune 0.44; small bodies 0.1) and a Lambert sphere's phase law.

use geoforge_core::body::{BodyClass, NodeKind, SolarSystem, AU_KM};
use geoforge_core::ephemeris::{self, Eclipse, Equatorial};
use geoforge_core::key::Key;
use geoforge_cosmos::sky::{self as stars, Neighbourhood, LY_PER_PC, NAKED_EYE};
use geoforge_cosmos::{Galaxy, GalaxyPreset};
use serde::Serialize;

/// Parsecs per AU.
const PC_PER_AU: f64 = 1.0 / 206_264.806;
/// A planet closer to the star than this in the sky is lost in its glare, degrees.
const GLARE_DEG: f64 = 15.0;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Kind {
    /// A star of the neighbourhood.
    Star,
    /// A star of the world's own system.
    Sun,
    Planet,
    Moon,
}

/// Something in the sky.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Object {
    pub name: Option<String>,
    pub kind: Kind,
    pub ra_deg: f64,
    pub dec_deg: f64,
    /// Visual magnitude.
    pub v: f64,
    /// Light-years for neighbourhood stars, AU for the system's own.
    pub distance: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub spectral_class: Option<String>,
    /// For colour: the star's, or the light a body reflects.
    pub temperature_k: f64,
    /// Angle from the brightest of the system's stars, degrees (system bodies).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub elongation_deg: Option<f64>,
    /// Share of the face lit (system bodies).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub lit: Option<f64>,
}

impl Object {
    /// Bright enough and far enough from the star to be seen at night.
    #[must_use]
    pub fn visible(&self) -> bool {
        self.v <= NAKED_EYE && self.elongation_deg.is_none_or(|e| e >= GLARE_DEG)
    }
}

/// The sky of a world on a date.
#[derive(Debug, Clone, Serialize)]
pub struct Sky {
    /// Earth days from the present.
    pub day: f64,
    pub galaxy: String,
    pub radius_kpc: f64,
    /// Neighbourhood stars, brightest first.
    pub stars: Vec<Object>,
    /// The system's stars and bodies, brightest first.
    pub bodies: Vec<Object>,
    pub eclipses: Vec<Eclipse>,
    /// The path of the system's brightest star through the year, as (RA, Dec).
    pub ecliptic: Vec<(f64, f64)>,
    #[serde(skip)]
    neighbourhood: Neighbourhood,
    #[serde(skip)]
    frame: Equatorial,
}

/// The sky of `world` on `day`, with the eclipses over the `eclipse_days` that follow.
///
/// # Errors
/// If `world` is not a body orbiting something in `system`.
pub fn sky(system: &SolarSystem, world: Key, day: f64, eclipse_days: f64) -> Result<Sky, String> {
    let frame = Equatorial::of(system, world).ok_or("the world orbits nothing")?;
    let (galaxy, radius_kpc) = galaxy_of(system);
    let neighbourhood = stars::neighbourhood(&galaxy, radius_kpc, system.root.seed.child(b"sky"));
    let stars = neighbourhood
        .stars
        .iter()
        .map(|s| {
            let (ra, dec) = frame.ra_dec(s.direction);
            Object {
                name: s.name.clone(),
                kind: Kind::Star,
                ra_deg: ra,
                dec_deg: dec,
                v: s.v,
                distance: s.distance_pc * LY_PER_PC,
                spectral_class: Some(s.star.spectral_class.clone()),
                temperature_k: s.star.temperature_k,
                elongation_deg: None,
                lit: None,
            }
        })
        .collect();
    let bodies = bodies(system, world, day, &frame);
    let ecliptic = bodies
        .iter()
        .find(|b| b.kind == Kind::Sun)
        .and_then(|sun| system.nodes.iter().find(|n| sun.name.as_deref() == Some(n.name.as_str())))
        .map(|star| {
            let year = ephemeris::period_days(system, world).unwrap_or(365.25);
            let here = |t: f64| ephemeris::position_of(system, world, t);
            (0..360)
                .map(|i| {
                    let t = day + year * f64::from(i) / 360.0;
                    frame.ra_dec(sub(ephemeris::position_of(system, star.key, t), here(t)))
                })
                .collect()
        })
        .unwrap_or_default();
    Ok(Sky { day, galaxy: galaxy.name.clone(), radius_kpc, stars, bodies, eclipses: ephemeris::eclipses(system, world, day, eclipse_days), ecliptic, neighbourhood, frame })
}

impl Sky {
    /// The galaxy's glow on an equirectangular `width`×`height` grid of RA (0–360, left to
    /// right, increasing eastward to the left as on a sky chart) and Dec (+90 at top),
    /// normalised to 1 at the brightest point.
    #[must_use]
    pub fn glow(&self, width: usize, height: usize) -> Vec<f32> {
        let mut out = vec![0.0f32; width * height];
        for (i, v) in out.iter_mut().enumerate() {
            let (ra, dec) = chart_coords(i % width, i / width, width, height);
            let (ra, dec) = (ra.to_radians(), dec.to_radians());
            let local = [dec.cos() * ra.cos(), dec.cos() * ra.sin(), dec.sin()];
            let f = &self.frame;
            let dir = [0, 1, 2].map(|k| f.x[k] * local[0] + f.y[k] * local[1] + f.z[k] * local[2]);
            *v = self.neighbourhood.glow(dir) as f32;
        }
        let top = out.iter().copied().fold(f32::MIN, f32::max).max(1e-9);
        out.iter_mut().for_each(|v| *v /= top);
        out
    }
}

/// RA and Dec at the centre of chart pixel (x, y): RA 360 at the left edge, 0 at the right.
#[must_use]
pub fn chart_coords(x: usize, y: usize, width: usize, height: usize) -> (f64, f64) {
    (360.0 * (1.0 - (x as f64 + 0.5) / width as f64), 90.0 - 180.0 * (y as f64 + 0.5) / height as f64)
}

/// The galaxy the system was drawn from and how far from its centre; the Milky Way at the
/// Sun's radius for a habitable-recipe system.
fn galaxy_of(system: &SolarSystem) -> (Galaxy, f64) {
    match &system.root.site {
        Some(site) => {
            let galaxy = GalaxyPreset::ALL.into_iter().map(GalaxyPreset::galaxy).find(|g| g.name == site.galaxy).unwrap_or_default();
            (galaxy, site.radius_kpc)
        }
        None => {
            let g = Galaxy::milky_way();
            let r = g.layout.reference_kpc;
            (g, r)
        }
    }
}

/// The system's stars, planets and the world's own moons (or, for a moon, its planet and
/// sister moons) as seen from `world` on `day`.
fn bodies(system: &SolarSystem, world: Key, day: f64, frame: &Equatorial) -> Vec<Object> {
    let at = |k: Key| ephemeris::position_of(system, k, day);
    let here = at(world);
    let parent = system.node(world).and_then(|n| n.parent);
    let star_nodes: Vec<_> = system.nodes.iter().filter_map(|n| if let NodeKind::Star(s) = &n.node { Some((n, s)) } else { None }).collect();
    // Each star's apparent magnitude from here.
    let suns: Vec<(Key, f64, [f64; 3])> = star_nodes
        .iter()
        .filter_map(|(n, s)| {
            let p = at(n.key);
            let d = norm(sub(p, here));
            stars::absolute_v(s).map(|m| (n.key, m + 5.0 * (d * PC_PER_AU / 10.0).log10(), p))
        })
        .collect();
    let brightest = suns.iter().min_by(|a, b| a.1.total_cmp(&b.1)).map(|s| s.2);
    let mut out = Vec::new();
    for (n, s) in &star_nodes {
        let Some(&(_, v, p)) = suns.iter().find(|x| x.0 == n.key) else {
            continue;
        };
        let to = sub(p, here);
        let (ra, dec) = frame.ra_dec(to);
        let elongation = brightest.filter(|b| *b != p).map(|b| angle(to, sub(b, here)));
        out.push(Object {
            name: Some(n.name.clone()),
            kind: Kind::Sun,
            ra_deg: ra,
            dec_deg: dec,
            v,
            distance: norm(to),
            spectral_class: Some(s.spectral_class.clone()),
            temperature_k: s.temperature_k,
            elongation_deg: elongation,
            lit: None,
        });
    }
    for n in &system.nodes {
        let NodeKind::Body(b) = &n.node else { continue };
        if n.key == world {
            continue;
        }
        let up = n.parent.and_then(|p| system.node(p));
        let is_planet = up.is_some_and(|p| matches!(p.node, NodeKind::Star(_) | NodeKind::Barycenter));
        let mine = n.parent == Some(world) || (parent.is_some() && (Some(n.key) == parent || n.parent == parent));
        if !is_planet && !mine {
            continue;
        }
        let p = at(n.key);
        // Lit by the star that looks brightest from the body.
        let Some(&(_, _, star)) = suns.iter().min_by(|a, c| {
            let f = |s: &(Key, f64, [f64; 3])| s.1 + 5.0 * (norm(sub(s.2, p)) / norm(sub(s.2, here))).log10();
            f(a).total_cmp(&f(c))
        }) else {
            continue;
        };
        let sun_v = suns.iter().find(|s| s.2 == star).map_or(0.0, |s| s.1);
        let (to, from_star) = (sub(p, here), sub(p, star));
        let (d, r) = (norm(to), norm(from_star));
        let phase = angle(sub(here, p), sub(star, p)).to_radians();
        let lambert = (phase.sin() + (std::f64::consts::PI - phase) * phase.cos()) / std::f64::consts::PI;
        let radius_au = b.radius_km / AU_KM;
        let kind = if is_planet && n.parent != Some(world) { Kind::Planet } else { Kind::Moon };
        let flux = albedo(b.class, kind) * (radius_au / d).powi(2) * lambert * (norm(sub(star, here)) / r).powi(2);
        let (ra, dec) = frame.ra_dec(to);
        let elongation = brightest.map(|s| angle(to, sub(s, here)));
        let lit = ephemeris::lit_fraction(p, here, star);
        let temperature_k = star_nodes.iter().find(|(m, _)| at(m.key) == star).map_or(5772.0, |(_, s)| s.temperature_k);
        out.push(Object {
            name: Some(n.name.clone()),
            kind,
            ra_deg: ra,
            dec_deg: dec,
            v: sun_v - 2.5 * flux.max(1e-300).log10(),
            distance: d,
            spectral_class: None,
            temperature_k,
            elongation_deg: elongation,
            lit: Some(lit),
        });
    }
    out.sort_by(|a, b| a.v.total_cmp(&b.v));
    out
}

/// Geometric albedo by class.
fn albedo(class: BodyClass, kind: Kind) -> f64 {
    match class {
        BodyClass::Terrestrial if kind == Kind::Moon => 0.15,
        BodyClass::Terrestrial => 0.3,
        BodyClass::GasGiant => 0.5,
        BodyClass::IceGiant => 0.45,
        BodyClass::SmallBody => 0.1,
    }
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn norm(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// Angle between two vectors, degrees.
fn angle(a: [f64; 3], b: [f64; 3]) -> f64 {
    let c = (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]) / (norm(a) * norm(b)).max(1e-300);
    c.clamp(-1.0, 1.0).acos().to_degrees()
}

#[cfg(test)]
mod tests {
    use super::*;
    use geoforge_core::Archetype;

    /// The Sun, Earth, Moon and Jupiter at opposition: Jupiter near −2.7, the full Moon near
    /// −12.7, the Sun −26.7, all in the right places.
    #[test]
    fn the_solar_system_seen_from_earth() {
        let earth = Archetype::Earth.context(Key(1), "Earth").exterior;
        let (mut moon, mut jupiter) = (earth, earth);
        (moon.mass_kg, moon.radius_km) = (7.342e22, 1737.4);
        (jupiter.class, jupiter.mass_kg, jupiter.radius_km) = (BodyClass::GasGiant, 1.898e27, 69_911.0);
        let body = |b| serde_json::to_value(NodeKind::Body(b)).unwrap();
        // At the equinox Earth is at −x from the Sun; Jupiter beyond it, the Moon behind Earth.
        let sys: SolarSystem = serde_json::from_value(serde_json::json!({
            "root": { "key": 1, "seed": 1, "name": "Sol", "epoch": { "system_age_gyr": 4.57 } },
            "nodes": [
                { "key": 2, "name": "Sun", "node": { "kind": "star", "mass_solar": 1.0, "luminosity_solar": 1.0, "radius_solar": 1.0, "temperature_k": 5772.0, "age_gyr": 4.57, "metallicity_dex": 0.0, "spectral_class": "G2V" } },
                { "key": 3, "name": "Earth", "parent": 2, "node": body(earth), "elements": { "a_au": 1.0, "e": 0.0, "i_deg": 0.0, "lan_deg": 0.0, "argp_deg": 0.0, "m0_deg": 180.0 } },
                { "key": 4, "name": "Moon", "parent": 3, "node": body(moon), "elements": { "a_au": 0.002_569_6, "e": 0.0, "i_deg": 5.0, "lan_deg": 0.0, "argp_deg": 0.0, "m0_deg": 180.0 } },
                { "key": 5, "name": "Jupiter", "parent": 2, "node": body(jupiter), "elements": { "a_au": 5.2, "e": 0.0, "i_deg": 0.0, "lan_deg": 0.0, "argp_deg": 0.0, "m0_deg": 180.0 } }
            ]
        }))
        .unwrap();
        let s = sky(&sys, Key(3), 0.0, 400.0).unwrap();
        let get = |n: &str| s.bodies.iter().find(|b| b.name.as_deref() == Some(n)).unwrap();
        let (sun, moon, jup) = (get("Sun"), get("Moon"), get("Jupiter"));
        assert!((sun.v + 26.75).abs() < 0.1, "Sun {}", sun.v);
        assert!(sun.ra_deg.min(360.0 - sun.ra_deg) < 1e-6 && sun.dec_deg.abs() < 1e-6);
        assert!((jup.v + 2.7).abs() < 0.3 && jup.elongation_deg.unwrap() > 179.0 && jup.visible(), "{jup:?}");
        assert!((moon.v + 12.7).abs() < 0.3 && moon.lit.unwrap() > 0.99, "{moon:?}");
        // A full Moon this close to the node: a lunar eclipse now, and some to come.
        assert!(!s.eclipses.is_empty());
        assert!(s.stars.len() > 4000 && s.ecliptic.len() == 360);
        let glow = s.glow(36, 18);
        assert!(glow.iter().copied().fold(0.0, f32::max) == 1.0);
    }
}
