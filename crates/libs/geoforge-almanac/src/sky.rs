//! A world's sky on a date: the naked-eye stars of its neighbourhood, the galaxy's glow, the
//! other galaxies to about the eighth magnitude (drawn into the glow by their surface
//! brightness, dimmed by the home galaxy's dust; [`geoforge_cosmos::universe::sky`]), the
//! other stars and planets of its system, its moons and their phases, the eclipses to come,
//! and the world's constellations (drawn from its brightest stars when it has none of its
//! own; [`geoforge_cosmos::constellations`]). The stars are the galaxy's ([`Starfield`]) around the
//! system's place, with the world's edits; all else is derived from the system file.
//!
//! Coordinates are the world's own right ascension and declination: its equator, and its
//! northern spring equinox (the present) at RA 0.
//!
//! Planets shine by reflected starlight with geometric albedos by class (rocky planets 0.3 —
//! Earth 0.37, Venus 0.69, Mars 0.17, Mercury 0.14; moons 0.15 — the Moon 0.12, Callisto
//! 0.20, Europa 0.67; gas giants 0.5 — Jupiter 0.54, Saturn 0.50; ice giants 0.45 — Uranus
//! 0.49, Neptune 0.44; small bodies 0.1) and a Lambert sphere's phase law.

use crate::ephemeris::{self, Eclipse, Equatorial};
use crate::Result;
use geoforge_core::body::{BodyClass, NodeKind, SolarSystem, AU_KM};
use geoforge_core::galactic::{Constellation, StarId};
use geoforge_core::key::Key;
use geoforge_cosmos::sky::{self as stars, Neighbourhood, Observer, LY_PER_PC, NAKED_EYE};
use geoforge_cosmos::universe::sky::SkyGalaxy;
use geoforge_cosmos::{Galaxy, Starfield};
use geoforge_store::World;
use serde::Serialize;

/// Parsecs per AU.
const PC_PER_AU: f64 = 1.0 / 206_264.806;
/// A planet closer to the star than this in the sky is lost in its glare, degrees.
const GLARE_DEG: f64 = 15.0;
/// Galaxies are listed to this magnitude, about a pair of binoculars'.
pub const GALAXY_LIMIT: f64 = 8.0;

/// Another galaxy in the sky.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct GalaxyObject {
    pub name: String,
    /// Hubble type.
    pub kind: String,
    pub ra_deg: f64,
    pub dec_deg: f64,
    /// Integrated visual magnitude, through the home galaxy's dust.
    pub v: f64,
    /// Millions of light-years.
    pub distance_mly: f64,
    /// Across its major axis to its catalogued edge, degrees.
    pub size_deg: f64,
    pub axis_ratio: f64,
    /// Position angle of its major axis, degrees east of north.
    pub pa_deg: f64,
    /// Surface brightness at its centre, V magnitudes per square arcsecond (a dark sky is
    /// about 22).
    pub surface_brightness: f64,
}

impl GalaxyObject {
    /// What it takes to see it from a dark site, by its magnitude and how faint its light is
    /// spread: "naked eye", "binoculars", "telescope" or "long exposures" (the Large
    /// Magellanic Cloud at 21.4 mag/arcsec² is plain to the eye; a dwarf spheroidal at 25 is
    /// not, however bright it is in all).
    #[must_use]
    pub fn seen_with(&self) -> &'static str {
        match self.surface_brightness {
            mu if mu <= 22.5 && self.v <= stars::NAKED_EYE => "naked eye",
            mu if mu <= 23.0 => "binoculars",
            mu if mu <= 24.0 => "telescope",
            _ => "long exposures",
        }
    }
}

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
    /// A neighbourhood star's identity.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub id: Option<StarId>,
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
    /// Other galaxies, brightest first.
    pub galaxies: Vec<GalaxyObject>,
    pub eclipses: Vec<Eclipse>,
    /// The path of the system's brightest star through the year, as (RA, Dec).
    pub ecliptic: Vec<(f64, f64)>,
    #[serde(skip)]
    neighbourhood: Neighbourhood,
    #[serde(skip)]
    seen: Vec<SkyGalaxy>,
    #[serde(skip)]
    seed: geoforge_core::seed::Seed,
    #[serde(skip)]
    pub frame: Frame,
}

/// The stars of `world`'s galaxy as it has them: its own galaxy model if it has one, its
/// edits over them.
///
/// # Errors
/// If the world has no system, or its galaxy documents do not read.
pub fn starfield(world: &World) -> Result<Starfield> {
    let system = world.load_system()?.ok_or_else(|| crate::Error::Invalid("this world has no solar system (a hand-made body), so no galaxy around it".into()))?;
    // The world's own model, else its address's galaxy's (else the preset its site names).
    let galaxy: Option<Galaxy> = world.load_galaxy_model()?;
    let galaxy = galaxy.or_else(|| system.root.placement.and_then(|a| geoforge_cosmos::resolve::galaxy(&a).ok()).map(|g| g.model));
    Ok(Starfield::of_system(&system, galaxy, world.load_star_edits()?))
}

/// A world's sky frame: where it sees from in its galaxy, and its equator and equinox.
#[derive(Debug, Clone)]
pub struct Frame {
    pub observer: Observer,
    pub equatorial: Equatorial,
}

impl Frame {
    /// The frame of `world` in `system` among `field`'s stars.
    ///
    /// # Errors
    /// If `world` is not a body orbiting something in `system`.
    pub fn of(system: &SolarSystem, world: Key, field: &Starfield) -> std::result::Result<Self, String> {
        let equatorial = Equatorial::of(system, world).ok_or("the world orbits nothing")?;
        Ok(Self { observer: Observer::of_system(system, field.galaxy().reference_kpc), equatorial })
    }

    /// Where `p` (pc, in the galaxy) lies in this sky: RA, Dec (degrees) and distance (pc).
    #[must_use]
    pub fn radec_of(&self, p: [f64; 3]) -> (f64, f64, f64) {
        let v = [0, 1, 2].map(|k| p[k] - self.observer.position_pc[k]);
        let d = norm(v).max(1e-9);
        let (ra, dec) = self.equatorial.ra_dec(self.observer.to_system_frame(v.map(|x| x / d)));
        (ra, dec, d)
    }

    /// The place in the galaxy (pc) at RA and Dec (degrees), `distance_pc` away.
    #[must_use]
    pub fn galactic_of(&self, ra_deg: f64, dec_deg: f64, distance_pc: f64) -> [f64; 3] {
        let (ra, dec) = (ra_deg.to_radians(), dec_deg.to_radians());
        let local = [dec.cos() * ra.cos(), dec.cos() * ra.sin(), dec.sin()];
        let f = &self.equatorial;
        let dir = [0, 1, 2].map(|k| f.x[k] * local[0] + f.y[k] * local[1] + f.z[k] * local[2]);
        let u = self.observer.to_galaxy_frame(dir);
        [0, 1, 2].map(|k| self.observer.position_pc[k] + distance_pc * u[k])
    }
}

/// The sky of `world` on `day` among `field`'s stars, with the eclipses over the
/// `eclipse_days` that follow.
///
/// # Errors
/// If `world` is not a body orbiting something in `system`.
pub fn sky(system: &SolarSystem, world: Key, day: f64, eclipse_days: f64, field: &Starfield) -> std::result::Result<Sky, String> {
    let place = Frame::of(system, world, field)?;
    let frame = place.equatorial;
    let neighbourhood = stars::neighbourhood(field, &place.observer, NAKED_EYE);
    let radius_kpc = neighbourhood.radius_kpc;
    let stars = neighbourhood
        .stars
        .iter()
        .map(|s| {
            let (ra, dec) = frame.ra_dec(s.direction);
            Object {
                id: Some(s.id),
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
    let observer = &place.observer;
    let seen: Vec<SkyGalaxy> = system
        .root
        .site
        .as_ref()
        .map(|site| geoforge_cosmos::universe::sky::seen(site, system.root.placement.as_ref(), system.root.seed, field.galaxy(), observer, GALAXY_LIMIT + 1.0))
        .unwrap_or_default()
        .into_iter()
        .map(|g| {
            let dust = neighbourhood.extinction_out(observer.to_system_frame(g.direction));
            g.dimmed(dust)
        })
        .filter(|g| g.v <= GALAXY_LIMIT)
        .collect();
    let galaxies = seen.iter().map(|g| galaxy_object(g, observer, &frame)).collect();
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
    Ok(Sky {
        day,
        galaxy: neighbourhood.galaxy.clone(),
        radius_kpc,
        stars,
        bodies,
        galaxies,
        eclipses: ephemeris::eclipses(system, world, day, eclipse_days),
        ecliptic,
        neighbourhood,
        seen,
        seed: system.root.seed,
        frame: place,
    })
}

/// Where `g` is in the sky, and how its major axis lies.
fn galaxy_object(g: &SkyGalaxy, observer: &Observer, frame: &Equatorial) -> GalaxyObject {
    let toward = observer.to_system_frame(g.direction);
    let (ra, dec) = frame.ra_dec(toward);
    let major = observer.to_system_frame(g.major);
    let local = [dot(major, frame.x), dot(major, frame.y), dot(major, frame.z)];
    let (a, d) = (ra.to_radians(), dec.to_radians());
    let north = [-d.sin() * a.cos(), -d.sin() * a.sin(), d.cos()];
    let east = [-a.sin(), a.cos(), 0.0];
    GalaxyObject {
        name: g.name.clone(),
        kind: g.kind.clone(),
        ra_deg: ra,
        dec_deg: dec,
        v: g.v,
        distance_mly: g.distance_pc * LY_PER_PC / 1e6,
        size_deg: g.size_deg(),
        axis_ratio: g.axis_ratio,
        pa_deg: dot(local, east).atan2(dot(local, north)).to_degrees().rem_euclid(180.0),
        surface_brightness: 26.4 - 2.5 * g.sigma0.max(1e-12).log10(),
    }
}

/// A constellation's lines on the sky, (RA, Dec) at each end, degrees.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Figure {
    pub name: String,
    pub lines: Vec<[(f64, f64); 2]>,
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
            let f = &self.frame.equatorial;
            let dir = [0, 1, 2].map(|k| f.x[k] * local[0] + f.y[k] * local[1] + f.z[k] * local[2]);
            let galaxies: f64 = self.seen.iter().map(|g| g.sigma(self.frame.observer.to_galaxy_frame(dir))).sum();
            *v = (self.neighbourhood.glow(dir) + galaxies / self.neighbourhood.pole_sigma()) as f32;
        }
        let top = out.iter().copied().fold(f32::MIN, f32::max).max(1e-9);
        out.iter_mut().for_each(|v| *v /= top);
        out
    }
}

impl Sky {
    /// Where `p` (pc, in the galaxy) lies in this sky: RA, Dec (degrees) and distance (pc).
    #[must_use]
    pub fn radec_of(&self, p: [f64; 3]) -> (f64, f64, f64) {
        self.frame.radec_of(p)
    }

    /// The place in the galaxy (pc) at RA and Dec (degrees), `distance_pc` away.
    #[must_use]
    pub fn galactic_of(&self, ra_deg: f64, dec_deg: f64, distance_pc: f64) -> [f64; 3] {
        self.frame.galactic_of(ra_deg, dec_deg, distance_pc)
    }

    /// Figures drawn from this sky's brightest stars, for a world with none of its own.
    #[must_use]
    pub fn automatic_constellations(&self) -> Vec<Constellation> {
        geoforge_cosmos::constellations::automatic(&self.neighbourhood.stars, self.seed)
    }

    /// The world's `constellations`, or if it has none the automatic ones.
    #[must_use]
    pub fn constellations_or_automatic(&self, constellations: Vec<Constellation>) -> Vec<Constellation> {
        if constellations.is_empty() {
            self.automatic_constellations()
        } else {
            constellations
        }
    }

    /// The figures of `constellations`, their stars found in this sky or else in `field`.
    #[must_use]
    pub fn figures(&self, constellations: &[Constellation], field: &Starfield) -> Vec<Figure> {
        let at = |id: &StarId| {
            self.stars.iter().find(|o| o.id == Some(*id)).map(|o| (o.ra_deg, o.dec_deg)).or_else(|| {
                field.star(id).map(|s| {
                    let (ra, dec, _) = self.radec_of(s.position_pc);
                    (ra, dec)
                })
            })
        };
        constellations.iter().map(|c| Figure { name: c.name.clone(), lines: c.lines.iter().filter_map(|[a, b]| Some([at(a)?, at(b)?])).collect() }).collect()
    }
}

/// RA and Dec at the centre of chart pixel (x, y): RA 360 at the left edge, 0 at the right.
#[must_use]
pub fn chart_coords(x: usize, y: usize, width: usize, height: usize) -> (f64, f64) {
    (360.0 * (1.0 - (x as f64 + 0.5) / width as f64), 90.0 - 180.0 * (y as f64 + 0.5) / height as f64)
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
            id: None,
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
            id: None,
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

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
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
        let field = Starfield::of_system(&sys, None, geoforge_core::galactic::StarEdits::default());
        let s = sky(&sys, Key(3), 0.0, 400.0, &field).unwrap();
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
        // Positions turn both ways, and a constellation of the brightest stars is drawn.
        let star = &s.stars[0];
        let p = s.galactic_of(star.ra_deg, star.dec_deg, star.distance / LY_PER_PC);
        let (ra, dec, d) = s.radec_of(p);
        assert!((ra - star.ra_deg).abs() < 1e-6 && (dec - star.dec_deg).abs() < 1e-6 && (d * LY_PER_PC - star.distance).abs() < 1e-6);
        let c = Constellation { name: "Kite".into(), lines: vec![[s.stars[0].id.unwrap(), s.stars[1].id.unwrap()], [s.stars[1].id.unwrap(), s.stars[2].id.unwrap()]] };
        let f = s.figures(&[c], &field);
        assert_eq!(f[0].lines.len(), 2);
        assert_eq!(f[0].lines[0][0], (s.stars[0].ra_deg, s.stars[0].dec_deg));
    }

    /// From the real Earth the Milky Way crosses the sky where it does: bright through
    /// Sagittarius and Cygnus, faint at the galactic poles; and the sky holds about the Sun's
    /// 9,110 naked-eye stars.
    #[test]
    fn the_milky_way_crosses_the_real_sky() {
        let (sys, body) = geoforge_cosmos::sol::archetype(geoforge_core::Seed(1), Archetype::Earth);
        let field = Starfield::of_system(&sys, None, geoforge_core::galactic::StarEdits::default());
        let s = sky(&sys, body.key, 0.0, 30.0, &field).unwrap();
        assert!((7_000..12_000).contains(&s.stars.len()), "{} naked-eye stars", s.stars.len());
        let glow = s.glow(360, 180);
        let at = |ra: f64, dec: f64| glow[((90.0 - dec).floor() as usize).min(179) * 360 + ((360.0 - ra).floor() as usize).min(359)];
        let best =
            |ra: f64, dec: f64| (-6..=6).flat_map(|i| (-6..=6).map(move |j| (ra + f64::from(i), dec + f64::from(j)))).map(|(r, d)| at(r.rem_euclid(360.0), d.clamp(-89.0, 89.0))).fold(0.0, f32::max);
        // Toward the centre (Sagittarius), in Cygnus, at the poles (Coma Berenices, Sculptor).
        let (centre, cygnus, north, south) = (best(270.0, -28.0), best(305.0, 40.0), at(192.9, 27.1), at(12.9, -27.1));
        assert!(centre > 5.0 * north && centre > 5.0 * south, "centre {centre}, poles {north} {south}");
        assert!(cygnus > 3.0 * north, "Cygnus {cygnus}, pole {north}");
        // The Local Group where it is: the Clouds as patches brighter than the sky around
        // them, Andromeda a dimmer smudge.
        let named = |n: &str| s.galaxies.iter().find(|g| g.name == n).unwrap_or_else(|| panic!("no {n}"));
        let (lmc, m31) = (named("Large Magellanic Cloud"), named("Andromeda Galaxy"));
        assert!((lmc.dec_deg + 69.76).abs() < 0.5 && (m31.dec_deg - 41.27).abs() < 0.5 && (m31.v - 3.44).abs() < 0.3, "{lmc:?} {m31:?}");
        let (cloud, beside) = (at(lmc.ra_deg, lmc.dec_deg), at(lmc.ra_deg + 40.0, lmc.dec_deg + 15.0));
        assert!(cloud > 2.0 * beside, "LMC {cloud}, beside it {beside}");
        assert!(at(m31.ra_deg, m31.dec_deg) > at(m31.ra_deg + 10.0, m31.dec_deg), "Andromeda");
        // At their centres the Large Cloud is 21.4 mag/arcsec² (observed about 21), Andromeda's
        // disk about 20; both plain to the eye, Centaurus A not.
        assert!((lmc.surface_brightness - 21.4).abs() < 0.5 && (m31.surface_brightness - 20.3).abs() < 0.7, "{} {}", lmc.surface_brightness, m31.surface_brightness);
        assert_eq!([lmc.seen_with(), m31.seen_with(), named("Centaurus A").seen_with()], ["naked eye", "naked eye", "binoculars"]);
        assert!(s.automatic_constellations().len() > 40);
    }
}
