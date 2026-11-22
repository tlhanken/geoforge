//! The Solar System from published data, set at a moment of its real calendar.
//!
//! - **Planets and Pluto**: JPL's Keplerian elements and their rates for approximate
//!   positions, 1800–2050 (Standish; ssd.jpl.nasa.gov/planets/approx_pos.html), on the J2000
//!   ecliptic. They hold heliocentric longitudes to 15–600″, so a planet's place in the sky to
//!   arcminutes. Earth's row is the Earth–Moon barycentre's.
//! - **The Moon**: Meeus's mean lunar elements (*Astronomical Algorithms*, ch. 47) at the
//!   moment, moved from the equinox of date to J2000, with its node regressing and perigee
//!   advancing, so eclipse seasons drift as they do. Its periodic terms (evection 1.3°,
//!   variation 0.7°) are left out: eclipses fall within hours of the real ones.
//! - **Other moons**: JPL's mean elements (epoch J2000, ssd.jpl.nasa.gov/sats/elem), measured
//!   on each moon's Laplace plane or its planet's equator, carried to the moment along their
//!   periods with their nodes regressing, then turned onto the ecliptic (and then fixed).
//! - **Bodies**: masses, radii, spins and tilts from NASA's planetary fact sheets and JPL's
//!   satellite physical parameters; north poles from the IAU's rotational elements
//!   (Archinal et al. 2018).
//! - **Belts and rings**: the main belt (2.1–3.3 AU, ~1.5 million asteroids over 1 km), the
//!   Kuiper belt (30–50 AU, ~10⁵ objects over 100 km) and Saturn's main rings (C to A,
//!   1.5×10¹⁹ kg of ice).
//!
//! Orbits here are fixed ellipses, so every orbit's size is set from its measured period
//! (the two-body period then is the real one): the Moon's from the sidereal month, which the
//! Sun's pull makes 0.1% longer than its true size gives, a drift of 6° a year otherwise.
//! The neighbourhood stars of the sky are still drawn from the galaxy model.
//!
//! A world's present is its northern spring equinox ([`archetype`]): for Earth the March
//! equinox of 2026, for another body its own last one before that.

use geoforge_core::body::{
    BodyClass, BodyContext, BodyExterior, Composition, Epoch, GalacticSite, NodeKind, OrbitalNode, Population, PowerLaw, Precession, SolarSystem, Species, Star, StellarStage, SystemRoot, AU_KM,
    SOLAR_MASS_KG,
};
use geoforge_core::derived::orbital_period_days;
use geoforge_core::{Archetype, Elements, Key, Seed};

/// The 2026 March equinox, 2026-03-20 14:46 UTC, as a Julian date (TT = UTC + 69.184 s).
pub const EQUINOX_2026_JD_TT: f64 = 2_461_120.116_08;
const J2000_JD: f64 = 2_451_545.0;
/// Obliquity of the ecliptic at J2000, degrees (the IAU's 84 381.406″).
const OBLIQUITY_J2000_DEG: f64 = 23.439_279;
/// General precession in longitude, degrees per Julian century (5028.796″).
const PRECESSION_DEG_CY: f64 = 1.396_888;
const SIDEREAL_MONTH_DAYS: f64 = 27.321_661;
/// The Moon's node regresses once in 18.6 years and its perigee advances once in 8.85
/// (Meeus's rates for Ω and L′ − M′, per day).
const MOON_PRECESSION: Precession = Precession { node_deg_day: -1_934.136_289_1 / 36525.0, periapsis_deg_day: (481_267.881_234_21 - 477_198.867_505_5) / 36525.0 };
/// The Moon's equator is tilted this far from the ecliptic, its pole on the far side of the
/// ecliptic pole from the orbit's (Cassini's laws), degrees.
const MOON_EQUATOR_TILT_DEG: f64 = 1.5427;

struct Planet {
    name: &'static str,
    /// a (AU), e, I, L, ϖ, Ω (degrees) at J2000.
    at: [f64; 6],
    /// Their rates per Julian century.
    rate: [f64; 6],
    class: BodyClass,
    mass_kg: f64,
    radius_km: f64,
    rotation_h: f64,
    obliquity_deg: f64,
    /// Spin pole (the right-hand rule's north), right ascension and declination on the ICRF
    /// (J2000 equator), degrees.
    pole: (f64, f64),
}

use BodyClass::{GasGiant, IceGiant, SmallBody, Terrestrial};

#[rustfmt::skip]
const PLANETS: [Planet; 9] = [
    Planet { name: "Mercury", at: [0.38709927, 0.20563593, 7.00497902, 252.25032350, 77.45779628, 48.33076593], rate: [0.00000037, 0.00001906, -0.00594749, 149472.67411175, 0.16047689, -0.12534081], class: Terrestrial, mass_kg: 3.3011e23, radius_km: 2439.7, rotation_h: 1407.6, obliquity_deg: 0.034, pole: (281.0103, 61.4155) },
    // Retrograde spin: tilts past 90°.
    Planet { name: "Venus", at: [0.72333566, 0.00677672, 3.39467605, 181.97909950, 131.60246718, 76.67984255], rate: [0.00000390, -0.00004107, -0.00078890, 58517.81538729, 0.00268329, -0.27769418], class: Terrestrial, mass_kg: 4.8675e24, radius_km: 6051.8, rotation_h: 5832.5, obliquity_deg: 177.36, pole: (92.76, -67.16) },
    Planet { name: "Earth", at: [1.00000261, 0.01671123, -0.00001531, 100.46457166, 102.93768193, 0.0], rate: [0.00000562, -0.00004392, -0.01294668, 35999.37244981, 0.32327364, 0.0], class: Terrestrial, mass_kg: 5.9722e24, radius_km: 6371.0, rotation_h: 23.9345, obliquity_deg: 23.44, pole: (0.0, 90.0) },
    Planet { name: "Mars", at: [1.52371034, 0.09339410, 1.84969142, -4.55343205, -23.94362959, 49.55953891], rate: [0.00001847, 0.00007882, -0.00813131, 19140.30268499, 0.44441088, -0.29257343], class: Terrestrial, mass_kg: 6.4171e23, radius_km: 3389.5, rotation_h: 24.6229, obliquity_deg: 25.19, pole: (317.681, 52.887) },
    Planet { name: "Jupiter", at: [5.20288700, 0.04838624, 1.30439695, 34.39644051, 14.72847983, 100.47390909], rate: [-0.00011607, -0.00013253, -0.00183714, 3034.74612775, 0.21252668, 0.20469106], class: GasGiant, mass_kg: 1.89819e27, radius_km: 69911.0, rotation_h: 9.925, obliquity_deg: 3.13, pole: (268.057, 64.495) },
    Planet { name: "Saturn", at: [9.53667594, 0.05386179, 2.48599187, 49.95424423, 92.59887831, 113.66242448], rate: [-0.00125060, -0.00050991, 0.00193609, 1222.49362201, -0.41897216, -0.28867794], class: GasGiant, mass_kg: 5.6834e26, radius_km: 58232.0, rotation_h: 10.656, obliquity_deg: 26.73, pole: (40.589, 83.537) },
    Planet { name: "Uranus", at: [19.18916464, 0.04725744, 0.77263783, 313.23810451, 170.95427630, 74.01692503], rate: [-0.00196176, -0.00004397, -0.00242939, 428.48202785, 0.40805281, 0.04240589], class: IceGiant, mass_kg: 8.6813e25, radius_km: 25362.0, rotation_h: 17.24, obliquity_deg: 97.77, pole: URANUS_POLE },
    Planet { name: "Neptune", at: [30.06992276, 0.00859048, 1.77004347, -55.12002969, 44.96476227, 131.78422574], rate: [0.00026291, 0.00005105, 0.00035372, 218.45945325, -0.32241464, -0.00508664], class: IceGiant, mass_kg: 1.02413e26, radius_km: 24622.0, rotation_h: 16.11, obliquity_deg: 28.32, pole: (299.36, 43.46) },
    Planet { name: "Pluto", at: [39.48211675, 0.24882730, 17.14001206, 238.92903833, 224.06891629, 110.30393684], rate: [-0.00031596, 0.00005170, 0.00004818, 145.20780515, -0.04062942, -0.01183482], class: Terrestrial, mass_kg: 1.303e22, radius_km: 1188.3, rotation_h: 153.29, obliquity_deg: 122.53, pole: (132.993, -6.163) },
];

struct Moon {
    name: &'static str,
    planet: &'static str,
    /// JPL mean elements at J2000: e, ω, M, i, node (degrees), sidereal period (days), nodal
    /// regression period (years; 0 for none).
    el: [f64; 7],
    /// Pole of the plane they are measured on (Laplace plane or the planet's equator), right
    /// ascension and declination on the ICRF, degrees.
    plane: (f64, f64),
    class: BodyClass,
    mass_kg: f64,
    radius_km: f64,
}

/// Uranus spins backward: its spin pole is the IAU north pole's antipode (as is Venus's).
const URANUS_POLE: (f64, f64) = (77.311, 15.175);
const NEPTUNE_POLE: (f64, f64) = (299.36, 43.46);

#[rustfmt::skip]
const MOONS: [Moon; 20] = [
    Moon { name: "Phobos", planet: "Mars", el: [0.015, 216.3, 189.7, 1.1, 169.2, 0.3187, 2.3], plane: (317.7, 52.9), class: SmallBody, mass_kg: 1.0659e16, radius_km: 11.08 },
    Moon { name: "Deimos", planet: "Mars", el: [0.000, 0.0, 205.0, 1.8, 54.3, 1.2625, 56.2], plane: (316.6, 53.5), class: SmallBody, mass_kg: 1.4762e15, radius_km: 6.2 },
    Moon { name: "Io", planet: "Jupiter", el: [0.004, 49.1, 330.9, 0.0, 0.0, 1.762732, 0.0], plane: (268.1, 64.5), class: Terrestrial, mass_kg: 8.9319e22, radius_km: 1821.5 },
    Moon { name: "Europa", planet: "Jupiter", el: [0.009, 45.0, 345.4, 0.5, 184.0, 3.525463, 30.202], plane: (268.1, 64.5), class: Terrestrial, mass_kg: 4.7998e22, radius_km: 1560.8 },
    Moon { name: "Ganymede", planet: "Jupiter", el: [0.001, 198.3, 324.8, 0.2, 58.5, 7.155588, 137.812], plane: (268.2, 64.6), class: Terrestrial, mass_kg: 1.4819e23, radius_km: 2631.2 },
    Moon { name: "Callisto", planet: "Jupiter", el: [0.007, 43.8, 87.4, 0.3, 309.1, 16.690440, 577.264], plane: (268.7, 64.8), class: Terrestrial, mass_kg: 1.0759e23, radius_km: 2410.3 },
    Moon { name: "Mimas", planet: "Saturn", el: [0.020, 160.4, 275.3, 1.6, 66.2, 0.942422, 0.986], plane: (40.6, 83.5), class: Terrestrial, mass_kg: 3.751e19, radius_km: 198.2 },
    Moon { name: "Enceladus", planet: "Saturn", el: [0.005, 119.5, 57.0, 0.0, 0.0, 1.370218, 0.0], plane: (40.6, 83.5), class: Terrestrial, mass_kg: 1.0805e20, radius_km: 252.1 },
    Moon { name: "Tethys", planet: "Saturn", el: [0.001, 335.3, 0.0, 1.1, 273.0, 1.887802, 4.982], plane: (40.6, 83.5), class: Terrestrial, mass_kg: 6.175e20, radius_km: 531.0 },
    Moon { name: "Dione", planet: "Saturn", el: [0.002, 116.0, 212.0, 0.0, 0.0, 2.736916, 0.0], plane: (40.6, 83.5), class: Terrestrial, mass_kg: 1.0955e21, radius_km: 561.4 },
    Moon { name: "Rhea", planet: "Saturn", el: [0.001, 44.3, 31.5, 0.3, 133.7, 4.517503, 35.775], plane: (40.6, 83.5), class: Terrestrial, mass_kg: 2.3065e21, radius_km: 763.5 },
    Moon { name: "Titan", planet: "Saturn", el: [0.029, 78.3, 11.7, 0.3, 78.6, 15.945448, 687.370], plane: (36.4, 84.0), class: Terrestrial, mass_kg: 1.3452e23, radius_km: 2574.7 },
    Moon { name: "Iapetus", planet: "Saturn", el: [0.028, 254.5, 74.8, 7.6, 86.5, 79.331002, 3130.302], plane: (288.7, 78.9), class: Terrestrial, mass_kg: 1.8057e21, radius_km: 734.4 },
    Moon { name: "Miranda", planet: "Uranus", el: [0.001, 154.8, 73.0, 4.4, 100.9, 1.413479, 17.787], plane: URANUS_POLE, class: Terrestrial, mass_kg: 6.44e19, radius_km: 235.8 },
    Moon { name: "Ariel", planet: "Uranus", el: [0.001, 9.6, 193.5, 0.0, 0.0, 2.520379, 0.0], plane: URANUS_POLE, class: Terrestrial, mass_kg: 1.251e21, radius_km: 578.9 },
    Moon { name: "Umbriel", planet: "Uranus", el: [0.004, 183.4, 253.0, 0.1, 174.8, 4.144177, 129.945], plane: URANUS_POLE, class: Terrestrial, mass_kg: 1.275e21, radius_km: 584.7 },
    Moon { name: "Titania", planet: "Uranus", el: [0.002, 184.0, 68.1, 0.1, 29.5, 8.705869, 1644.649], plane: URANUS_POLE, class: Terrestrial, mass_kg: 3.400e21, radius_km: 788.9 },
    Moon { name: "Oberon", planet: "Uranus", el: [0.002, 132.2, 143.6, 0.1, 76.8, 13.463237, 192.798], plane: URANUS_POLE, class: Terrestrial, mass_kg: 3.076e21, radius_km: 761.4 },
    // Retrograde: inclined 157° to Neptune's equator.
    Moon { name: "Triton", planet: "Neptune", el: [0.000, 0.0, 63.0, 157.3, 177.7, 5.8792, 299.4], plane: NEPTUNE_POLE, class: Terrestrial, mass_kg: 2.139e22, radius_km: 1352.6 },
    Moon { name: "Charon", planet: "Pluto", el: [0.000, 0.0, 304.1, 0.0, 0.0, 6.387222, 0.0], plane: (132.993, -6.163), class: Terrestrial, mass_kg: 1.586e21, radius_km: 606.0 },
];

/// The Solar System `jd_tt` (a Julian date, TT) as its present. `seed` only draws the
/// neighbourhood stars of its sky.
#[must_use]
pub fn at(seed: Seed, jd_tt: f64) -> SolarSystem {
    let t = centuries(jd_tt);
    let mut nodes = Vec::new();
    let mut next = 1u64;
    let mut key = || {
        next += 1;
        Key(next)
    };
    let sun = key();
    let star = Star {
        mass_solar: 1.0,
        luminosity_solar: 1.0,
        radius_solar: 1.0,
        temperature_k: 5772.0,
        age_gyr: 4.57,
        metallicity_dex: 0.0,
        spectral_class: "G2V".into(),
        stage: StellarStage::MainSequence,
        initial_mass_solar: None,
    };
    nodes.push(OrbitalNode { key: sun, name: "Sun".into(), parent: None, node: NodeKind::Star(star), elements: None, precession: None });
    let sun_kg = SOLAR_MASS_KG;
    for p in &PLANETS {
        let planet = key();
        let ext = BodyExterior { class: p.class, mass_kg: p.mass_kg, radius_km: p.radius_km, rotation_h: p.rotation_h, obliquity_deg: p.obliquity_deg };
        let mut el = planet_elements(p, t);
        el.a_au = size_for_period(360.0 * 36525.0 / p.rate[3], sun_kg + p.mass_kg);
        nodes.push(OrbitalNode { key: planet, name: p.name.into(), parent: Some(sun), node: NodeKind::Body(ext), elements: Some(el), precession: None });
        if p.name == "Earth" {
            let (mass_kg, radius_km) = (7.342e22, 1737.4);
            let ext = BodyExterior { class: Terrestrial, mass_kg, radius_km, rotation_h: SIDEREAL_MONTH_DAYS * 24.0, obliquity_deg: 6.68 };
            let mut el = moon_elements(t);
            el.a_au = size_for_period(SIDEREAL_MONTH_DAYS, p.mass_kg + mass_kg);
            nodes.push(OrbitalNode { key: key(), name: "Moon".into(), parent: Some(planet), node: NodeKind::Body(ext), elements: Some(el), precession: Some(MOON_PRECESSION) });
        }
        if p.name == "Saturn" {
            // C ring's inner edge to A ring's outer edge. Particles 1 cm–10 m, dN/dD ∝ D⁻³,
            // holding the rings' 1.54×10¹⁹ kg (Iess et al. 2019) as ice: ~1.6×10¹⁹ of them.
            let rings = Population {
                inner_au: 74_658.0 / AU_KM,
                outer_au: 136_775.0 / AU_KM,
                count_estimate: 1.6e19,
                size_law: PowerLaw { min_km: 1e-5, max_km: 0.01, exponent: 3.0 },
                composition: Composition::from_weights(&[(Species::Ice, 0.95), (Species::Silicate, 0.05)]),
            };
            nodes.push(OrbitalNode { key: key(), name: "Saturn's rings".into(), parent: Some(planet), node: NodeKind::Ring(rings), elements: None, precession: None });
        }
        for m in MOONS.iter().filter(|m| m.planet == p.name) {
            let ext = BodyExterior { class: m.class, mass_kg: m.mass_kg, radius_km: m.radius_km, rotation_h: m.el[5] * 24.0, obliquity_deg: 0.0 };
            let mut el = satellite_elements(m, jd_tt);
            el.a_au = size_for_period(m.el[5], p.mass_kg + m.mass_kg);
            nodes.push(OrbitalNode { key: key(), name: m.name.into(), parent: Some(planet), node: NodeKind::Body(ext), elements: Some(el), precession: None });
        }
        if p.name == "Mars" {
            let belt = Population {
                inner_au: 2.1,
                outer_au: 3.3,
                count_estimate: 1.5e6,
                size_law: PowerLaw { min_km: 1.0, max_km: 939.4, exponent: 3.0 },
                composition: Composition::from_weights(&[(Species::Silicate, 0.75), (Species::Iron, 0.15), (Species::Ice, 0.10)]),
            };
            nodes.push(OrbitalNode { key: key(), name: "Main Belt".into(), parent: Some(sun), node: NodeKind::Belt(belt), elements: None, precession: None });
        }
    }
    // ~10⁵ over 100 km; with dN/dD ∝ D⁻³ that is ~10⁹ over 1 km. Pluto is the largest.
    let kuiper = Population {
        inner_au: 30.0,
        outer_au: 50.0,
        count_estimate: 1e9,
        size_law: PowerLaw { min_km: 1.0, max_km: 2376.6, exponent: 3.0 },
        composition: Composition::from_weights(&[(Species::Ice, 0.6), (Species::Silicate, 0.4)]),
    };
    nodes.push(OrbitalNode { key: key(), name: "Kuiper Belt".into(), parent: Some(sun), node: NodeKind::Belt(kuiper), elements: None, precession: None });
    let root = SystemRoot {
        key: Key(1),
        seed,
        name: "Solar System".into(),
        epoch: Epoch { system_age_gyr: 4.568 },
        placement: None,
        site: Some(GalacticSite { galaxy: "Milky Way".into(), radius_kpc: 8.2 }),
        present_jd_tt: Some(jd_tt),
    };
    SolarSystem { root, nodes }
}

/// The Solar System and the archetype's body in it, at that body's present: Earth's 2026
/// March equinox for Earth, otherwise the body's last northern spring equinox before it. The
/// context is the archetype's, with the body's key and its real orbit.
#[must_use]
pub fn archetype(seed: Seed, a: Archetype) -> (SolarSystem, BodyContext) {
    let name = a.real_body();
    let jd = if name == "Earth" { EQUINOX_2026_JD_TT } else { northern_equinox_before(name, EQUINOX_2026_JD_TT).unwrap_or(EQUINOX_2026_JD_TT) };
    let mut system = at(seed, jd);
    let key = system.by_name(name).expect("every archetype is a Solar System body").key;
    let body = known_context(&system, key).expect("an archetype's body has its context");
    if let Some(node) = system.nodes.iter_mut().find(|n| n.key == key) {
        node.node = NodeKind::Body(body.exterior);
    }
    (system, body)
}

/// The measured context of body `key` of the real Solar System, for a body with an
/// archetype (Earth, Mars, Venus, Mercury, the Moon, Titan, Europa, Jupiter, Neptune);
/// `None` for any other body or system.
#[must_use]
pub fn known_context(system: &SolarSystem, key: Key) -> Option<BodyContext> {
    system.root.present_jd_tt?;
    let node = system.node(key)?;
    let a = Archetype::ALL.into_iter().find(|a| a.real_body() == node.name)?;
    let mut body = a.context(key, &node.name);
    body.orbit = node.elements?;
    Some(body)
}

/// The last northern spring equinox of `name` (a planet, Pluto or a moon) at or before
/// `jd_tt`: when the Sun crosses its equator northward. `None` for an unknown body.
#[must_use]
pub fn northern_equinox_before(name: &str, jd_tt: f64) -> Option<f64> {
    let host = MOONS.iter().find(|m| m.name == name).map_or(name, |m| m.planet);
    let host = if name == "Moon" { "Earth" } else { host };
    let planet = PLANETS.iter().find(|p| p.name == host)?;
    // The Sun's height above the body's equator, as a sine.
    let f = |jd: f64| -> Option<f64> {
        let r = heliocentric(planet, jd);
        Some(-dot(r, pole_ecliptic(name, jd)?) / norm(r))
    };
    let year = 360.0 * 36525.0 / planet.rate[3];
    // A moon's equator meets the Sun twice a month too (the Moon's), so step finely enough.
    let step = if name == "Moon" { 0.5 } else { year / 400.0 };
    let (mut hi, mut lo) = (jd_tt, jd_tt - step);
    let limit = jd_tt - 1.05 * year;
    while !(f(lo)? < 0.0 && f(hi)? >= 0.0) {
        hi = lo;
        lo -= step;
        if lo < limit {
            return None;
        }
    }
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if f(mid)? < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Some(0.5 * (lo + hi))
}

fn centuries(jd_tt: f64) -> f64 {
    (jd_tt - J2000_JD) / 36525.0
}

fn planet_elements(p: &Planet, t: f64) -> Elements {
    let v: Vec<f64> = (0..6).map(|k| p.at[k] + p.rate[k] * t).collect();
    let (a_au, e, i_deg, l, peri, node) = (v[0], v[1], v[2], v[3], v[4], v[5]);
    Elements { a_au, e, i_deg, lan_deg: node.rem_euclid(360.0), argp_deg: (peri - node).rem_euclid(360.0), m0_deg: (l - peri).rem_euclid(360.0) }
}

fn heliocentric(p: &Planet, jd_tt: f64) -> [f64; 3] {
    let e = planet_elements(p, centuries(jd_tt));
    e.position_at(e.m0_deg)
}

/// Meeus's mean lunar elements (eq. 47.1–47.5) at `t` centuries, moved to the J2000 equinox.
fn moon_elements(t: f64) -> Elements {
    let p = PRECESSION_DEG_CY * t;
    let lon = 218.316_447_7 + 481_267.881_234_21 * t - 0.001_578_6 * t * t - p;
    let anomaly = 134.963_396_4 + 477_198.867_505_5 * t + 0.008_741_4 * t * t;
    let node = moon_node(t);
    let peri = lon - anomaly;
    Elements { a_au: 0.0, e: 0.0549, i_deg: 5.145_396, lan_deg: node.rem_euclid(360.0), argp_deg: (peri - node).rem_euclid(360.0), m0_deg: anomaly.rem_euclid(360.0) }
}

/// The Moon's ascending node on the J2000 ecliptic, degrees.
fn moon_node(t: f64) -> f64 {
    125.044_547_9 - 1_934.136_289_1 * t - PRECESSION_DEG_CY * t
}

/// JPL's mean elements carried to `jd_tt` (mean longitude along the period, the node
/// regressing) and turned from the moon's reference plane onto the J2000 ecliptic.
fn satellite_elements(m: &Moon, jd_tt: f64) -> Elements {
    let [e, w0, m0, i, node0, period, pnode] = m.el;
    let days = jd_tt - J2000_JD;
    let lon = node0 + w0 + m0 + 360.0 * days / period;
    let node = if pnode > 0.0 { node0 - 360.0 * days / (pnode * 365.25) } else { node0 };
    let peri = node0 + w0;
    let local = Elements { a_au: 1.0, e, i_deg: i, lan_deg: node, argp_deg: peri - node, m0_deg: lon - peri };
    to_ecliptic(&local, m.plane)
}

/// Elements measured on a plane with pole `(ra, dec)` (from its ascending node on the ICRF
/// equator) as elements on the J2000 ecliptic.
fn to_ecliptic(local: &Elements, (ra, dec): (f64, f64)) -> Elements {
    let z = icrf(ra, dec);
    let x = [-ra.to_radians().sin(), ra.to_radians().cos(), 0.0];
    let y = cross(z, x);
    let turn = |v: [f64; 3]| ecliptic([x[0] * v[0] + y[0] * v[1] + z[0] * v[2], x[1] * v[0] + y[1] * v[1] + z[1] * v[2], x[2] * v[0] + y[2] * v[1] + z[2] * v[2]]);
    let (node, inc) = (local.lan_deg.to_radians(), local.i_deg.to_radians());
    let normal = turn([node.sin() * inc.sin(), -node.cos() * inc.sin(), inc.cos()]);
    // The periapsis (true anomaly 0) direction.
    let unit = Elements { a_au: 1.0, e: 0.0, m0_deg: 0.0, ..*local };
    let peri = turn(unit.position_at(0.0));
    let i_deg = normal[2].clamp(-1.0, 1.0).acos().to_degrees();
    let lan = if normal[0].hypot(normal[1]) < 1e-12 { 0.0 } else { normal[0].atan2(-normal[1]) };
    let ascending = [lan.cos(), lan.sin(), 0.0];
    let argp = dot(cross(ascending, peri), normal).atan2(dot(ascending, peri));
    Elements { a_au: local.a_au, e: local.e, i_deg, lan_deg: lan.to_degrees().rem_euclid(360.0), argp_deg: argp.to_degrees().rem_euclid(360.0), m0_deg: local.m0_deg.rem_euclid(360.0) }
}

/// A body's north pole on the J2000 ecliptic at `jd_tt`: Earth's precessing, the Moon's by
/// Cassini's laws from its node, the rest fixed.
fn pole_ecliptic(name: &str, jd_tt: f64) -> Option<[f64; 3]> {
    let t = centuries(jd_tt);
    match name {
        "Earth" => {
            let (lon, tilt) = ((90.0 - PRECESSION_DEG_CY * t).to_radians(), OBLIQUITY_J2000_DEG.to_radians());
            Some([tilt.sin() * lon.cos(), tilt.sin() * lon.sin(), tilt.cos()])
        }
        "Moon" => {
            let (lon, tilt) = ((moon_node(t) + 90.0).to_radians(), MOON_EQUATOR_TILT_DEG.to_radians());
            Some([tilt.sin() * lon.cos(), tilt.sin() * lon.sin(), tilt.cos()])
        }
        _ => {
            let (ra, dec) = PLANETS.iter().find(|p| p.name == name).map(|p| p.pole).or_else(|| MOONS.iter().find(|m| m.name == name).map(|m| m.plane))?;
            Some(ecliptic(icrf(ra, dec)))
        }
    }
}

/// The semi-major axis (AU) whose two-body period around `mass_kg` is `days`.
fn size_for_period(days: f64, mass_kg: f64) -> f64 {
    (days / orbital_period_days(1.0, mass_kg)).powf(2.0 / 3.0)
}

fn icrf(ra: f64, dec: f64) -> [f64; 3] {
    let (a, d) = (ra.to_radians(), dec.to_radians());
    [d.cos() * a.cos(), d.cos() * a.sin(), d.sin()]
}

/// ICRF (J2000 equator) to J2000 ecliptic.
fn ecliptic(v: [f64; 3]) -> [f64; 3] {
    let (s, c) = OBLIQUITY_J2000_DEG.to_radians().sin_cos();
    [v[0], c * v[1] + s * v[2], -s * v[1] + c * v[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn node<'a>(s: &'a SolarSystem, name: &str) -> &'a OrbitalNode {
        s.by_name(name).unwrap_or_else(|| panic!("{name}"))
    }

    #[test]
    fn earth_equinox_falls_on_the_published_moment() {
        // The table's Earth finds its own 2026 March equinox within the minutes its
        // approximations (barycentre, aberration, nutation) allow.
        let jd = northern_equinox_before("Earth", EQUINOX_2026_JD_TT + 1.0).unwrap();
        assert!((jd - EQUINOX_2026_JD_TT).abs() * 24.0 * 60.0 < 30.0, "{:.1} min off", (jd - EQUINOX_2026_JD_TT) * 1440.0);
    }

    #[test]
    fn mars_year_39_begins_in_september_2026_and_38_in_november_2024() {
        // Mars Year 38 began 2024-11-12 (Clancy et al.'s count; Ls 0): the last equinox
        // before March 2026.
        let jd = northern_equinox_before("Mars", EQUINOX_2026_JD_TT).unwrap();
        let nov_12_2024 = 2_460_626.5;
        assert!((jd - nov_12_2024).abs() < 2.0, "{:.1} days off", jd - nov_12_2024);
    }

    #[test]
    fn periods_and_planes_are_the_real_ones() {
        let s = at(Seed(1), EQUINOX_2026_JD_TT);
        let period = |n: &str| {
            let b = node(&s, n);
            orbital_period_days(b.elements.unwrap().a_au, s.mass_of(b.parent.unwrap()) + s.mass_of(b.key))
        };
        for (n, want) in [("Earth", 365.256), ("Jupiter", 4332.6), ("Moon", 27.3217), ("Io", 1.7627), ("Titan", 15.945), ("Triton", 5.8792)] {
            assert!((period(n) / want - 1.0).abs() < 2e-4, "{n}: {}", period(n));
        }
        let el = |n: &str| node(&s, n).elements.unwrap();
        // Uranus's moons orbit in its equator, 98° from its orbit; Triton retrograde.
        assert!((el("Titania").i_deg - 97.8).abs() < 1.0, "{}", el("Titania").i_deg);
        assert!(el("Triton").i_deg > 90.0);
        // The Galileans lie near Jupiter's equator, ~2° from the ecliptic.
        assert!(el("Ganymede").i_deg < 3.0);
        assert!((el("Moon").a_au * AU_KM - 384_400.0).abs() < 1000.0);
        assert_eq!(s.nodes.iter().filter(|n| n.parent == Some(node(&s, "Sun").key) && matches!(n.node, NodeKind::Body(_))).count(), 9);
    }

    #[test]
    fn archetypes_land_on_their_bodies() {
        for a in Archetype::ALL {
            let (s, body) = archetype(Seed(1), a);
            let n = s.node(body.key).unwrap();
            assert_eq!(n.name, a.real_body());
            assert_eq!(n.elements, Some(body.orbit));
            assert_eq!(n.node, NodeKind::Body(body.exterior));
        }
    }
}
