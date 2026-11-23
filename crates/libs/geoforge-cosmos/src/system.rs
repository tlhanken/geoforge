//! The orbital tree: a star, its planets on a geometric ladder through the habitable zone,
//! belts in the gaps, and moons.

use crate::draw::Draw;
use crate::dynamics::{self, HasOrbit};
use crate::planets::Formed;
use crate::stars::spectral_class;
use geoforge_core::body::{
    BodyClass, BodyExterior, Composition, Elements, Epoch, NodeKind, OrbitalNode, Population, PowerLaw, SolarSystem, Species, Star, StellarStage, SystemRoot, AU_KM, EARTH_MASS_KG, EARTH_RADIUS_KM,
    SOLAR_MASS_KG,
};
use geoforge_core::derived::{effective_temperature_k, frost_line_au, habitable_zone_au, main_sequence_luminosity, main_sequence_radius, orbital_period_days};
use geoforge_core::key::{Key, KeyAllocator};
use geoforge_core::seed::Seed;

/// A generated system and the key of its focus world.
#[derive(Debug, Clone, PartialEq)]
pub struct Generated {
    pub system: SolarSystem,
    pub focus: Key,
}

pub(crate) const ROMAN: [&str; 12] = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII"];

/// Radius of a solid body of `mass_kg` and mean density `density_g_cc`, km.
pub(crate) fn solid_radius_km(mass_kg: f64, density_g_cc: f64) -> f64 {
    (3.0 * mass_kg / (4.0 * std::f64::consts::PI * density_g_cc * 1000.0)).cbrt() / 1000.0
}

/// Generate a solar system around a Sun-like star with an Earth-sized world in the
/// habitable zone (the focus), siblings inside and outside it, belts in the gaps, moon
/// systems around the giants and, half the time, a large moon for the focus.
#[must_use]
pub fn generate_system(seed: Seed) -> Generated {
    let mut d = Draw::new(seed.child(b"system"));
    let mut keys = KeyAllocator::starting_at(0);
    let root_key = keys.fresh();
    let sys_name = geoforge_names::name(seed.0);

    // Star: F8–K5 main sequence, so the habitable zone holds a stable, long-lived world.
    let m = d.range(0.70, 1.25);
    let l = main_sequence_luminosity(m);
    let r = main_sequence_radius(m);
    let teff = effective_temperature_k(l, r);
    let age = d.range(2.5, 7.5).min(10.0 * m.powf(-2.5) * 0.9);
    let star = Star {
        mass_solar: m,
        luminosity_solar: l,
        radius_solar: r,
        temperature_k: teff,
        age_gyr: age,
        metallicity_dex: d.normal(0.0, 0.18).clamp(-0.6, 0.4),
        spectral_class: spectral_class(teff),
        stage: StellarStage::MainSequence,
        initial_mass_solar: None,
    };
    let star_key = keys.fresh();
    let mut nodes = vec![OrbitalNode { key: star_key, name: sys_name.clone(), parent: None, node: NodeKind::Star(star), elements: None, precession: None }];

    // Orbits: a geometric ladder through the habitable zone. The focus receives 0.85–1.15×
    // Earth's starlight: inside the habitable zone, toward its warm edge, as Earth does.
    let (hz_in, _hz_out) = habitable_zone_au(l);
    let focus_a = (l / d.range(0.85, 1.15)).sqrt().max(hz_in);
    let frost = frost_line_au(l);
    let mut orbits = vec![focus_a];
    let mut a = focus_a;
    while a > 0.25 * l.sqrt() + 0.05 {
        a /= d.range(1.45, 1.9);
        if a > 0.08 {
            orbits.push(a);
        }
    }
    a = focus_a;
    while a < 45.0 * l.sqrt() {
        a *= d.range(1.45, 2.1);
        orbits.push(a);
    }
    orbits.sort_by(f64::total_cmp);

    let mut drafts = Vec::new();
    for &a in &orbits {
        let is_focus = (a - focus_a).abs() < 1e-12;
        let (class, mass_e) = if is_focus {
            (BodyClass::Terrestrial, d.range(0.6, 1.6))
        } else if a < frost {
            (BodyClass::Terrestrial, d.range(0.05, 2.5))
        } else if a < frost * 2.5 && d.u() < 0.65 {
            (BodyClass::GasGiant, d.range(40.0, 700.0))
        } else if d.u() < 0.7 {
            (BodyClass::IceGiant, d.range(8.0, 25.0))
        } else {
            (BodyClass::SmallBody, d.range(0.001, 0.02))
        };
        let radius_e = planet_radius_earth(class, mass_e);
        let ext = BodyExterior {
            class,
            mass_kg: mass_e * EARTH_MASS_KG,
            radius_km: radius_e * EARTH_RADIUS_KM,
            rotation_h: match class {
                BodyClass::GasGiant | BodyClass::IceGiant => d.range(9.0, 18.0),
                _ => d.range(14.0, 40.0),
            },
            obliquity_deg: if is_focus { d.range(8.0, 32.0) } else { d.range(0.0, 35.0) },
        };
        let el = Elements {
            a_au: a,
            e: if is_focus { d.range(0.0, 0.04) } else { d.range(0.0, 0.12) },
            i_deg: d.range(0.0, 3.5),
            lan_deg: d.range(0.0, 360.0),
            argp_deg: d.range(0.0, 360.0),
            m0_deg: d.range(0.0, 360.0),
        };
        drafts.push(Draft { orbit: Formed { a_au: a, e: el.e, class, mass_earth: mass_e }, ext, el });
    }
    // Crowded neighbours move apart or scatter; the focus stays where it is.
    dynamics::settle(&mut Draw::new(seed.child(b"dynamics")), &mut drafts, m, Some(focus_a));

    let mut focus = Key(0);
    let mut planets = Vec::new();
    for (i, draft) in drafts.iter().enumerate() {
        let key = keys.fresh();
        let nm = format!("{sys_name} {}", ROMAN.get(i).copied().unwrap_or("+"));
        let el = Elements { a_au: draft.orbit.a_au, e: draft.orbit.e, ..draft.el };
        nodes.push(OrbitalNode { key, name: nm.clone(), parent: Some(star_key), node: NodeKind::Body(draft.ext), elements: Some(el), precession: None });
        planets.push((key, nm, draft.ext, draft.orbit.a_au));
        if (draft.orbit.a_au - focus_a).abs() < 1e-12 {
            focus = key;
        }
    }
    let giants: Vec<f64> = drafts.iter().filter(|p| matches!(p.orbit.class, BodyClass::GasGiant | BodyClass::IceGiant)).map(|p| p.orbit.a_au).collect();
    let last_rock = drafts.iter().filter(|p| p.orbit.class == BodyClass::Terrestrial && giants.first().is_none_or(|&g| p.orbit.a_au < g)).map(|p| p.orbit.a_au).fold(0.0, f64::max);

    // An asteroid belt inside the first giant, between its 4:1 and 2:1 resonances.
    if let Some((inner, outer)) = giants.first().and_then(|&g| dynamics::asteroid_belt(g, last_rock)) {
        nodes.push(OrbitalNode {
            key: keys.fresh(),
            name: format!("{sys_name} Belt"),
            parent: Some(star_key),
            node: NodeKind::Belt(Population {
                inner_au: inner,
                outer_au: outer,
                count_estimate: 10f64.powf(d.range(5.5, 6.5)),
                size_law: PowerLaw { min_km: 1.0, max_km: d.range(400.0, 950.0), exponent: 3.5 },
                composition: Composition::from_weights(&[(Species::Silicate, 0.75), (Species::Iron, 0.15), (Species::Ice, 0.10)]),
            }),
            elements: None,
            precession: None,
        });
    }
    // An outer icy belt: the outermost giant's 3:2 to 2:1 resonances, else beyond the planets.
    let last = drafts.last().map_or(30.0, |p| p.orbit.a_au);
    let (outer_in, outer_out) = giants.last().map_or((last * 1.3, last * 1.9), |&g| dynamics::kuiper_belt(g));
    nodes.push(OrbitalNode {
        key: keys.fresh(),
        name: format!("{sys_name} Outer Belt"),
        parent: Some(star_key),
        node: NodeKind::Belt(Population {
            inner_au: outer_in,
            outer_au: outer_out,
            count_estimate: 10f64.powf(d.range(8.0, 9.0)),
            size_law: PowerLaw { min_km: 1.0, max_km: d.range(1000.0, 2500.0), exponent: 3.0 },
            composition: Composition::from_weights(&[(Species::Ice, 0.6), (Species::Silicate, 0.4)]),
        }),
        elements: None,
        precession: None,
    });

    // Moons, inside each planet's Hill sphere.
    for (key, nm, ext, a_au) in &planets {
        let hill = hill_km(*a_au, ext.mass_kg, m);
        let moons = match ext.class {
            // Galilean-like: 2–5 large moons holding ~2×10⁻⁴ of the planet's mass.
            BodyClass::GasGiant => {
                let count = 2 + (d.u() * 4.0) as usize;
                moon_system(&mut d, &mut keys, *key, nm, ext, count, 2e-4)
            }
            // Uranus- or Neptune-like: none to two.
            BodyClass::IceGiant => {
                let count = (d.u() * 3.0) as usize;
                moon_system(&mut d, &mut keys, *key, nm, ext, count, 1e-4)
            }
            // The focus world: a large moon half the time.
            BodyClass::Terrestrial if *key == focus && d.u() < 0.55 => {
                let mass_kg = ext.mass_kg * d.range(0.004, 0.02);
                let drawn = d.range(30.0, 80.0) * ext.radius_km;
                big_moon_orbit_km(drawn, ext, hill)
                    .map(|a_km| moon(&mut d, keys.fresh(), *key, ext.mass_kg, MoonSpec { name: format!("{nm} a"), mass_kg, a_km, density_g_cc: 3.3 }))
                    .into_iter()
                    .collect()
            }
            _ => Vec::new(),
        };
        nodes.extend(moons.into_iter().filter(|n| n.elements.is_some_and(|e| e.a_au * AU_KM < MOON_HILL_SHARE * hill)));
    }

    let root = SystemRoot { key: root_key, seed, name: sys_name, epoch: Epoch { system_age_gyr: age }, placement: None, site: Some(sunlike_site(seed)), present_jd_tt: None };
    Generated { system: SolarSystem { root, nodes }, focus }
}

/// A place in the Milky Way like the Sun's: within a kiloparsec or so of its radius, at any
/// azimuth, in the old thin disk (an exponential of 250 pc in height).
fn sunlike_site(seed: Seed) -> geoforge_core::body::GalacticSite {
    let g = crate::galaxy::Galaxy::milky_way();
    let mut d = Draw::new(seed.child(b"site"));
    let r = (g.reference_kpc + d.normal(0.0, 0.7)).clamp(g.reference_kpc - 1.5, g.reference_kpc + 1.5) * 1000.0;
    let phi = std::f64::consts::TAU * d.u();
    let z = -250.0 * d.u().max(1e-12).ln() * if d.u() < 0.5 { 1.0 } else { -1.0 };
    crate::population::site(&g, [r * phi.cos(), r * phi.sin(), z], seed.child(b"galaxy"), None, &mut d)
}

/// A planet being placed: its orbit, settled with its neighbours before it gets a key.
struct Draft {
    orbit: Formed,
    ext: BodyExterior,
    el: Elements,
}

impl HasOrbit for Draft {
    fn orbit(&self) -> &Formed {
        &self.orbit
    }
    fn orbit_mut(&mut self) -> &mut Formed {
        &mut self.orbit
    }
}

/// Moons stay inside this share of their planet's Hill radius: prograde orbits are stable to
/// about half of it (Domingos et al. 2006).
pub(crate) const MOON_HILL_SHARE: f64 = 0.4;

/// Hill radius of a planet of `planet_kg` at `a_au` from a star of `star_solar` M☉, km.
pub(crate) fn hill_km(a_au: f64, planet_kg: f64, star_solar: f64) -> f64 {
    a_au * AU_KM * (planet_kg / (3.0 * star_solar * SOLAR_MASS_KG)).cbrt()
}

/// Fluid Roche limit of a moon of density `moon_g_cc` around `p`, km (the Moon's around
/// Earth: 18,000 km).
pub(crate) fn roche_km(p: &BodyExterior, moon_g_cc: f64) -> f64 {
    let planet_g_cc = p.mass_kg / (4.0 / 3.0 * std::f64::consts::PI * (p.radius_km * 1e5).powi(3)) * 1000.0;
    2.44 * p.radius_km * (planet_g_cc / moon_g_cc).cbrt()
}

/// Orbit of a big moon drawn at `drawn_km` from a rocky planet: brought inside the Hill
/// sphere if it lies beyond, `None` if that leaves it near the Roche limit.
pub(crate) fn big_moon_orbit_km(drawn_km: f64, p: &BodyExterior, hill_km: f64) -> Option<f64> {
    let a = drawn_km.min(0.9 * MOON_HILL_SHARE * hill_km);
    (a > 2.0 * roche_km(p, 3.3)).then_some(a)
}

/// Radius (Earth radii) of a planet of `class` and `mass_e` Earth masses.
pub(crate) fn planet_radius_earth(class: BodyClass, mass_e: f64) -> f64 {
    match class {
        BodyClass::Terrestrial | BodyClass::SmallBody => mass_e.powf(0.27),
        BodyClass::IceGiant => 3.9 * (mass_e / 17.0).powf(0.25),
        BodyClass::GasGiant => 11.2 * (mass_e / 318.0).powf(-0.04),
    }
}

/// `count` moons sharing `mass_share` of the planet's mass, on a widening ladder of orbits.
pub(crate) fn moon_system(d: &mut Draw, keys: &mut KeyAllocator, parent: Key, parent_name: &str, p: &BodyExterior, count: usize, mass_share: f64) -> Vec<OrbitalNode> {
    let weights: Vec<f64> = (0..count).map(|_| d.range(0.3, 1.0)).collect();
    let total: f64 = weights.iter().sum();
    let budget = p.mass_kg * mass_share * d.range(0.5, 1.5);
    let mut a_km = p.radius_km * d.range(5.5, 7.0);
    let letters = "abcdefgh";
    weights
        .iter()
        .enumerate()
        .map(|(i, w)| {
            let spec = MoonSpec { name: format!("{parent_name} {}", &letters[i..=i]), mass_kg: budget * w / total, a_km, density_g_cc: 1.9 };
            let node = moon(d, keys.fresh(), parent, p.mass_kg, spec);
            a_km *= d.range(1.5, 2.0);
            node
        })
        .collect()
}

/// A moon to place.
pub(crate) struct MoonSpec {
    pub name: String,
    pub mass_kg: f64,
    /// Orbit radius, km.
    pub a_km: f64,
    /// Bulk density: ~1.9 for ice-rich moons, ~3.3 for rocky ones.
    pub density_g_cc: f64,
}

/// A moon on a near-circular orbit, spinning once per orbit (tidally locked).
pub(crate) fn moon(d: &mut Draw, key: Key, parent: Key, parent_mass_kg: f64, spec: MoonSpec) -> OrbitalNode {
    let MoonSpec { name, mass_kg, a_km, density_g_cc } = spec;
    let a_au = a_km / AU_KM;
    let class = if mass_kg >= 5e-5 * EARTH_MASS_KG { BodyClass::Terrestrial } else { BodyClass::SmallBody };
    let ext = BodyExterior { class, mass_kg, radius_km: solid_radius_km(mass_kg, density_g_cc), rotation_h: orbital_period_days(a_au, parent_mass_kg) * 24.0, obliquity_deg: d.range(0.0, 2.0) };
    let el = Elements { a_au, e: d.range(0.0, 0.03), i_deg: d.range(0.0, 2.0), lan_deg: d.range(0.0, 360.0), argp_deg: d.range(0.0, 360.0), m0_deg: d.range(0.0, 360.0) };
    OrbitalNode { key, name, parent: Some(parent), node: NodeKind::Body(ext), elements: Some(el), precession: None }
}
