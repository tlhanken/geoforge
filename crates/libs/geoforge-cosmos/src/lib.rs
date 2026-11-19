//! Thin cosmology: one solar system as the world's root.
//!
//! The system is a pure function of its own key and seed (never of ancestors), so it can
//! later be attached to a galaxy sector without changing anything inside it. Its job is to
//! produce a [`BodyContext`] for the focus world; the surface pipeline never sees this crate.

use geoforge_core::body::{
    Archetype, BodyClass, BodyContext, BodyExterior, BodyInterior, Composition, Elements, Epoch, Irradiance, NodeKind, OrbitalNode,
    Phase, Population, PowerLaw, Shell, ShellKind, SolarSystem, Species, Star, SystemRoot, TectonicRegime, TidalPartner,
    EARTH_MASS_KG, EARTH_RADIUS_KM,
};
use geoforge_core::derived::{
    effective_temperature_k, flux_w_m2, frost_line_au, habitable_zone_au, main_sequence_luminosity, main_sequence_radius,
};
use geoforge_core::key::{Key, KeyAllocator};
use geoforge_core::seed::{range, unit, Seed};

/// Deterministic draws from a seed: each call uses a fresh labelled child.
struct Draw {
    seed: Seed,
    n: u64,
}

impl Draw {
    fn new(seed: Seed) -> Self {
        Self { seed, n: 0 }
    }
    fn u(&mut self) -> f64 {
        self.n += 1;
        unit(self.seed.child_u64(b"draw", self.n))
    }
    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        self.n += 1;
        range(self.seed.child_u64(b"draw", self.n), lo, hi)
    }
    fn normal(&mut self, mean: f64, sd: f64) -> f64 {
        let u1 = self.u().max(1e-12);
        let u2 = self.u();
        mean + sd * (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
}

/// A pronounceable name from a seed.
#[must_use]
pub fn name(seed: Seed) -> String {
    const ON: [&str; 16] = ["k", "t", "v", "s", "r", "m", "n", "l", "d", "th", "z", "b", "h", "sh", "g", "p"];
    const NU: [&str; 8] = ["a", "e", "i", "o", "u", "ae", "ai", "or"];
    let mut d = Draw::new(seed.child(b"name"));
    let syll = 2 + (d.u() * 2.0) as usize;
    let mut s = String::new();
    for _ in 0..syll {
        s.push_str(ON[(d.u() * ON.len() as f64) as usize % ON.len()]);
        s.push_str(NU[(d.u() * NU.len() as f64) as usize % NU.len()]);
    }
    let mut c = s.chars();
    c.next().map(|f| f.to_uppercase().collect::<String>() + c.as_str()).unwrap_or_default()
}

fn spectral_class(teff: f64) -> String {
    let (letter, lo, hi) = match teff {
        t if t >= 30_000.0 => ("O", 30_000.0, 50_000.0),
        t if t >= 10_000.0 => ("B", 10_000.0, 30_000.0),
        t if t >= 7_500.0 => ("A", 7_500.0, 10_000.0),
        t if t >= 6_000.0 => ("F", 6_000.0, 7_500.0),
        t if t >= 5_200.0 => ("G", 5_200.0, 6_000.0),
        t if t >= 3_700.0 => ("K", 3_700.0, 5_200.0),
        _ => ("M", 2_400.0, 3_700.0),
    };
    let sub = ((hi - teff) / (hi - lo) * 10.0).clamp(0.0, 9.0) as u32;
    format!("{letter}{sub}V")
}

const ROMAN: [&str; 12] = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII"];

/// A generated system and the key of its focus world.
pub struct Generated {
    pub system: SolarSystem,
    pub focus: Key,
}

/// Generate a solar system around a Sun-like star with an Earth-sized world in the
/// habitable zone (the focus), siblings inside and outside it, and belts in the gaps.
#[must_use]
pub fn generate_system(seed: Seed) -> Generated {
    let mut d = Draw::new(seed.child(b"system"));
    let mut keys = KeyAllocator::starting_at(0);
    let root_key = keys.fresh();
    let sys_name = name(seed);

    // Star: F8–K5 main sequence, so the habitable zone holds a stable, long-lived world.
    let m = d.range(0.70, 1.25);
    let l = main_sequence_luminosity(m);
    let r = main_sequence_radius(m);
    let teff = effective_temperature_k(l, r);
    let age = d.range(2.5, 7.5).min(10.0 * m.powf(-2.5) * 0.9);
    let star = Star { mass_solar: m, luminosity_solar: l, radius_solar: r, temperature_k: teff, age_gyr: age, metallicity_dex: d.normal(0.0, 0.18).clamp(-0.6, 0.4), spectral_class: spectral_class(teff) };
    let star_key = keys.fresh();
    let mut nodes = vec![OrbitalNode { key: star_key, name: sys_name.clone(), parent: None, node: NodeKind::Star(star), elements: None }];

    // Orbits: a geometric ladder through the habitable zone.
    // The focus sits where it receives 0.85–1.15× Earth's insolation: inside the habitable
    // zone, toward its warm edge, as Earth does.
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

    let mut focus = Key(0);
    let mut prev_a = 0.0;
    let mut first_giant_a = None;
    let mut last_rock_a: f64 = 0.0;
    for (i, &a) in orbits.iter().enumerate() {
        let key = keys.fresh();
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
        if matches!(class, BodyClass::GasGiant | BodyClass::IceGiant) && first_giant_a.is_none() {
            first_giant_a = Some(a);
        }
        if class == BodyClass::Terrestrial {
            last_rock_a = last_rock_a.max(a);
        }
        let radius_e = match class {
            BodyClass::Terrestrial | BodyClass::SmallBody => mass_e.powf(0.27),
            BodyClass::IceGiant => 3.9 * (mass_e / 17.0).powf(0.25),
            BodyClass::GasGiant => 11.2 * (mass_e / 318.0).powf(-0.04),
        };
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
        let nm = format!("{sys_name} {}", ROMAN.get(i).copied().unwrap_or("+"));
        nodes.push(OrbitalNode { key, name: nm, parent: Some(star_key), node: NodeKind::Body(ext), elements: Some(el) });
        if is_focus {
            focus = key;
        }
        prev_a = a;
    }
    let _ = prev_a;

    // Asteroid belt between the last rocky planet and the first giant, if the gap is wide.
    if let Some(g) = first_giant_a {
        if g / last_rock_a.max(1e-9) > 2.2 {
            let inner = last_rock_a * 1.35;
            let outer = g * 0.72;
            if outer > inner {
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
                });
            }
        }
    }
    // Outer icy belt beyond the last planet.
    let last = *orbits.last().unwrap_or(&30.0);
    nodes.push(OrbitalNode {
        key: keys.fresh(),
        name: format!("{sys_name} Outer Belt"),
        parent: Some(star_key),
        node: NodeKind::Belt(Population {
            inner_au: last * 1.3,
            outer_au: last * 1.9,
            count_estimate: 10f64.powf(d.range(8.0, 9.0)),
            size_law: PowerLaw { min_km: 1.0, max_km: d.range(1000.0, 2500.0), exponent: 3.0 },
            composition: Composition::from_weights(&[(Species::Ice, 0.6), (Species::Silicate, 0.4)]),
        }),
        elements: None,
    });

    // A large moon for the focus world (half the time), as a node of the tree.
    if d.u() < 0.55 {
        let fm = nodes.iter().find(|n| n.key == focus).and_then(|n| if let NodeKind::Body(b) = n.node { Some(b) } else { None });
        if let Some(fb) = fm {
            let moon_mass = fb.mass_kg * d.range(0.004, 0.02);
            let dist_r = d.range(30.0, 80.0);
            nodes.push(OrbitalNode {
                key: keys.fresh(),
                name: format!("{} a", nodes.iter().find(|n| n.key == focus).map_or(String::new(), |n| n.name.clone())),
                parent: Some(focus),
                node: NodeKind::Body(BodyExterior {
                    class: BodyClass::Terrestrial,
                    mass_kg: moon_mass,
                    radius_km: (moon_mass / EARTH_MASS_KG).powf(0.3) * EARTH_RADIUS_KM,
                    rotation_h: 0.0,
                    obliquity_deg: d.range(0.0, 7.0),
                }),
                elements: Some(Elements { a_au: dist_r * fb.radius_km / geoforge_core::body::AU_KM, e: d.range(0.0, 0.06), i_deg: d.range(0.0, 6.0), lan_deg: d.range(0.0, 360.0), argp_deg: d.range(0.0, 360.0), m0_deg: d.range(0.0, 360.0) }),
            });
        }
    }

    let root = SystemRoot { key: root_key, seed, name: sys_name, epoch: Epoch { system_age_gyr: age }, placement: None };
    Generated { system: SolarSystem { root, nodes }, focus }
}

/// The interior of a body (its own layer), from its exterior and orbit.
#[must_use]
pub fn interior_for(seed: Seed, ext: &BodyExterior, flux_w_m2: f64) -> BodyInterior {
    use Phase::{Liquid, Solid};
    use ShellKind::{Core, Crust, Mantle, Ocean};
    let mut d = Draw::new(seed.child(b"interior"));
    match ext.class {
        BodyClass::GasGiant => Archetype::GasGiant.interior(),
        BodyClass::IceGiant => Archetype::IceGiant.interior(),
        BodyClass::SmallBody => BodyInterior {
            shells: vec![
                Shell::new(Mantle, Solid, ext.radius_km * 0.8, Composition::single(Species::Silicate)),
                Shell::new(Crust, Solid, ext.radius_km * 0.2, Composition::from_weights(&[(Species::Silicate, 0.6), (Species::Ice, 0.4)])),
            ],
            regime: TectonicRegime::Inactive,
        },
        BodyClass::Terrestrial => {
            let r = ext.radius_km;
            let core = r * d.range(0.5, 0.58);
            let crust = d.range(25.0, 45.0);
            let temperate = (80.0..2500.0).contains(&flux_w_m2);
            let mut shells = vec![
                Shell::new(Core, Solid, core, Composition::from_weights(&[(Species::Iron, 0.85), (Species::Nickel, 0.05), (Species::Sulfur, 0.10)])),
                Shell::new(Mantle, Solid, r - core - crust, Composition::single(Species::Silicate)),
                Shell::new(Crust, Solid, crust, Composition::from_weights(&[(Species::Basalt, 0.6), (Species::Granite, 0.4)])),
            ];
            let wet = temperate && ext.mass_kg > 0.3 * EARTH_MASS_KG;
            if wet {
                shells.push(Shell::new(Ocean, Liquid, d.range(2.35, 2.9), Composition::single(Species::Water)));
                let air = Composition::from_moles(&[(Species::Nitrogen, 0.78), (Species::Oxygen, 0.21), (Species::Argon, 0.0093), (Species::CarbonDioxide, 0.0004)]);
                shells.push(Shell::atmosphere(d.range(0.7, 1.6), 100.0, air));
            }
            let regime = if wet { TectonicRegime::Plates } else { TectonicRegime::StagnantLid };
            BodyInterior { shells, regime }
        }
    }
}

/// The handoff to the surface pipeline for body `key` of `system`.
///
/// # Panics
/// If `key` is not a body orbiting a star.
#[must_use]
pub fn body_context(system: &SolarSystem, key: Key) -> BodyContext {
    let node = system.node(key).expect("body exists");
    let NodeKind::Body(ext) = node.node else { panic!("{key} is not a body") };
    let el = node.elements.expect("bodies have orbits");
    let star = system
        .nodes
        .iter()
        .find_map(|n| if let NodeKind::Star(s) = &n.node { Some(s.clone()) } else { None })
        .expect("system has a star");
    let flux = flux_w_m2(star.luminosity_solar, el.a_au);
    let tides = system
        .children(key)
        .iter()
        .filter_map(|m| {
            let NodeKind::Body(b) = m.node else { return None };
            Some(TidalPartner { name: m.name.clone(), mass_kg: b.mass_kg, distance_km: m.elements.map_or(0.0, |e| e.a_au * geoforge_core::body::AU_KM) })
        })
        .collect();
    BodyContext {
        key,
        name: node.name.clone(),
        exterior: ext,
        interior: interior_for(system.root.seed.child_u64(b"body", key.0), &ext, flux),
        orbit: el,
        irradiance: Irradiance { flux_w_m2: flux, star_temperature_k: star.temperature_k, eccentricity: el.e, obliquity_deg: ext.obliquity_deg },
        tides,
        epoch: system.root.epoch,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use geoforge_core::derived::{orbital_period_days, solar_masses};

    #[test]
    fn systems_are_plausible_and_deterministic() {
        for s in 0..40u64 {
            let g = generate_system(Seed(s));
            let g2 = generate_system(Seed(s));
            assert_eq!(g.system, g2.system);
            let ctx = body_context(&g.system, g.focus);
            // Focus world: Earth-ish flux, solid surface with plate tectonics and an ocean.
            assert!((0.8..1.25).contains(&(ctx.irradiance.flux_w_m2 / 1361.0)), "flux {}", ctx.irradiance.flux_w_m2);
            assert!(ctx.interior.has_lithosphere() && ctx.interior.surface_ocean_gel_km() > 1.0);
            let g_ms2 = ctx.gravity_m_s2();
            assert!((5.0..15.0).contains(&g_ms2), "gravity {g_ms2}");
            // Orbits are ordered and periods follow Kepler with the star's mass.
            let star = g.system.nodes.iter().find_map(|n| if let NodeKind::Star(s) = &n.node { Some(s.clone()) } else { None }).unwrap();
            let p = orbital_period_days(ctx.orbit.a_au, star.mass_solar * geoforge_core::body::SOLAR_MASS_KG);
            assert!(p > 50.0 && p < 1000.0, "period {p}");
            let _ = solar_masses(1.0);
        }
    }
}
