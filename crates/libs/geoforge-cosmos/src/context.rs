//! The handoff to the surface pipeline: a [`BodyContext`] for any body of a system, including
//! belt members generated on demand.

use crate::asteroids::{self, Asteroid};
use crate::draw::Draw;
use crate::worlds::{interior_for, regulate_co2, Setting};
use geoforge_core::body::{
    BodyClass, BodyContext, BodyExterior, Elements, Irradiance, NodeKind, OrbitalNode, Population, SolarSystem, Species, Star, StellarStage, TidalPartner, AU_KM, SOLAR_MASS_KG,
};
use geoforge_core::derived::{flux_w_m2, frost_line_au, main_sequence_luminosity, orbital_period_days};
use geoforge_core::key::Key;

/// A moon this close to its planet (in planet radii) is tidally heated.
const TIDAL_HEATING_RADII: f64 = 12.0;

/// What the top-level bodies under `parent` orbit: one star, or a barycentre's stars as one.
struct Host {
    luminosity_solar: f64,
    mass_solar: f64,
    /// Of the brightest star.
    temperature_k: f64,
    /// Where the frost line was when the planets formed, in today's orbits, AU.
    frost_au: f64,
    /// The year's mean light from the system's other stars, W/m² ([`others_flux_w_m2`]).
    others_flux_w_m2: f64,
}

/// The year's mean light that the stars a host does not include give its planets, W/m²: a wide
/// companion's on the primary's planets, the primary's on the companion's, an outer star's on
/// all. Each star's light is taken at its separation from the host, averaged over that orbit,
/// L/(4π a²√(1−e²)): the planets' own orbits are small beside it (they keep inside a third of
/// it to stay stable). Alpha Centauri A gives B's planets ~1/300 of Earth's sunlight on
/// average (1.5 L☉ at 23.4 AU, e 0.52).
fn others_flux_w_m2(system: &SolarSystem, host_key: Key, own: &[Key]) -> f64 {
    let Some(host) = system.node(host_key) else { return 0.0 };
    system
        .nodes
        .iter()
        .filter(|n| !own.contains(&n.key))
        .filter_map(|n| {
            let NodeKind::Star(s) = &n.node else { return None };
            let sep = if n.parent == Some(host_key) {
                n.elements?
            } else if host.parent == Some(n.key) {
                host.elements?
            } else {
                // Both orbit the same centre: the outer of the two orbits is the separation.
                [n.elements, host.elements].into_iter().flatten().max_by(|a, b| a.a_au.total_cmp(&b.a_au))?
            };
            Some(flux_w_m2(s.luminosity_solar, sep.a_au) / (1.0 - sep.e * sep.e).max(1e-6).sqrt())
        })
        .sum()
}

fn host(system: &SolarSystem, parent: Key) -> Option<Host> {
    let node = system.node(parent)?;
    let stars: Vec<(Key, &Star)> = match &node.node {
        NodeKind::Star(s) => vec![(parent, s)],
        NodeKind::Barycenter => system.children(parent).iter().filter_map(|c| if let NodeKind::Star(s) = &c.node { Some((c.key, s)) } else { None }).collect(),
        _ => return None,
    };
    let own: Vec<Key> = stars.iter().map(|(k, _)| *k).collect();
    let stars: Vec<&Star> = stars.into_iter().map(|(_, s)| s).collect();
    let brightest = stars.iter().max_by(|a, b| a.luminosity_solar.total_cmp(&b.luminosity_solar))?;
    let mass_solar: f64 = stars.iter().map(|s| s.mass_solar).sum();
    let luminosity_solar = stars.iter().map(|s| s.luminosity_solar).sum();
    // Stars without a recorded birth mass are as they were born; others formed their planets
    // in their main-sequence light, on orbits that have widened as they lost mass. A neutron
    // star's planets formed after the supernova, from rock.
    let frost_au = if stars.iter().all(|s| s.initial_mass_solar.is_none()) {
        frost_line_au(luminosity_solar)
    } else if stars.iter().any(|s| s.stage == StellarStage::NeutronStar) {
        f64::INFINITY
    } else {
        let birth: f64 = stars.iter().map(|s| s.initial_mass()).sum();
        frost_line_au(stars.iter().map(|s| main_sequence_luminosity(s.initial_mass())).sum()) * birth / mass_solar
    };
    Some(Host { luminosity_solar, mass_solar, temperature_k: brightest.temperature_k, frost_au, others_flux_w_m2: others_flux_w_m2(system, parent, &own) })
}

/// The ancestor of `key` that orbits the star directly (the body itself for a planet).
fn planet_level(system: &SolarSystem, key: Key) -> Option<&OrbitalNode> {
    let mut n = system.node(key)?;
    while let Some(p) = n.parent.and_then(|p| system.node(p)) {
        if matches!(p.node, NodeKind::Star(_) | NodeKind::Barycenter) {
            return Some(n);
        }
        n = p;
    }
    None
}

/// The context of body `key`, or `None` if it is not a body orbiting in this system. A body of
/// the real Solar System that has an archetype is that archetype, measured, not drawn.
#[must_use]
pub fn body_context(system: &SolarSystem, key: Key) -> Option<BodyContext> {
    if let Some(known) = crate::sol::known_context(system, key) {
        return Some(known);
    }
    let node = system.node(key)?;
    let NodeKind::Body(ext) = node.node else { return None };
    let orbit = node.elements?;
    let planet = planet_level(system, key)?;
    let star = host(system, planet.parent?)?;
    let year = planet.elements?;
    let parent = node.parent.and_then(|p| system.node(p)).and_then(|p| if let NodeKind::Body(b) = p.node { Some((p, b)) } else { None });

    let flux = flux_w_m2(star.luminosity_solar, year.a_au) + star.others_flux_w_m2;
    // A moon of a giant spins in its planet's equatorial plane, so it shares the planet's seasons.
    let obliquity_deg = match parent {
        Some((_, pb)) if matches!(pb.class, BodyClass::GasGiant | BodyClass::IceGiant) => pb.obliquity_deg,
        _ => ext.obliquity_deg,
    };
    let irradiance =
        Irradiance { flux_w_m2: flux, star_temperature_k: star.temperature_k, eccentricity: year.e, obliquity_deg, year_days: orbital_period_days(year.a_au, star.mass_solar * SOLAR_MASS_KG) };
    let setting = Setting { flux_w_m2: flux, icy: year.a_au > star.frost_au, tidally_heated: parent.is_some_and(|(_, pb)| orbit.a_au * AU_KM < TIDAL_HEATING_RADII * pb.radius_km) };
    // Tides: a moon feels its planet, a planet its moons; both feel the star.
    let star_tide = TidalPartner { name: system.node(planet.parent?).map_or_else(String::new, |n| n.name.clone()), mass_kg: star.mass_solar * SOLAR_MASS_KG, distance_km: year.a_au * AU_KM };
    let mut tides: Vec<TidalPartner> = match parent {
        Some((p, pb)) => vec![TidalPartner { name: p.name.clone(), mass_kg: pb.mass_kg, distance_km: orbit.a_au * AU_KM }],
        None => system
            .children(key)
            .iter()
            .filter_map(|m| {
                let NodeKind::Body(b) = m.node else { return None };
                Some(TidalPartner { name: m.name.clone(), mass_kg: b.mass_kg, distance_km: m.elements.map_or(0.0, |e| e.a_au * AU_KM) })
            })
            .collect(),
    };
    tides.push(star_tide);
    let mut ctx = BodyContext {
        key,
        name: node.name.clone(),
        exterior: ext,
        interior: interior_for(system.root.seed.child_u64(b"body", key.0), &ext, &setting),
        orbit,
        irradiance,
        tides,
        epoch: system.root.epoch,
    };
    regulate_co2(&mut ctx);
    Some(ctx)
}

/// The `index`-th largest member of belt `belt` (0 = the largest), generated on demand: belts
/// are stored as populations, never as lists of bodies.
#[must_use]
pub fn belt_member(system: &SolarSystem, belt: Key, index: u64) -> Option<BodyContext> {
    belt_member_made(system, belt, index).map(|(ctx, _)| ctx)
}

/// A belt member with what it is made of (see [`crate::asteroids`]).
#[must_use]
pub fn belt_member_made(system: &SolarSystem, belt: Key, index: u64) -> Option<(BodyContext, Asteroid)> {
    let node = system.node(belt)?;
    let NodeKind::Belt(pop) = &node.node else { return None };
    let star = host(system, node.parent?)?;
    let key = Key((belt.0 << 32) | (index + 1));
    let mut d = Draw::new(system.root.seed.child_u64(b"belt-member", key.0));
    let diameter_km = member_diameter_km(pop, index);
    let icy = pop.composition.fraction(Species::Ice) > 0.3;
    let rotation_h = d.range(4.0, 20.0);
    let obliquity_deg = d.range(0.0, 60.0);
    let mut orbit = Elements {
        a_au: d.range(pop.inner_au, pop.outer_au),
        e: d.range(0.0, 0.2),
        i_deg: d.range(0.0, 15.0),
        lan_deg: d.range(0.0, 360.0),
        argp_deg: d.range(0.0, 360.0),
        m0_deg: d.range(0.0, 360.0),
    };
    let made = asteroids::member(system.root.seed.child_u64(b"belt", belt.0), pop, index, diameter_km, star.frost_au, &mut orbit);
    let radius_km = diameter_km / 2.0;
    let mass_kg = 4.0 / 3.0 * std::f64::consts::PI * (radius_km * 1000.0).powi(3) * made.density_g_cc * 1000.0;
    // Round (dwarf planets) from about 800 km across.
    let class = if diameter_km >= 800.0 { BodyClass::Terrestrial } else { BodyClass::SmallBody };
    let ext = BodyExterior { class, mass_kg, radius_km, rotation_h, obliquity_deg };
    let a_au = orbit.a_au;
    let flux = flux_w_m2(star.luminosity_solar, a_au) + star.others_flux_w_m2;
    let setting = Setting { flux_w_m2: flux, icy: icy || a_au > star.frost_au, tidally_heated: false };
    let ctx = BodyContext {
        key,
        name: format!("{} #{}", node.name, index + 1),
        exterior: ext,
        interior: interior_for(system.root.seed.child_u64(b"body", key.0), &ext, &setting),
        orbit,
        irradiance: Irradiance {
            flux_w_m2: flux,
            star_temperature_k: star.temperature_k,
            eccentricity: orbit.e,
            obliquity_deg: ext.obliquity_deg,
            year_days: orbital_period_days(a_au, star.mass_solar * SOLAR_MASS_KG),
        },
        tides: Vec::new(),
        epoch: system.root.epoch,
    };
    Some((ctx, made))
}

/// Diameter of the `index`-th largest member: from N(>D) = count · (D/min)^(1−q).
fn member_diameter_km(pop: &Population, index: u64) -> f64 {
    let law = pop.size_law;
    let rank = index as f64 + 1.0;
    (law.min_km * (pop.count_estimate / rank).powf(1.0 / (law.exponent - 1.0))).min(law.max_km)
}

/// Any body by key: a node of the tree, or a belt member (as [`belt_member`] keys them).
#[must_use]
pub fn any_body(system: &SolarSystem, key: Key) -> Option<BodyContext> {
    match key.0 >> 32 {
        0 => body_context(system, key),
        belt => belt_member(system, Key(belt), (key.0 & 0xffff_ffff).checked_sub(1)?),
    }
}

/// A body by name: a node of the tree, or a belt member as `"<belt name> #<n>"` (1 = largest).
#[must_use]
pub fn find_body(system: &SolarSystem, name: &str) -> Option<BodyContext> {
    if let Some(n) = system.by_name(name) {
        return body_context(system, n.key);
    }
    let (belt, n) = name.rsplit_once('#')?;
    let index: u64 = n.trim().parse().ok().filter(|&i| i >= 1)?;
    belt_member(system, system.by_name(belt.trim())?.key, index - 1)
}
