//! What is shown at each level, as plain data for the page to draw.
//!
//! - **Universe**: the galaxies around the home galaxy (the focus body's, or the default),
//!   heavier than 10⁸ M☉ within 10 Mpc, by position, type and size; the real universe's are
//!   the Local Group and the bright galaxies near it.
//! - **Galaxy**: its shape as a cloud of stars drawn from its model (colour by age), and the
//!   systems made in it.
//! - **Stars**: the stars around one, within 15 pc, by colour and brightness.
//! - **System**: its stars, bodies (position on the day, orbit as a ring of points) and
//!   belts.
//! - **Body**: what it is, whether it is made and generated, and the maps it can show.

use crate::{star_id, Result, ViewError, Viewer};
use geoforge_analysis::appearance::star_colour;
use geoforge_core::address::{Address, GalaxyRef, Level, StarRef};
use geoforge_core::body::{BodyClass, NodeKind, SolarSystem, StellarStage};
use geoforge_core::galactic::UniverseSite;
use geoforge_cosmos::galaxy::Morphology;
use geoforge_cosmos::remnants::notes;
use geoforge_cosmos::universe::Universe;
use geoforge_cosmos::worlds::{climate_estimate, habitability, Habitability};
use serde::Serialize;

/// The universe around its home galaxy, out to this far, Mpc.
const UNIVERSE_RADIUS_MPC: f64 = 10.0;
/// Galaxies lighter than this are left out of the universe's view, M☉ of stars.
const UNIVERSE_LIGHTEST_MSUN: f64 = 1e8;
/// Stars drawn to show a galaxy's shape: this many by mass (its old body), and half as many
/// young ones (under `YOUNG_GYR`), which hold much of its light and trace its arms.
const CLOUD_STARS: usize = 20_000;
const YOUNG_GYR: f64 = 0.1;
/// The stars shown around one, out to this far, pc.
const STARS_RADIUS_PC: f64 = 15.0;
/// Points on a drawn orbit.
const ORBIT_POINTS: usize = 96;

/// A scene at one level.
#[derive(Debug, Serialize)]
#[serde(tag = "level", rename_all = "snake_case")]
pub enum Scene {
    Universe(UniverseScene),
    Galaxy(GalaxyScene),
    Stars(StarsScene),
    System(SystemScene),
    Body(BodyScene),
}

#[derive(Debug, Serialize)]
pub struct UniverseScene {
    pub address: String,
    pub name: String,
    pub radius_mpc: f64,
    pub galaxies: Vec<GalaxyDot>,
}

#[derive(Debug, Serialize)]
pub struct GalaxyDot {
    pub address: Option<String>,
    pub name: String,
    pub kind: String,
    /// From the home galaxy, Mpc.
    pub position_mpc: [f32; 3],
    pub stellar_msun: f64,
    pub radius_kpc: f64,
    pub colour: [u8; 3],
    pub home: bool,
    /// What its active nucleus is, if it is active ("Seyfert nucleus", "quasar").
    pub active: Option<String>,
}

#[derive(Debug, Serialize)]
pub struct GalaxyScene {
    pub address: String,
    pub name: String,
    pub kind: String,
    pub stellar_msun: f64,
    /// Where its disks end, kpc.
    pub radius_kpc: f64,
    /// Its shape: x, y, z (kpc) and r, g, b (0–1) of each drawn star, one after another;
    /// its old stars and, apart, its young ones (which trace the arms).
    pub cloud: Vec<f32>,
    pub young: Vec<f32>,
    /// The systems made in it, and the star it centres on.
    pub systems: Vec<SystemDot>,
    /// What sits at its centre ([`crate::nebulae::centre_line`]).
    pub centre: String,
    pub globulars: Vec<crate::nebulae::DeepSkyDot>,
}

#[derive(Debug, Serialize)]
pub struct SystemDot {
    pub address: String,
    pub name: String,
    pub position_kpc: [f32; 3],
    pub made: bool,
}

#[derive(Debug, Serialize)]
pub struct StarsScene {
    pub address: String,
    pub galaxy: String,
    pub radius_pc: f64,
    pub stars: Vec<StarDot>,
}

#[derive(Debug, Serialize)]
pub struct StarDot {
    pub address: String,
    pub name: String,
    pub spectral: String,
    /// "star", "white_dwarf", "brown_dwarf", "neutron_star" or "black_hole".
    pub kind: &'static str,
    /// From the centre star, pc.
    pub position_pc: [f32; 3],
    pub colour: [u8; 3],
    /// `None` for a dark remnant.
    pub absolute_v: Option<f64>,
    pub made: bool,
    pub centre: bool,
}

/// A star's kind as the page draws it.
fn kind_of(stage: StellarStage) -> &'static str {
    match stage {
        StellarStage::WhiteDwarf => "white_dwarf",
        StellarStage::BrownDwarf => "brown_dwarf",
        StellarStage::NeutronStar => "neutron_star",
        StellarStage::BlackHole => "black_hole",
        _ => "star",
    }
}

#[derive(Debug, Serialize)]
pub struct SystemScene {
    pub address: String,
    pub name: String,
    pub galaxy: String,
    pub age_gyr: f64,
    pub stars: Vec<StarBody>,
    pub bodies: Vec<BodyDot>,
    pub belts: Vec<Belt>,
    pub rings: Vec<Ring>,
    /// Where its fictional materials are: material, place, and the body to go to.
    pub exotics: Vec<ExoticFind>,
}

/// A fictional material found in a system.
#[derive(Debug, Serialize)]
pub struct ExoticFind {
    pub material: String,
    pub place: String,
    pub address: Option<String>,
}

/// A fictional material's deposit on a body's surface.
#[derive(Debug, Serialize)]
pub struct Marker {
    pub material: String,
    pub colour: [u8; 3],
    pub lat: f64,
    pub lon: f64,
}

/// A planet's rings, in its radii, tilted with its equator.
#[derive(Debug, Clone, Serialize)]
pub struct Ring {
    pub parent: u64,
    pub inner_r: f64,
    pub outer_r: f64,
    /// Normal optical depth: ~0.4 for Saturn's, ~10⁻⁶ for Jupiter's dust.
    pub depth: f64,
    /// Ice (bright) rather than dust.
    pub icy: bool,
    pub tilt_deg: f64,
}

/// The rings of `system`'s node `parent`, if any.
fn ring_of(system: &SolarSystem, parent: geoforge_core::Key) -> Option<Ring> {
    let NodeKind::Body(planet) = &system.node(parent)?.node else { return None };
    system.nodes.iter().filter(|n| n.parent == Some(parent)).find_map(|n| match &n.node {
        NodeKind::Ring(p) => {
            let r_au = planet.radius_km / geoforge_core::body::AU_KM;
            Some(Ring {
                parent: parent.0,
                inner_r: p.inner_au / r_au,
                outer_r: p.outer_au / r_au,
                depth: p.optical_depth(),
                icy: p.composition.fraction(geoforge_core::body::Species::Ice) > 0.5,
                tilt_deg: planet.obliquity_deg,
            })
        }
        _ => None,
    })
}

#[derive(Debug, Serialize)]
pub struct StarBody {
    pub key: u64,
    pub name: String,
    pub spectral: String,
    pub kind: &'static str,
    pub colour: [u8; 3],
    pub radius_solar: f64,
    pub temperature_k: f64,
    pub position_au: [f32; 3],
    /// A remnant's spin, field, horizon or feeding ([`notes`]).
    pub notes: Vec<String>,
}

#[derive(Debug, Serialize)]
pub struct BodyDot {
    pub key: u64,
    pub address: String,
    pub name: String,
    pub class: String,
    pub parent: Option<u64>,
    pub radius_km: f64,
    pub position_au: [f32; 3],
    /// Its orbit around its parent, AU from the parent.
    pub orbit_au: Vec<[f32; 3]>,
    pub made: bool,
    pub focus: bool,
    /// "temperate" or "habitable" for a world with open water ([`habitability`]).
    pub habitable: Option<&'static str>,
}

#[derive(Debug, Serialize)]
pub struct Belt {
    pub name: String,
    pub parent: Option<u64>,
    pub inner_au: f64,
    pub outer_au: f64,
    /// Its members' types, "C 45%, S 30%, …".
    pub mix: Option<String>,
}

#[derive(Debug, Serialize)]
pub struct BodyScene {
    pub address: String,
    pub name: String,
    pub class: String,
    pub system: String,
    pub made: bool,
    /// The maps it can show now (true colour for every body).
    pub layers: Vec<String>,
    pub facts: Vec<[String; 2]>,
    /// Generation under way or done, if any.
    pub job: Option<String>,
    pub ring: Option<Ring>,
    /// Its fictional materials' deposits.
    pub deposits: Vec<Marker>,
}

/// One step of the way down to a scene: its address, name, and whether it is the stars
/// around a system rather than the system.
#[derive(Debug, Serialize)]
pub struct Crumb {
    pub address: String,
    pub name: String,
    pub stars: bool,
}

impl Viewer {
    /// The way down to `address` from its universe, named.
    ///
    /// # Errors
    /// If a level names nothing.
    pub fn crumbs(&self, address: &Address) -> Result<Vec<Crumb>> {
        let crumb = |a: Address, name: String, stars: bool| Crumb { address: a.to_string(), name, stars };
        let mut out = vec![crumb(address.up_to(Level::Universe), format!("Universe {}", address.up_to(Level::Universe)), false)];
        if address.level() >= Level::Galaxy {
            out.push(crumb(address.up_to(Level::Galaxy), self.galaxy(address)?.at.model.name.clone(), false));
        }
        if address.level() >= Level::System {
            let sys = self.system(address)?.system;
            out.push(crumb(address.up_to(Level::System), "Neighbourhood".into(), true));
            out.push(crumb(address.up_to(Level::System), sys.root.name.clone(), false));
            if let Some(key) = address.body {
                out.push(crumb(*address, sys.node(key).map_or_else(|| key.to_string(), |n| n.name.clone()), false));
            }
        }
        Ok(out)
    }

    /// The scene at `address`; with `stars`, the stars around the system at `address`.
    ///
    /// # Errors
    /// If the address names nothing.
    pub fn scene(&self, address: &Address, stars: bool) -> Result<Scene> {
        Ok(match (address.level(), stars) {
            (Level::Universe, _) => Scene::Universe(self.universe_scene()?),
            (Level::Galaxy, _) => Scene::Galaxy(self.galaxy_scene(address)?),
            (Level::System, true) | (Level::Body, true) => Scene::Stars(self.stars_scene(address)?),
            (Level::System, false) => Scene::System(self.system_scene(address)?),
            (Level::Body, false) => Scene::Body(self.body_scene(address)?),
        })
    }

    fn universe_scene(&self) -> Result<UniverseScene> {
        let (universe, _) = self.start();
        let home = self.home();
        let galaxies = match universe.universe {
            UniverseSite::Real => real_galaxies(&universe),
            UniverseSite::Seeded(seed) => {
                let u = Universe::new(seed);
                let Some(GalaxyRef::Generated(id)) = home.galaxy else { return Err(ViewError("no home galaxy".into())) };
                let h = u.galaxy(&id).ok_or_else(|| ViewError(format!("no galaxy {id}")))?;
                u.galaxies_seen(h.position_mpc, UNIVERSE_RADIUS_MPC, 15.5)
                    .into_iter()
                    .filter(|g| g.stellar_msun >= UNIVERSE_LIGHTEST_MSUN || g.id == h.id)
                    .map(|g| GalaxyDot {
                        address: Some(universe.galaxy(GalaxyRef::Generated(g.id)).to_string()),
                        name: g.name(),
                        kind: g.kind(),
                        position_mpc: [0, 1, 2].map(|k| (g.position_mpc[k] - h.position_mpc[k]) as f32),
                        stellar_msun: g.stellar_msun,
                        radius_kpc: g.half_light_kpc * 2.0,
                        colour: type_colour(g.morphology),
                        home: g.id == h.id,
                        active: active(&g.nucleus()),
                    })
                    .collect()
            }
        };
        Ok(UniverseScene { address: universe.to_string(), name: format!("Universe {universe}"), radius_mpc: UNIVERSE_RADIUS_MPC, galaxies })
    }

    fn galaxy_scene(&self, address: &Address) -> Result<GalaxyScene> {
        let g = self.galaxy(address)?;
        let model = &g.at.model;
        let (cloud, young) = g
            .cloud
            .get_or_init(|| {
                let flat = |stars: Vec<([f64; 3], f64)>| -> Vec<f32> {
                    stars
                        .into_iter()
                        .flat_map(|(p, age)| {
                            let c = age_colour(age);
                            [p[0] / 1000.0, p[1] / 1000.0, p[2] / 1000.0, c[0], c[1], c[2]].map(|v| v as f32)
                        })
                        .collect()
                };
                let seed = g.at.stars_seed.child(b"view");
                (flat(model.sample(seed, CLOUD_STARS)), flat(model.sample_young(seed, CLOUD_STARS / 2, YOUNG_GYR)))
            })
            .clone();
        let galaxy = address.up_to(Level::Galaxy);
        let mut systems: Vec<SystemDot> = self
            .worlds()?
            .into_iter()
            .filter_map(|w| {
                let a = w.manifest.address?;
                (a.up_to(Level::Galaxy) == galaxy).then_some(())?;
                let site = w.load_system().ok()??.root.site?;
                Some(SystemDot {
                    address: a.up_to(Level::System).to_string(),
                    name: w.manifest.focus_name.clone(),
                    position_kpc: site.position_pc.unwrap_or_default().map(|v| (v / 1000.0) as f32),
                    made: true,
                })
            })
            .collect();
        if systems.is_empty() {
            if let Some(star) = self.home_star(&galaxy).and_then(|s| star_id(&s).and_then(|id| g.field.star(&id))) {
                systems.push(SystemDot {
                    address: galaxy.star(StarRef::Star(star.id)).to_string(),
                    name: g.field.name_of(&star.id),
                    position_kpc: star.position_pc.map(|v| (v / 1000.0) as f32),
                    made: false,
                });
            }
        }
        let home = systems.first().map_or([0.0; 3], |s| s.position_kpc.map(|v| f64::from(v) * 1000.0));
        let globulars = crate::nebulae::globular_dots(&g.at, home);
        Ok(GalaxyScene {
            address: galaxy.to_string(),
            name: model.name.clone(),
            kind: model.morphology.label(),
            stellar_msun: model.mass_msun(),
            radius_kpc: model.edge_kpc,
            cloud,
            young,
            systems,
            centre: crate::nebulae::centre_line(&g.at.nucleus()),
            globulars,
        })
    }

    fn stars_scene(&self, address: &Address) -> Result<StarsScene> {
        let g = self.galaxy(address)?;
        let galaxy = address.up_to(Level::Galaxy);
        let system = address.up_to(Level::System);
        let centre = match system.star {
            Some(StarRef::Star(id)) => g.field.star(&id).ok_or_else(|| ViewError(format!("no star {id}")))?.position_pc,
            _ => self.system(&system)?.system.root.site.and_then(|s| s.position_pc).ok_or_else(|| ViewError("the system has no place".into()))?,
        };
        let made: Vec<Address> = self.worlds()?.into_iter().filter_map(|w| w.manifest.address.map(|a| a.up_to(Level::System))).collect();
        let mut stars: Vec<StarDot> = g
            .field
            .near(centre, STARS_RADIUS_PC, f64::INFINITY)
            .into_iter()
            .map(|s| {
                let a = galaxy.star(StarRef::Star(s.id));
                StarDot {
                    address: a.to_string(),
                    name: g.field.name_of(&s.id),
                    spectral: s.star.spectral_class.clone(),
                    kind: kind_of(s.star.stage),
                    position_pc: [0, 1, 2].map(|k| (s.position_pc[k] - centre[k]) as f32),
                    colour: star_colour(&s.star),
                    absolute_v: s.absolute_v.is_finite().then_some(s.absolute_v),
                    made: made.contains(&a),
                    centre: a == system,
                }
            })
            .collect();
        if system.star == Some(StarRef::Sun) {
            let sys = self.system(&system)?.system;
            if let Some(sun) = sys.nodes.iter().find_map(|n| if let NodeKind::Star(s) = &n.node { Some(s.clone()) } else { None }) {
                stars.push(StarDot {
                    address: system.to_string(),
                    name: "Sun".into(),
                    spectral: sun.spectral_class.clone(),
                    kind: "star",
                    position_pc: [0.0; 3],
                    colour: star_colour(&sun),
                    absolute_v: Some(4.83),
                    made: made.contains(&system),
                    centre: true,
                });
            }
        }
        Ok(StarsScene { address: system.to_string(), galaxy: galaxy.to_string(), radius_pc: STARS_RADIUS_PC, stars })
    }

    fn system_scene(&self, address: &Address) -> Result<SystemScene> {
        let at = address.up_to(Level::System);
        let generated = self.system(&at)?;
        let sys = &generated.system;
        let placed = geoforge_almanac::ephemeris::positions(sys, 0.0);
        let where_is = |k| placed.iter().find(|p| p.key == k).map_or([0.0; 3], |p| p.position_au.map(|v| v as f32));
        let made: Vec<Address> = self.worlds()?.into_iter().filter_map(|w| w.manifest.address).collect();
        let mut stars = Vec::new();
        let mut bodies = Vec::new();
        let mut belts = Vec::new();
        for n in &sys.nodes {
            match &n.node {
                NodeKind::Star(s) => stars.push(StarBody {
                    key: n.key.0,
                    name: n.name.clone(),
                    spectral: s.spectral_class.clone(),
                    kind: kind_of(s.stage),
                    colour: star_colour(s),
                    radius_solar: s.radius_solar,
                    temperature_k: s.temperature_k,
                    position_au: where_is(n.key),
                    notes: notes(sys, n.key),
                }),
                NodeKind::Body(b) => {
                    let a = at.body(n.key);
                    let habitable = body_of(&a, sys).ok().and_then(|c| habitable_name(habitability(&c)));
                    bodies.push(BodyDot {
                        key: n.key.0,
                        address: a.to_string(),
                        name: n.name.clone(),
                        class: class_name(b.class).into(),
                        parent: n.parent.map(|p| p.0),
                        radius_km: b.radius_km,
                        position_au: where_is(n.key),
                        orbit_au: n.elements.map(|e| (0..ORBIT_POINTS).map(|k| e.position_at(360.0 * k as f64 / ORBIT_POINTS as f64).map(|v| v as f32)).collect()).unwrap_or_default(),
                        made: made.contains(&a),
                        focus: Some(n.key) == generated.focus,
                        habitable,
                    });
                }
                NodeKind::Belt(p) => {
                    belts.push(Belt { name: n.name.clone(), parent: n.parent.map(|k| k.0), inner_au: p.inner_au, outer_au: p.outer_au, mix: geoforge_cosmos::asteroids::mix_line(sys, n.key) })
                }
                _ => {}
            }
        }
        let rings = sys.nodes.iter().filter(|n| matches!(n.node, NodeKind::Ring(_))).filter_map(|n| n.parent.and_then(|p| ring_of(sys, p))).collect();
        let exotics = geoforge_cosmos::exotics::in_system(&self.exotics()?, sys, true)
            .into_iter()
            .map(|f| ExoticFind { material: f.material, place: f.place, address: f.body.map(|k| at.up_to(Level::System).body(k).to_string()) })
            .collect();
        Ok(SystemScene {
            address: at.to_string(),
            name: sys.root.name.clone(),
            galaxy: at.up_to(Level::Galaxy).to_string(),
            age_gyr: sys.root.epoch.system_age_gyr,
            stars,
            bodies,
            belts,
            rings,
            exotics,
        })
    }

    fn body_scene(&self, address: &Address) -> Result<BodyScene> {
        let world = self.world(address)?;
        let sys = self.system(address)?.system;
        let body = match &world {
            Some(w) => w.load_body()?,
            None => body_of(address, &sys)?,
        };
        let mut layers = vec!["appearance".to_string()];
        if let Some(w) = &world {
            layers.extend(geoforge_pipeline::map_names(w).into_iter().filter(|n| *n != "appearance").map(str::to_string));
        }
        let i = &body.interior;
        let facts = [
            ("Kind", class_name(body.exterior.class).to_string()),
            ("Radius", format!("{:.0} km", body.exterior.radius_km)),
            ("Gravity", format!("{:.2} g", body.gravity_m_s2() / 9.80665)),
            ("Day", format!("{:.1} h", body.exterior.rotation_h)),
            ("Year", format!("{:.1} days", body.irradiance.year_days)),
            ("Starlight", format!("{:.2}× Earth's", body.irradiance.flux_w_m2 / 1361.0)),
            ("Air", if i.has_atmosphere() { format!("{:.2} bar", i.surface_pressure_bar()) } else { "none".into() }),
            ("Water", if i.has_surface_liquid_water() { format!("seas, {:.1} km deep if spread", i.surface_ocean_gel_km()) } else { "no seas".into() }),
        ]
        .map(|(k, v)| [k.to_string(), v])
        .to_vec();
        let mut facts = facts;
        if let Some(a) = geoforge_cosmos::asteroids::of(&sys, body.key) {
            facts.extend(a.facts(body.exterior.mass_kg).into_iter().map(|(k, v)| [k.to_string(), v]));
        }
        if let Some(c) = climate_estimate(&body) {
            let h = habitable_name(habitability(&body)).unwrap_or("no open water");
            facts.push(["Climate".into(), format!("~{:.0} °C, {:.0}% thawed: {h}", c.mean_c, (100.0 * c.thawed_share).max(0.0))]);
        }
        let found = geoforge_cosmos::exotics::on_body(&self.exotics()?, &sys, body.key);
        let mut deposits = Vec::new();
        for d in found {
            let n = d.sites.len();
            let note = if d.note.is_empty() { String::new() } else { format!(" — {}", d.note) };
            facts.push([format!("{} (fictional)", d.material), format!("{n} deposit{}, ~{:.0e} t{note}", if n == 1 { "" } else { "s" }, d.tonnes)]);
            deposits.extend(d.sites.iter().map(|&(lat, lon)| Marker { material: d.material.clone(), colour: d.colour, lat, lon }));
        }
        Ok(BodyScene {
            address: address.to_string(),
            name: body.name.clone(),
            class: class_name(body.exterior.class).into(),
            system: address.up_to(Level::System).to_string(),
            made: world.is_some(),
            layers,
            facts,
            job: self.jobs.lock().expect("jobs").get(address).cloned(),
            ring: address.body.and_then(|k| ring_of(&sys, k)),
            deposits,
        })
    }

    /// Body `address`'s map `layer` as a PNG `width` pixels wide: any map of a generated
    /// world; true colour of any body (drawn from what it is, if not generated).
    ///
    /// # Errors
    /// If the body or map does not exist.
    pub fn texture(&self, address: &Address, layer: &str, width: u32) -> Result<Vec<u8>> {
        let width = width.clamp(64, 4096);
        let img = match self.world(address)? {
            Some(w) => geoforge_pipeline::map_image(&w, layer)?,
            None => None,
        };
        let img = match img {
            Some(i) => image::imageops::resize(&i, width, width / 2, image::imageops::FilterType::Triangle),
            None if layer == "appearance" => {
                let sys = self.system(address)?.system;
                let body = body_of(address, &sys)?;
                let pixels = geoforge_analysis::appearance::body(&body, sys.root.seed.child_u64(b"look", body.key.0), width as usize);
                image::RgbImage::from_fn(width, width / 2, |x, y| image::Rgb(pixels[(y * width + x) as usize]))
            }
            None => return Err(ViewError(format!("{address} has no {layer} map yet"))),
        };
        let mut png = Vec::new();
        image::DynamicImage::ImageRgb8(img).write_to(&mut std::io::Cursor::new(&mut png), image::ImageOutputFormat::Png)?;
        Ok(png)
    }
}

/// A habitability as the page names it, `None` for none.
fn habitable_name(h: Habitability) -> Option<&'static str> {
    match h {
        Habitability::Temperate => Some("temperate"),
        Habitability::Habitable => Some("habitable"),
        Habitability::None => None,
    }
}

pub(crate) fn body_of(address: &Address, system: &SolarSystem) -> Result<geoforge_core::BodyContext> {
    Ok(geoforge_cosmos::resolve::body(address, system)?)
}

/// The real universe's galaxies: the Milky Way and those its sky lists, placed from the
/// Sun's sky (distance and direction).
fn real_galaxies(universe: &Address) -> Vec<GalaxyDot> {
    let mut out = vec![GalaxyDot {
        address: Some(universe.galaxy(GalaxyRef::MilkyWay).to_string()),
        name: "Milky Way".into(),
        kind: "SBbc".into(),
        position_mpc: [0.0; 3],
        stellar_msun: geoforge_cosmos::Galaxy::milky_way().mass_msun(),
        radius_kpc: 15.0,
        colour: type_colour(geoforge_cosmos::Galaxy::milky_way().morphology),
        home: true,
        active: None,
    }];
    let sol = geoforge_cosmos::sol::archetype(geoforge_cosmos::resolve::REAL_SEED, geoforge_core::Archetype::Earth).0;
    let observer = geoforge_cosmos::sky::Observer::of_system(&sol, 8.2);
    if let Some(site) = &sol.root.site {
        for g in geoforge_cosmos::universe::sky::seen(site, None, sol.root.seed, &geoforge_cosmos::Galaxy::milky_way(), &observer, 12.0) {
            let p = [0, 1, 2].map(|k| ((observer.position_pc[k] + g.direction[k] * g.distance_pc) / 1e6) as f32);
            let luminosity = 10f64.powf(-0.4 * (g.absolute_v - 4.83));
            out.push(GalaxyDot {
                address: None,
                name: g.name,
                kind: g.kind,
                position_mpc: p,
                stellar_msun: 2.0 * luminosity,
                radius_kpc: g.scale_deg.to_radians() * g.distance_pc * 7.0 / 1000.0,
                colour: [200, 210, 230],
                home: false,
                // The real galaxies' nuclei are not modelled.
                active: None,
            });
        }
    }
    out
}

/// A galaxy's colour by type: spirals blue-white, early types gold, dwarfs grey-blue.
fn type_colour(m: Morphology) -> [u8; 3] {
    match m {
        Morphology::Spiral { .. } => [190, 205, 255],
        Morphology::Starburst => [170, 190, 255],
        Morphology::Lenticular => [240, 220, 180],
        Morphology::Elliptical => [255, 210, 150],
        Morphology::DwarfIrregular => [170, 190, 230],
        Morphology::DwarfSpheroidal => [190, 180, 170],
    }
}

/// Starlight's colour by its age: young stars blue-white, old ones gold, Gyr.
fn age_colour(age_gyr: f64) -> [f64; 3] {
    let t = (age_gyr.max(0.01).log10() + 2.0) / 3.1;
    let t = t.clamp(0.0, 1.0);
    [0.62 + 0.38 * t, 0.72 + 0.12 * t, 1.0 - 0.45 * t]
}

fn class_name(c: BodyClass) -> &'static str {
    match c {
        BodyClass::Terrestrial => "rocky",
        BodyClass::GasGiant => "gas giant",
        BodyClass::IceGiant => "ice giant",
        BodyClass::SmallBody => "small body",
    }
}

/// A nucleus's activity, if it is active.
fn active(n: &geoforge_cosmos::nuclei::Nucleus) -> Option<String> {
    (n.activity != geoforge_cosmos::nuclei::Activity::Quiet).then(|| format!("{}{}", n.activity.label(), if n.jets { ", radio jets" } else { "" }))
}
