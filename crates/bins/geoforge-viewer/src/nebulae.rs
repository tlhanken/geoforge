//! Nebulae as the viewer shows them ([`geoforge_cosmos::nebulae`]): round a galaxy's home star
//! for its view, round a system for the stars around it. They are not places of their own:
//! each is named by what makes it and leads to the star that lights or left it. A galaxy's
//! globular clusters are drawn the same way ([`globular_dots`]).

use crate::{star_id, Result, ViewError, Viewer};
use geoforge_core::address::{Address, Level, StarRef};
use geoforge_cosmos::nebulae::{self, NebulaKind};
use serde::Serialize;

/// How far round a galaxy's home star its view shows nebulae, pc.
pub const GALAXY_REACH_PC: f64 = 3000.0;
/// How far round a system the stars around it list nebulae, pc.
pub const STARS_REACH_PC: f64 = 500.0;

/// A nebula or globular cluster to draw and tell about.
#[derive(Debug, Clone, Serialize)]
pub struct DeepSkyDot {
    pub id: String,
    /// emission, reflection, dark, planetary, supernova_remnant or globular.
    pub kind: &'static str,
    pub name: String,
    /// Where it is in the galaxy, kpc, and from the centre asked about, pc.
    pub position_kpc: [f64; 3],
    pub offset_pc: [f64; 3],
    pub distance_pc: f64,
    pub radius_pc: f64,
    pub facts: Vec<[String; 2]>,
    /// The system of the star that lights it or left it.
    pub star: Option<String>,
    /// The fictional materials it holds.
    pub materials: Vec<String>,
}

impl Viewer {
    /// The nebulae round `address`: within [`GALAXY_REACH_PC`] of a galaxy's home star for a
    /// galaxy, within [`STARS_REACH_PC`] of a system (or body) for anything deeper; nearest
    /// first.
    ///
    /// # Errors
    /// If the galaxy or system cannot be generated, or has no place.
    pub fn nebulae(&self, address: &Address) -> Result<Vec<DeepSkyDot>> {
        let g = self.galaxy(address)?;
        let galaxy = address.up_to(Level::Galaxy);
        let (centre, reach) = if address.level() == Level::Galaxy {
            let Some(home) = self.home_star(&galaxy) else { return Ok(Vec::new()) };
            (self.place_of(&home)?, GALAXY_REACH_PC)
        } else {
            (self.place_of(&address.up_to(Level::System))?, STARS_REACH_PC)
        };
        let found = if reach == GALAXY_REACH_PC { g.nebulae.get_or_init(|| nebulae::near(g.field.sectors(), centre, reach)).clone() } else { nebulae::near(g.field.sectors(), centre, reach) };
        let exotics = self.exotics().unwrap_or_default();
        Ok(found
            .iter()
            .map(|n| {
                let star = n.star.map(|id| galaxy.star(StarRef::Star(id)));
                let mut facts = vec![["Kind".to_string(), n.kind.label().to_string()], ["Across".into(), format!("{:.1} pc ({:.1} ly)", 2.0 * n.radius_pc, 2.0 * n.radius_pc * 3.2616)]];
                if let Some(age) = n.age_myr {
                    let what = if matches!(n.kind, NebulaKind::Planetary | NebulaKind::SupernovaRemnant) { "Since the star died" } else { "Its stars' age" };
                    facts.push([what.into(), if age < 1.0 { format!("{:.0} thousand years", age * 1000.0) } else { format!("{age:.1} million years") }]);
                }
                if let Some(m) = n.mass_msun {
                    facts.push(["Gas".into(), format!("{m:.1e} M☉")]);
                }
                if let Some(q) = n.ionising_per_s {
                    facts.push(["Ionising light".into(), format!("10^{:.1} photons/s", q.log10())]);
                }
                if let (Some(id), Some(stage)) = (n.star, n.star_stage) {
                    facts.push(["Its star".into(), format!("{} ({})", g.field.name_of(&id), stage_name(stage))]);
                }
                let materials: Vec<String> = geoforge_cosmos::exotics::in_nebula(&exotics, n).into_iter().map(|m| m.name.clone()).collect();
                for m in &materials {
                    facts.push([format!("{m} (fictional)"), "found here".into()]);
                }
                DeepSkyDot {
                    id: n.id.to_string(),
                    kind: n.kind.code(),
                    name: n.name.clone(),
                    position_kpc: n.centre_pc.map(|v| v / 1000.0),
                    offset_pc: [0, 1, 2].map(|k| n.centre_pc[k] - centre[k]),
                    distance_pc: n.distance_pc(centre),
                    radius_pc: n.radius_pc,
                    facts,
                    star: star.map(|a| a.to_string()),
                    materials,
                }
            })
            .collect())
    }

    /// Where system `address`'s star is in its galaxy, pc.
    fn place_of(&self, address: &Address) -> Result<[f64; 3]> {
        let g = self.galaxy(address)?;
        match star_id(address) {
            Some(id) => Ok(g.field.star(&id).ok_or_else(|| ViewError(format!("no star {id}")))?.position_pc),
            None => self.system(address)?.system.root.site.and_then(|s| s.position_pc).ok_or_else(|| ViewError("the system has no place".into())),
        }
    }
}

fn stage_name(stage: geoforge_core::body::StellarStage) -> &'static str {
    use geoforge_core::body::StellarStage as S;
    match stage {
        S::MainSequence => "a dwarf",
        S::Subgiant => "a subgiant",
        S::Giant => "a giant",
        S::Supergiant => "a supergiant",
        S::WhiteDwarf => "a white dwarf",
        S::NeutronStar => "a neutron star",
        S::BlackHole => "a black hole",
        S::BrownDwarf => "a brown dwarf",
    }
}

/// The nebulae of `found` whose names (or fictional materials) hold `q`, as search hits
/// leading to their stars (a dark cloud, which has none, to `galaxy`).
pub(crate) fn named(found: &[DeepSkyDot], q: &str, galaxy: &Address) -> Vec<crate::search::Hit> {
    found
        .iter()
        .filter(|n| n.name.to_lowercase().contains(q) || n.materials.iter().any(|m| m.to_lowercase().contains(q)))
        .map(|n| {
            let note = if n.materials.is_empty() {
                format!("{}, {:.0} ly off", kind_label(n.kind), n.distance_pc * 3.2616)
            } else {
                format!("{}: holds {} (fictional)", kind_label(n.kind), n.materials.join(", "))
            };
            crate::search::Hit { name: n.name.clone(), kind: "nebula".into(), address: n.star.clone().unwrap_or_else(|| galaxy.to_string()), note }
        })
        .collect()
}

fn kind_label(code: &str) -> &'static str {
    NebulaKind::ALL.into_iter().find(|k| k.code() == code).map_or("nebula", NebulaKind::label)
}

/// `galaxy`'s globular clusters as dots, their distances from `centre` (pc).
#[must_use]
pub fn globular_dots(galaxy: &geoforge_cosmos::resolve::GalaxyAt, centre: [f64; 3]) -> Vec<DeepSkyDot> {
    galaxy
        .globulars()
        .into_iter()
        .map(|g| DeepSkyDot {
            id: format!("gc{}", g.index),
            kind: "globular",
            name: g.name.clone(),
            position_kpc: g.position_pc.map(|v| v / 1000.0),
            offset_pc: [0, 1, 2].map(|k| g.position_pc[k] - centre[k]),
            distance_pc: g.distance_pc(centre),
            radius_pc: 3.0 * g.half_light_pc,
            facts: vec![
                ["Kind".into(), "globular cluster".into()],
                ["Stars".into(), format!("{:.1e} M☉", g.mass_msun)],
                ["Age".into(), format!("{:.1} billion years", g.age_gyr)],
                ["Metals [Fe/H]".into(), format!("{:+.2}", g.metallicity_dex)],
                ["Half its light within".into(), format!("{:.1} pc", g.half_light_pc)],
                ["From the centre".into(), format!("{:.1} kpc", g.distance_pc([0.0; 3]) / 1000.0)],
            ],
            star: None,
            materials: Vec::new(),
        })
        .collect()
}

/// What sits at `galaxy`'s centre, in a line: "black hole of 4.3×10⁶ M☉, quiet".
#[must_use]
pub fn centre_line(n: &geoforge_cosmos::nuclei::Nucleus) -> String {
    match n.black_hole_msun {
        None => "no central black hole".into(),
        Some(m) => format!("black hole of {m:.1e} M☉, {}{}", n.activity.label(), if n.jets { ", with radio jets" } else { "" }),
    }
}
