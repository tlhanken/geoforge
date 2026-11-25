//! Any item of a universe from its address: the galaxy, its stars, the system at a star, a
//! body of that system. Each function generates its item from the address alone (and, for
//! stars, the edits a world keeps over them); none knows about directories, worlds or what
//! is shown. A universe's default item ([`default_address`]) is just an address chosen by a
//! fixed rule.
//!
//! The real universe holds the Milky Way and the Sun's system as they are; the Milky Way's
//! other stars are generated, from [`REAL_SEED`].

use crate::galactic::{system_at, Barren};
use crate::galaxy::Galaxy;
use crate::starfield::Starfield;
use crate::system::Generated;
use crate::universe::{Universe, UniverseGalaxy};
use crate::worlds::{habitability, Habitability};
use crate::{any_body, sol};
use geoforge_core::address::{Address, GalaxyRef, Level, StarRef};
use geoforge_core::body::{BodyContext, SolarSystem, StellarStage};
use geoforge_core::galactic::{StarEdits, StarId, UniverseSite};
use geoforge_core::seed::Seed;
use geoforge_core::{Archetype, Key};
use std::fmt;

/// The seed of what the real universe generates: the Milky Way's stars and the Solar
/// System's details.
pub const REAL_SEED: Seed = Seed(0x5013_A12E);

/// Why an address names nothing.
#[derive(Debug)]
pub enum Unresolved {
    /// The address stops above the level asked for.
    Shallow { needs: Level },
    /// No such galaxy, star or body.
    Missing(String),
    /// The star's system has no solid body (its bodies cannot be generated).
    Barren(Barren),
}

impl fmt::Display for Unresolved {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Shallow { needs } => write!(f, "the address stops above a {needs:?}"),
            Self::Missing(what) => write!(f, "no {what}"),
            Self::Barren(b) => write!(f, "{b}"),
        }
    }
}

impl std::error::Error for Unresolved {}

/// A galaxy, ready to draw stars from.
#[derive(Debug, Clone)]
pub struct GalaxyAt {
    /// Its address (to the galaxy's level).
    pub address: Address,
    /// The generated galaxy; `None` for the Milky Way.
    pub generated: Option<UniverseGalaxy>,
    /// Its model: where its stars are, how old, how metal-rich.
    pub model: Galaxy,
    /// The seed of its stars.
    pub stars_seed: Seed,
}

impl GalaxyAt {
    /// Its stars, with `edits` over them.
    #[must_use]
    pub fn starfield(&self, edits: StarEdits) -> Starfield {
        Starfield::new(&self.model, self.stars_seed, edits)
    }
}

/// The galaxy at `address`.
///
/// # Errors
/// If the address names no galaxy.
pub fn galaxy(address: &Address) -> Result<GalaxyAt, Unresolved> {
    let address = address.up_to(Level::Galaxy);
    match (address.universe, address.galaxy) {
        (_, None) => Err(Unresolved::Shallow { needs: Level::Galaxy }),
        (UniverseSite::Real, Some(GalaxyRef::MilkyWay)) => {
            Ok(GalaxyAt { address, generated: None, model: Galaxy::milky_way(), stars_seed: sol_system().root.site.and_then(|s| s.galaxy_seed).unwrap_or(REAL_SEED) })
        }
        (UniverseSite::Seeded(seed), Some(GalaxyRef::Generated(id))) => {
            let g = Universe::new(seed).galaxy(&id).ok_or_else(|| Unresolved::Missing(format!("galaxy {id} in universe {}", seed.0)))?;
            Ok(GalaxyAt { address, model: g.model(), stars_seed: g.stars_seed(), generated: Some(g) })
        }
        _ => Err(Unresolved::Missing(format!("galaxy at {address}"))),
    }
}

/// The system at `address`, its stars `field` (the address's galaxy's, as
/// [`GalaxyAt::starfield`] makes them): the system born with its star, placed at its
/// address.
///
/// # Errors
/// If the address stops above a system, the star is not in the field, or its system has no
/// solid body.
pub fn system(address: &Address, field: &Starfield) -> Result<Generated, Unresolved> {
    let address = address.up_to(Level::System);
    let star = address.star.ok_or(Unresolved::Shallow { needs: Level::System })?;
    let mut g = match star {
        StarRef::Sun => {
            let system = sol_system();
            let focus = system.by_name("Earth").expect("the Sun's system has the Earth").key;
            Generated { system, focus }
        }
        StarRef::Star(id) => system_at(field, &id).ok_or_else(|| Unresolved::Missing(format!("star {id} in {}", field.galaxy().name)))?.map_err(Unresolved::Barren)?,
    };
    place(&mut g.system, address, field);
    Ok(g)
}

/// The system at `address` as [`system`] makes it, a barren one too: its tree, and its focus
/// if it has a solid body (a lone black hole's system is its hole alone).
///
/// # Errors
/// If the address stops above a system or the star is not in the field.
pub fn system_or_barren(address: &Address, field: &Starfield) -> Result<(SolarSystem, Option<Key>), Unresolved> {
    match system(address, field) {
        Ok(g) => Ok((g.system, Some(g.focus))),
        Err(Unresolved::Barren(b)) => {
            let mut sys = *b.system;
            place(&mut sys, address.up_to(Level::System), field);
            Ok((sys, None))
        }
        Err(e) => Err(e),
    }
}

/// Mark `system` as the one at `address` of `field`'s galaxy.
fn place(system: &mut SolarSystem, address: Address, field: &Starfield) {
    system.root.placement = Some(address);
    if let Some(site) = system.root.site.as_mut() {
        site.universe = Some(address.universe);
        site.galaxy.clone_from(&field.galaxy().name);
    }
}

/// The body at `address` of its `system` (the address's, as [`system`] makes it): for the
/// real Solar System's bodies with measured contexts, those.
///
/// # Errors
/// If the address stops above a body or the system has no such body.
pub fn body(address: &Address, system: &SolarSystem) -> Result<BodyContext, Unresolved> {
    let key = address.body.ok_or(Unresolved::Shallow { needs: Level::Body })?;
    let missing = || Unresolved::Missing(format!("body {key} in {}", system.root.name));
    if address.star == Some(StarRef::Sun) {
        if let Some(a) = system.node(key).and_then(|n| Archetype::ALL.into_iter().find(|a| a.real_body() == n.name)) {
            return Ok(sol::archetype(REAL_SEED, a).1);
        }
    }
    // A node of the tree, or a belt's member.
    any_body(system, key).ok_or_else(missing)
}

/// The real Solar System as `real/milky-way/sun` holds it: at the 2026 March equinox.
fn sol_system() -> SolarSystem {
    sol::archetype(REAL_SEED, Archetype::Earth).0
}

/// Stars searched around a galaxy's reference point for a Sun-like one with a wet world,
/// nearest first, before settling for any wet world, then any.
const SEARCH_RADII_PC: [f64; 4] = [30.0, 60.0, 120.0, 240.0];

/// A universe's default body ([`complete`] from the universe down).
#[must_use]
pub fn default_address(seed: Seed) -> Address {
    complete(&Address::universe(seed)).expect("every universe has a default body")
}

/// `address` taken down to a body by the default rule at each level it stops above:
/// - a universe: its galaxy most like the Milky Way nearest its origin
///   ([`Universe::home`]); the real universe: the Sun's system;
/// - a galaxy: around a point on its reference circle (at an azimuth from its stars' seed),
///   the nearest single main-sequence Sun-like star (5,200–6,200 K) whose system's most
///   Earth-like body is temperate ([`Habitability::Temperate`]); failing that within 240 pc,
///   the nearest such star whose focus has surface water, then the nearest star of any kind
///   whose focus has, then any star with a solid body;
/// - a system: its most Earth-like body (its focus).
///
/// # Errors
/// If the address names nothing, or a galaxy has no star with a solid body near its
/// reference circle.
pub fn complete(address: &Address) -> Result<Address, Unresolved> {
    let mut a = *address;
    if a.level() == Level::Universe {
        a = match a.universe {
            UniverseSite::Real => Address::sol(),
            UniverseSite::Seeded(seed) => a.galaxy(GalaxyRef::Generated(Universe::new(seed).home(&Galaxy::milky_way()).id)),
        };
    }
    let g = galaxy(&a)?;
    let field = g.starfield(StarEdits::default());
    if a.level() == Level::Galaxy {
        a = default_star(&g, &field).ok_or_else(|| Unresolved::Missing(format!("star with a solid body near {}'s reference circle", g.model.name)))?;
    }
    if a.level() == Level::System {
        a = a.body(system(&a, &field)?.focus);
    }
    Ok(a)
}

/// The default star's body of galaxy `g` (see [`complete`]).
fn default_star(g: &GalaxyAt, field: &Starfield) -> Option<Address> {
    let phi = std::f64::consts::TAU * (g.stars_seed.child(b"default star").0 as f64 / u64::MAX as f64);
    let r = g.model.reference_kpc * 1000.0;
    let centre = [r * phi.cos(), r * phi.sin(), 0.0];
    let sunlike = |s: &crate::sectors::GalacticStar| s.kind.companion_mass.is_none() && s.star.stage == StellarStage::MainSequence && (5200.0..6200.0).contains(&s.star.temperature_k);
    // The star's focus, and how habitable it is (`None` for a dry one).
    let focus = |id: StarId| -> Option<(Address, Option<Habitability>)> {
        let a = g.address.star(StarRef::Star(id));
        let sys = system(&a, field).ok()?;
        let at = a.body(sys.focus);
        let b = body(&at, &sys.system).ok()?;
        Some((at, b.interior.has_surface_liquid_water().then(|| habitability(&b))))
    };
    let mut seen = std::collections::HashSet::new();
    let (mut wet_sunlike, mut wet_any) = (None, None);
    for radius in SEARCH_RADII_PC {
        let mut stars = field.near(centre, radius, 16.0);
        stars.retain(|s| seen.insert(s.id));
        stars.sort_by(|a, b| a.distance_pc(centre).total_cmp(&b.distance_pc(centre)));
        for s in stars.iter().filter(|s| sunlike(s)) {
            match focus(s.id) {
                Some((a, Some(Habitability::Temperate))) => return Some(a),
                Some((a, Some(_))) if wet_sunlike.is_none() => wet_sunlike = Some(a),
                _ => {}
            }
        }
        if wet_any.is_none() {
            wet_any = stars.iter().take(400).find_map(|s| focus(s.id).and_then(|(a, h)| h.map(|_| a)));
        }
    }
    wet_sunlike.or(wet_any).or_else(|| field.near(centre, SEARCH_RADII_PC[0], 16.0).iter().find_map(|s| focus(s.id).map(|(a, _)| a)))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Universes' default bodies are temperate worlds of Sun-like stars.
    #[test]
    fn default_bodies_are_temperate() {
        for seed in [1, 2, 3, 42, 257, 525, 820, 12345] {
            let a = default_address(Seed(seed));
            let field = galaxy(&a).unwrap().starfield(StarEdits::default());
            let sys = system(&a, &field).unwrap();
            let b = body(&a, &sys.system).unwrap();
            assert_eq!(habitability(&b), Habitability::Temperate, "seed {seed}: {a}");
        }
    }

    /// Every level resolves from the address alone, the same every time; the default body
    /// is a wet world of a Sun-like star in a galaxy like the Milky Way.
    #[test]
    fn a_universe_resolves_down_to_a_body() {
        let a = default_address(Seed(7));
        assert_eq!(a.level(), Level::Body);
        assert_eq!(a, default_address(Seed(7)), "deterministic");
        let g = galaxy(&a).unwrap();
        assert_eq!(g.model.morphology, Galaxy::milky_way().morphology);
        let field = g.starfield(StarEdits::default());
        let sys = system(&a, &field).unwrap();
        assert_eq!(sys.system.root.placement, Some(a.up_to(Level::System)));
        let site = sys.system.root.site.as_ref().unwrap();
        assert_eq!((site.universe, site.galaxy.as_str()), (Some(UniverseSite::Seeded(Seed(7))), g.model.name.as_str()));
        let b = body(&a, &sys.system).unwrap();
        assert!(b.interior.has_surface_liquid_water());
        let star = &sys.system.nodes.iter().find_map(|n| if let geoforge_core::body::NodeKind::Star(s) = &n.node { Some(s.clone()) } else { None }).unwrap();
        assert!((5200.0..6200.0).contains(&star.temperature_k), "{} K", star.temperature_k);
        // Another star of the same galaxy is another system of the same universe.
        let other = field
            .near(
                field
                    .star(&match a.star {
                        Some(StarRef::Star(id)) => id,
                        _ => unreachable!(),
                    })
                    .unwrap()
                    .position_pc,
                20.0,
                12.0,
            )
            .into_iter()
            .find_map(|s| system(&a.star(StarRef::Star(s.id)), &field).ok())
            .unwrap();
        assert_eq!(other.system.root.site.as_ref().unwrap().galaxy, g.model.name);
        assert_ne!(other.system.root.seed, sys.system.root.seed);
    }

    #[test]
    fn the_real_universe_has_the_sun() {
        let a = Address::sol();
        let g = galaxy(&a).unwrap();
        assert_eq!(g.model.name, "Milky Way");
        let sys = system(&a, &g.starfield(StarEdits::default())).unwrap();
        let earth = sys.system.by_name("Earth").unwrap().key;
        assert_eq!(body(&a.body(earth), &sys.system).unwrap().name, "Earth");
        assert!(matches!(galaxy(&Address::universe(Seed(1))), Err(Unresolved::Shallow { needs: Level::Galaxy })));
        assert_eq!(complete(&Address::real()).unwrap(), a.body(earth));
    }
}
