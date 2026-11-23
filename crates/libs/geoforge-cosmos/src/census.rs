//! A census of systems drawn from a galaxy: what their stars are, what orbits them, how often
//! the focus is a wet world, and example seeds of each kind, to find one worth making.

use crate::context::any_body;
use crate::galactic::galactic_system;
use crate::galaxy::Galaxy;
use crate::population::draw;
use geoforge_core::body::{BodyClass, NodeKind, SolarSystem, Star, StellarStage};
use geoforge_core::seed::Seed;
use std::ops::Range;

/// Kinds of primary star, in the order a census lists them.
pub const KINDS: [&str; 14] = ["O", "B", "A", "F", "G", "K", "M", "brown dwarf", "subgiant", "giant", "supergiant", "white dwarf", "neutron star", "black hole"];

/// The census kind of a primary: its spectral letter on the main sequence, else its stage.
#[must_use]
pub fn kind_of(star: &Star) -> &'static str {
    match star.stage {
        StellarStage::MainSequence => KINDS.iter().take(7).find(|k| star.spectral_class.starts_with(**k)).copied().unwrap_or("M"),
        StellarStage::Subgiant => "subgiant",
        StellarStage::Giant => "giant",
        StellarStage::Supergiant => "supergiant",
        StellarStage::WhiteDwarf => "white dwarf",
        StellarStage::NeutronStar => "neutron star",
        StellarStage::BlackHole => "black hole",
        StellarStage::BrownDwarf => "brown dwarf",
    }
}

/// Counts over the systems of one kind of primary.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Kind {
    pub name: &'static str,
    pub systems: u64,
    /// Systems of more than one star now (a supernova unbinds most companions).
    pub multiples: u64,
    /// Planets (not moons), by class.
    pub rocky: u64,
    pub giants: u64,
    pub ice_giants: u64,
    /// Systems with no solid body at all.
    pub barren: u64,
    /// Systems whose focus has liquid water on its surface.
    pub wet: u64,
    /// Systems whose focus is a planet tidally locked to its star.
    pub locked: u64,
    /// The first few seeds of this kind that make a world (not barren), and of this kind with a
    /// wet focus.
    pub examples: Vec<u64>,
    pub wet_examples: Vec<u64>,
}

/// A census of the systems of `seeds` in `galaxy`.
#[derive(Debug, Clone, PartialEq)]
pub struct Census {
    pub galaxy: String,
    pub seeds: Range<u64>,
    /// One entry per kind in [`KINDS`] order, including kinds with no systems.
    pub kinds: Vec<Kind>,
}

impl Census {
    #[must_use]
    pub fn systems(&self) -> u64 {
        self.kinds.iter().map(|k| k.systems).sum()
    }
    /// Sum of `f` over kinds.
    pub fn total(&self, f: impl Fn(&Kind) -> u64) -> u64 {
        self.kinds.iter().map(f).sum()
    }
    #[must_use]
    pub fn kind(&self, name: &str) -> Option<&Kind> {
        self.kinds.iter().find(|k| k.name == name)
    }
}

/// How many example seeds a census keeps per kind.
const EXAMPLES: usize = 4;

/// Draw the systems of `seeds` from `galaxy` and count them.
#[must_use]
pub fn census(galaxy: &Galaxy, seeds: Range<u64>) -> Census {
    let mut kinds: Vec<Kind> = KINDS.iter().map(|&name| Kind { name, ..Kind::default() }).collect();
    for s in seeds.clone() {
        let ph = draw(galaxy, Seed(s));
        let k = &mut kinds[KINDS.iter().position(|n| *n == kind_of(&ph.primary())).unwrap_or(6)];
        k.systems += 1;
        let stars = |sys: &SolarSystem| sys.nodes.iter().filter(|n| matches!(n.node, NodeKind::Star(_))).count();
        let g = match galactic_system(Seed(s), galaxy) {
            Ok(g) => g,
            Err(b) => {
                k.barren += 1;
                k.multiples += u64::from(stars(&b.system) > 1);
                continue;
            }
        };
        k.multiples += u64::from(stars(&g.system) > 1);
        if k.examples.len() < EXAMPLES {
            k.examples.push(s);
        }
        for n in &g.system.nodes {
            let top = n.parent.and_then(|p| g.system.node(p)).is_some_and(|p| !matches!(p.node, NodeKind::Body(_)));
            if let (true, NodeKind::Body(b)) = (top, &n.node) {
                match b.class {
                    BodyClass::Terrestrial => k.rocky += 1,
                    BodyClass::GasGiant => k.giants += 1,
                    BodyClass::IceGiant => k.ice_giants += 1,
                    BodyClass::SmallBody => {}
                }
            }
        }
        let Some(focus) = any_body(&g.system, g.focus) else { continue };
        if focus.interior.has_surface_liquid_water() {
            k.wet += 1;
            if k.wet_examples.len() < EXAMPLES {
                k.wet_examples.push(s);
            }
        }
        let is_planet = g.system.node(g.focus).and_then(|n| n.parent).and_then(|p| g.system.node(p)).is_some_and(|p| !matches!(p.node, NodeKind::Body(_)));
        let spin = focus.exterior.rotation_h / (focus.irradiance.year_days * 24.0);
        if is_planet && ((spin - 1.0).abs() < 0.01 || (spin - 2.0 / 3.0).abs() < 0.01) {
            k.locked += 1;
        }
    }
    Census { galaxy: galaxy.name.clone(), seeds, kinds }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn milky_way_systems() {
        let c = census(&Galaxy::milky_way(), 0..3000);
        let n = c.systems() as f64;
        assert_eq!(c.systems(), 3000);
        let barren = c.total(|k| k.barren) as f64 / n;
        assert!(barren < 0.25, "barren {barren}");
        // Red dwarf worlds in the habitable zone are tidally locked.
        let m = c.kind("M").unwrap();
        assert!(m.wet > 0 && m.locked as f64 > 0.8 * m.wet as f64, "M: wet {} locked {}", m.wet, m.locked);
        // White dwarfs lost their inner planets; some kept outer ones.
        let wd = c.kind("white dwarf").unwrap();
        assert!(wd.systems > 100 && wd.giants + wd.ice_giants > 0);
        assert_eq!(c.kind("black hole").unwrap().systems, c.kind("black hole").unwrap().barren);
    }
}
