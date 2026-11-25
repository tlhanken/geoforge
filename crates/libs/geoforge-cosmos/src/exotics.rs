//! Where a universe's fictional materials ([`geoforge_core::exotics`]) are: deposits on its
//! worlds and in its asteroids, drawn from each material's own stream (keyed by its name),
//! so adding, removing or editing one moves no other and nothing else that is made.

use crate::asteroids;
use crate::context::{belt_member_made, body_context};
use crate::draw::Draw;
use geoforge_core::body::{BodyClass, NodeKind, ShellKind, SolarSystem, Species, StellarStage};
use geoforge_core::exotics::{Exotics, Found, Host, Material};
use geoforge_core::key::Key;
use geoforge_core::seed::{unit, Seed};
use geoforge_core::BodyContext;

/// A fictional material on one body.
#[derive(Debug, Clone, PartialEq)]
pub struct Deposit {
    pub material: String,
    pub note: String,
    pub colour: [u8; 3],
    /// How much there is.
    pub tonnes: f64,
    /// Where on the surface: latitude and longitude, degrees.
    pub sites: Vec<(f64, f64)>,
}

/// A material found in a system, and where.
#[derive(Debug, Clone, PartialEq)]
pub struct Find {
    pub material: String,
    /// "Resai V: 3 deposits" or "Main Belt: ~40 M-type bodies, e.g. #312".
    pub place: String,
    /// The body to go to: the world, or a belt's first example.
    pub body: Option<Key>,
}

/// The kinds a body is, as [`Found::kinds`] names them: a planet with seas is both rocky and
/// an ocean world.
#[must_use]
pub fn world_kinds(ctx: &BodyContext) -> Vec<&'static str> {
    match ctx.exterior.class {
        BodyClass::GasGiant | BodyClass::IceGiant => vec!["giant"],
        BodyClass::SmallBody => vec!["small"],
        BodyClass::Terrestrial => {
            // An ice shell, or a frozen-through crust mostly of ice (Callisto's).
            let icy = ctx.interior.surface_shell().is_some_and(|s| s.kind == ShellKind::IceShell || s.composition.fraction(Species::Ice) > 0.5);
            let mut kinds = vec![if icy { "icy" } else { "rocky" }];
            if ctx.interior.has_surface_liquid_water() {
                kinds.push("ocean");
            }
            kinds
        }
    }
}

/// The fictional deposits on body `key` of `system`: a node of the tree, or a belt member.
#[must_use]
pub fn on_body(exotics: &Exotics, system: &SolarSystem, key: Key) -> Vec<Deposit> {
    if exotics.is_empty() {
        return Vec::new();
    }
    let stars = stages(system);
    let member = (key.0 >> 32 != 0).then(|| belt_member_made(system, Key(key.0 >> 32), (key.0 & 0xFFFF_FFFF).saturating_sub(1))).flatten();
    let world = if member.is_some() { None } else { body_context(system, key) };
    let mut out = Vec::new();
    for m in &exotics.materials {
        for (r, f) in m.found.iter().enumerate() {
            if !f.admits_stars(&stars) || unit(stream(system, m, r).child_u64(b"body", key.0)) >= f.chance {
                continue;
            }
            let mut d = Draw::new(stream(system, m, r).child_u64(b"deposit", key.0));
            let deposit = match (f.host, &member, &world) {
                (Host::Asteroid, Some((ctx, made)), _) if f.admits(made.taxon.code()) => {
                    // 10 ppm to 1% of its mass, at one or two sites.
                    let share = 10f64.powf(d.range(-5.0, -2.0));
                    let n = 1 + usize::from(d.u() < 0.3);
                    Some((share * ctx.exterior.mass_kg / 1000.0, sites(&mut d, n)))
                }
                (Host::World, None, Some(ctx)) if world_kinds(ctx).iter().any(|k| f.admits(k)) => {
                    // A few deposits, each 10⁶–10¹⁰ t.
                    let n = (1 + d.poisson(1.5)).min(8);
                    let tonnes = (0..n).map(|_| 10f64.powf(d.range(6.0, 10.0))).sum();
                    Some((tonnes, sites(&mut d, n)))
                }
                _ => None,
            };
            if let Some((tonnes, sites)) = deposit {
                out.push(Deposit { material: m.name.clone(), note: m.note.clone(), colour: m.rgb(), tonnes, sites });
                break;
            }
        }
    }
    out
}

/// Where each material is in `system`: its worlds' deposits, and for each belt the bodies
/// expected to hold it (with up to three examples found among its first members when
/// `examples`; finding them takes a scan, so leave it off to flag many systems).
#[must_use]
pub fn in_system(exotics: &Exotics, system: &SolarSystem, examples: bool) -> Vec<Find> {
    if exotics.is_empty() {
        return Vec::new();
    }
    let mut out = Vec::new();
    let stars = stages(system);
    for n in &system.nodes {
        match &n.node {
            NodeKind::Body(_) => {
                for d in on_body(exotics, system, n.key) {
                    let what = if d.sites.len() == 1 { "1 deposit".to_string() } else { format!("{} deposits", d.sites.len()) };
                    out.push(Find { material: d.material, place: format!("{}: {what}", n.name), body: Some(n.key) });
                }
            }
            NodeKind::Belt(pop) => {
                let Some(frost_au) = crate::context::frost_line_of(system, n.key) else { continue };
                for m in &exotics.materials {
                    for (r, f) in m.found.iter().enumerate().filter(|(_, f)| f.host == Host::Asteroid && f.admits_stars(&stars)) {
                        let expected = f.chance * pop.count_estimate * asteroids::kind_share(pop, frost_au, |t| f.admits(t.code()));
                        if expected < 0.5 {
                            continue;
                        }
                        let kinds = kinds_named(f);
                        let found = if examples { members_with(system, n.key, m, r, f, pop.count_estimate) } else { Vec::new() };
                        let eg = if found.is_empty() { String::new() } else { format!(", e.g. {}", found.iter().map(|k| format!("#{}", k.0 & 0xFFFF_FFFF)).collect::<Vec<_>>().join(", ")) };
                        out.push(Find { material: m.name.clone(), place: format!("{}: ~{} {kinds}bodies{eg}", n.name, round(expected)), body: found.first().copied() });
                    }
                }
            }
            _ => {}
        }
    }
    out
}

/// The names of every material found in `system` (no scan of belts).
#[must_use]
pub fn names_in(exotics: &Exotics, system: &SolarSystem) -> Vec<String> {
    let mut names: Vec<String> = in_system(exotics, system, false).into_iter().map(|f| f.material).collect();
    names.sort();
    names.dedup();
    names
}

/// Members of belt `belt` holding material `m` by rule `r`, among its first members.
fn members_with(system: &SolarSystem, belt: Key, m: &Material, r: usize, f: &Found, count: f64) -> Vec<Key> {
    const SCAN: u64 = 200_000;
    const TRIES: usize = 200;
    let s = stream(system, m, r);
    (0..(count as u64).min(SCAN))
        .map(|i| Key((belt.0 << 32) | (i + 1)))
        .filter(|k| unit(s.child_u64(b"body", k.0)) < f.chance)
        .take(TRIES)
        .filter(|k| belt_member_made(system, belt, (k.0 & 0xFFFF_FFFF) - 1).is_some_and(|(_, a)| f.admits(a.taxon.code())))
        .take(3)
        .collect()
}

/// Material `m`'s stream for its rule `r` in `system`.
fn stream(system: &SolarSystem, m: &Material, r: usize) -> Seed {
    system.root.seed.child(b"exotic").child(m.name.trim().to_lowercase().as_bytes()).child_u64(b"rule", r as u64)
}

/// The stages of `system`'s stars.
fn stages(system: &SolarSystem) -> Vec<StellarStage> {
    system.nodes.iter().filter_map(|n| if let NodeKind::Star(s) = &n.node { Some(s.stage) } else { None }).collect()
}

/// "M-type " or "" for any.
fn kinds_named(f: &Found) -> String {
    if f.kinds.is_empty() {
        String::new()
    } else {
        format!("{} ", f.kinds.iter().map(|k| if k.eq_ignore_ascii_case("icy") { "icy".to_string() } else { format!("{}-type", k.to_uppercase()) }).collect::<Vec<_>>().join("/"))
    }
}

/// A count rounded to two figures.
fn round(x: f64) -> String {
    if x < 100.0 {
        format!("{:.0}", x.max(1.0))
    } else {
        let p = 10f64.powf(x.log10().floor() - 1.0);
        format!("{:.0}", (x / p).round() * p)
    }
}

/// `n` places spread evenly over a sphere.
fn sites(d: &mut Draw, n: usize) -> Vec<(f64, f64)> {
    (0..n).map(|_| ((2.0 * d.u() - 1.0).asin().to_degrees(), d.range(-180.0, 180.0))).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::asteroids::Taxon;

    fn list(found: Vec<Found>) -> Exotics {
        Exotics { materials: vec![Material { name: "Xylarium".into(), note: "a metastable alloy".into(), colour: None, found }] }
    }

    /// Asteroid types are named the same in the list as in the generator.
    #[test]
    fn kinds_agree() {
        assert!(Taxon::ALL.iter().map(|t| t.code()).eq(geoforge_core::exotics::ASTEROID_KINDS));
    }

    /// In the Sun's belt: a material of M-types is expected on about chance × count × the M
    /// share, examples are M-types holding it; the same list finds the same; another
    /// material leaves the first where it was; a rule for other stars finds nothing.
    #[test]
    fn found_in_belts_and_on_worlds() {
        let sol = crate::sol::at(Seed(1), 2_451_545.0);
        let ex = list(vec![Found { host: Host::Asteroid, kinds: vec!["M".into()], stars: vec![], chance: 0.01 }]);
        let finds = in_system(&ex, &sol, true);
        let main = finds.iter().find(|f| f.place.starts_with("Main Belt")).expect("in the main belt");
        let key = main.body.expect("an example");
        let made = asteroids::of(&sol, key).unwrap();
        assert_eq!(made.taxon, Taxon::M);
        let deposit = &on_body(&ex, &sol, key)[0];
        assert!(deposit.tonnes > 0.0 && !deposit.sites.is_empty() && deposit.sites.iter().all(|s| s.0.abs() <= 90.0));
        assert_eq!(finds, in_system(&ex, &sol, true), "deterministic");
        let mut two = ex.clone();
        two.materials.insert(0, Material { name: "Aurorite".into(), note: String::new(), colour: None, found: vec![Found { host: Host::World, kinds: vec![], stars: vec![], chance: 1.0 }] });
        assert!(in_system(&two, &sol, true).iter().filter(|f| f.material == "Xylarium").eq(finds.iter()), "another material moves none");
        // Every world holds Aurorite at chance 1: rocky, icy and giant alike.
        let worlds = in_system(&two, &sol, false).into_iter().filter(|f| f.material == "Aurorite").count();
        assert!(worlds >= 8, "{worlds} worlds");
        let elsewhere = list(vec![Found { host: Host::World, kinds: vec![], stars: vec![StellarStage::NeutronStar], chance: 1.0 }]);
        assert!(in_system(&elsewhere, &sol, false).is_empty());
        assert!(in_system(&Exotics::default(), &sol, true).is_empty());
    }

    /// Earth is both rocky and an ocean world, a giant a giant.
    #[test]
    fn kinds_of_worlds() {
        let sol = crate::sol::at(Seed(1), 2_451_545.0);
        let kinds = |name: &str| world_kinds(&crate::sol::known_context(&sol, sol.by_name(name).unwrap().key).unwrap());
        assert_eq!(kinds("Earth"), ["rocky", "ocean"]);
        assert_eq!(kinds("Jupiter"), ["giant"]);
        assert_eq!(kinds("Europa"), ["icy"]);
    }
}
