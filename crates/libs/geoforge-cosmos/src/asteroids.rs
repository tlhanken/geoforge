//! What a belt's small bodies are made of: their spectral types, the elements in them, how
//! dense and porous they are, and the families collisions leave.
//!
//! - **Types** follow the share of mass in each class by distance from the star in the
//!   Solar System's belt (DeMeo & Carry 2014), measured against the frost line (2.7 AU for the
//!   Sun) so a dimmer star's belt is laid out the same: stony S-types and basaltic V-types
//!   inside, carbonaceous C-types in the middle, metal M-types a few per cent throughout, dark
//!   P- and D-types outward, D-types most of a giant's Trojans. A belt mostly of ice (a Kuiper
//!   belt, a giant's irregular moons beyond the frost line) holds icy bodies.
//! - **Elements** are each type's meteorite: CI chondrites for C (Lodders 2003), L chondrites
//!   for S (Wasson & Kallemeyn 1988), iron meteorites for M (the IAB–IIIAB mean), eucrites for
//!   V, Tagish Lake for P and D (Brown et al. 2000), and for icy bodies 67P's mix of ices,
//!   organics and dust; each body's spread ±15% (rescaled to one). Metal is the free iron–nickel: 7% of an
//!   L chondrite, nearly all of an iron.
//! - **Density** is the meteorite's grain density less the voids of a rubble pile: bodies
//!   from 0.2 to 100 km across are piles of fragments 20–50% empty (Carry 2012), comets
//!   40–75%; bodies of a few hundred km keep some voids, the largest and the smallest
//!   (monoliths) few.
//! - **Families**: about a third of a belt's members belong to one of 5–30 families
//!   (Nesvorný et al. 2015), fragments of a broken parent sharing its type and, within a
//!   hundredth of an AU and a few hundredths of eccentricity, its orbit.

use crate::draw::Draw;
use geoforge_core::body::{Elements, Population, SolarSystem, Species};
use geoforge_core::key::Key;
use geoforge_core::seed::Seed;

/// A small body's spectral type.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Taxon {
    /// Stony: olivine and pyroxene with some metal (ordinary chondrites).
    S,
    /// Carbonaceous: clays, organics and water bound in minerals.
    C,
    /// Metallic: iron–nickel, a broken core.
    M,
    /// Basaltic: crust from a differentiated parent (Vesta's).
    V,
    /// Primitive, dark and reddish: carbon-rich, some ice.
    P,
    /// Primitive, darker and redder still: organics and ice.
    D,
    /// Water ice, dust and frozen volatiles: a comet's or a Kuiper belt object's.
    Icy,
}

impl Taxon {
    pub const ALL: [Self; 7] = [Self::S, Self::C, Self::M, Self::V, Self::P, Self::D, Self::Icy];

    /// Its letter ("icy" for icy bodies).
    #[must_use]
    pub fn code(self) -> &'static str {
        match self {
            Self::S => "S",
            Self::C => "C",
            Self::M => "M",
            Self::V => "V",
            Self::P => "P",
            Self::D => "D",
            Self::Icy => "icy",
        }
    }

    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Self::S => "S-type (stony)",
            Self::C => "C-type (carbonaceous)",
            Self::M => "M-type (metallic)",
            Self::V => "V-type (basaltic)",
            Self::P => "P-type (primitive)",
            Self::D => "D-type (primitive, organic-rich)",
            Self::Icy => "icy",
        }
    }

    /// Grain density (no voids), g/cm³: meteorites' (Consolmagno et al. 2008); for icy
    /// bodies a rock–ice mix's (67P's nucleus, Pluto's and Eris's bulk).
    #[must_use]
    pub fn grain_density(self) -> f64 {
        match self {
            Self::S => 3.6,
            Self::C => 2.5,
            Self::M => 7.8,
            Self::V => 3.2,
            Self::P => 2.0,
            Self::D => 1.8,
            Self::Icy => 1.9,
        }
    }

    /// Mass fractions of its meteorite analog, by element.
    #[must_use]
    pub fn analog(self) -> &'static [(Element, f64)] {
        use Element as E;
        match self {
            // CI chondrite (Lodders 2003), its water bound as H and O.
            Self::C => &[
                (E::O, 0.460),
                (E::Fe, 0.185),
                (E::Si, 0.107),
                (E::Mg, 0.096),
                (E::S, 0.054),
                (E::C, 0.035),
                (E::H, 0.0197),
                (E::Ni, 0.011),
                (E::Ca, 0.0092),
                (E::Al, 0.0085),
                (E::Na, 0.005),
                (E::Cr, 0.0026),
                (E::Mn, 0.0019),
                (E::P, 0.0009),
                (E::Co, 0.0005),
                (E::N, 0.0003),
                (E::Pgm, 3.4e-6),
                (E::Au, 1.5e-7),
            ],
            // L chondrite.
            Self::S => &[
                (E::O, 0.377),
                (E::Fe, 0.218),
                (E::Si, 0.185),
                (E::Mg, 0.149),
                (E::S, 0.022),
                (E::Ca, 0.013),
                (E::Ni, 0.012),
                (E::Al, 0.012),
                (E::Na, 0.007),
                (E::Cr, 0.0039),
                (E::Mn, 0.0026),
                (E::Co, 0.0006),
                (E::C, 0.0015),
                (E::Pgm, 2.9e-6),
                (E::Au, 1.6e-7),
            ],
            // Iron meteorites.
            Self::M => &[(E::Fe, 0.910), (E::Ni, 0.080), (E::Co, 0.005), (E::P, 0.002), (E::S, 0.002), (E::C, 0.001), (E::Pgm, 3.5e-5), (E::Au, 1.0e-6)],
            // Eucrite: basalt.
            Self::V => &[(E::O, 0.430), (E::Si, 0.230), (E::Fe, 0.145), (E::Ca, 0.074), (E::Al, 0.067), (E::Mg, 0.042), (E::Ti, 0.004), (E::Mn, 0.004), (E::Na, 0.003), (E::Cr, 0.002), (E::Pgm, 2e-8)],
            // Tagish Lake.
            Self::P | Self::D => &[
                (E::O, 0.440),
                (E::Fe, 0.180),
                (E::Si, 0.110),
                (E::Mg, 0.105),
                (E::C, 0.054),
                (E::S, 0.045),
                (E::H, 0.020),
                (E::Ni, 0.010),
                (E::Ca, 0.010),
                (E::Al, 0.009),
                (E::Na, 0.005),
                (E::N, 0.002),
                (E::Pgm, 3e-6),
            ],
            // A comet's or a Kuiper belt object's: 40% water ice, 10% other ices (CO, CO₂,
            // CH₄, NH₃), 15% refractory organics, 35% CI-like dust.
            Self::Icy => &[
                (E::O, 0.596),
                (E::C, 0.130),
                (E::H, 0.065),
                (E::Fe, 0.065),
                (E::Si, 0.037),
                (E::Mg, 0.034),
                (E::N, 0.020),
                (E::S, 0.019),
                (E::Ni, 0.004),
                (E::Ca, 0.003),
                (E::Al, 0.003),
                (E::Na, 0.002),
                (E::Pgm, 1.2e-6),
            ],
        }
    }

    /// Share of its mass that is free iron–nickel metal.
    #[must_use]
    pub fn metal_share(self) -> f64 {
        match self {
            Self::S => 0.07,
            Self::M => 0.98,
            _ => 0.0,
        }
    }
}

/// The elements tallied in a small body.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Element {
    H,
    C,
    N,
    O,
    Na,
    Mg,
    Al,
    Si,
    P,
    S,
    Ca,
    Ti,
    Cr,
    Mn,
    Fe,
    Co,
    Ni,
    /// The platinum-group metals together: Ru, Rh, Pd, Os, Ir, Pt.
    Pgm,
    Au,
}

impl Element {
    #[must_use]
    pub fn symbol(self) -> &'static str {
        match self {
            Self::H => "H",
            Self::C => "C",
            Self::N => "N",
            Self::O => "O",
            Self::Na => "Na",
            Self::Mg => "Mg",
            Self::Al => "Al",
            Self::Si => "Si",
            Self::P => "P",
            Self::S => "S",
            Self::Ca => "Ca",
            Self::Ti => "Ti",
            Self::Cr => "Cr",
            Self::Mn => "Mn",
            Self::Fe => "Fe",
            Self::Co => "Co",
            Self::Ni => "Ni",
            Self::Pgm => "Pt-group",
            Self::Au => "Au",
        }
    }
}

/// A collisional family.
#[derive(Debug, Clone, PartialEq)]
pub struct Family {
    pub index: u32,
    pub name: String,
    pub taxon: Taxon,
    /// Its members' shared orbit.
    pub a_au: f64,
    pub e: f64,
    pub i_deg: f64,
}

/// A small body of a belt as it is made.
#[derive(Debug, Clone, PartialEq)]
pub struct Asteroid {
    pub taxon: Taxon,
    /// Mass fractions by element, summing to one.
    pub makeup: Vec<(Element, f64)>,
    /// Share of its volume that is empty.
    pub porosity: f64,
    /// Bulk density, g/cm³.
    pub density_g_cc: f64,
    pub family: Option<Family>,
}

/// Share of a belt's members in families.
const FAMILY_SHARE: f64 = 0.35;

impl Asteroid {
    /// Share of its mass that is `element`.
    #[must_use]
    pub fn share(&self, element: Element) -> f64 {
        self.makeup.iter().find(|x| x.0 == element).map_or(0.0, |x| x.1)
    }

    /// Its `n` leading elements by mass, as "O 46%, Fe 19%, …".
    #[must_use]
    pub fn leading(&self, n: usize) -> String {
        let mut by = self.makeup.clone();
        by.sort_by(|a, b| b.1.total_cmp(&a.1));
        by.iter().take(n).map(|(e, w)| format!("{} {:.0}%", e.symbol(), 100.0 * w)).collect::<Vec<_>>().join(", ")
    }

    /// What to show of it, for a body of `mass_kg`: type, elements, metal, density, family.
    #[must_use]
    pub fn facts(&self, mass_kg: f64) -> Vec<(&'static str, String)> {
        let tonnes = |share: f64| format!("{:.1e} t", share * mass_kg / 1000.0);
        let mut out = vec![("Type", self.taxon.name().to_string()), ("Makeup", self.leading(5))];
        if self.taxon.metal_share() > 0.0 {
            out.push(("Metal", format!("{:.0}% iron–nickel, {}", 100.0 * self.taxon.metal_share(), tonnes(self.taxon.metal_share()))));
        }
        out.push(("Platinum group", format!("{:.2} ppm, {}", 1e6 * self.share(Element::Pgm), tonnes(self.share(Element::Pgm)))));
        out.push(("Density", format!("{:.2} g/cm³, {:.0}% empty", self.density_g_cc, 100.0 * self.porosity)));
        if let Some(f) = &self.family {
            out.push(("Family", format!("{} family ({})", f.name, f.taxon.name())));
        }
        out
    }
}

/// The share of each type at `x` = distance / frost line (the Sun's belt from 0.78 to 1.22,
/// its Hildas 1.45, its Trojans 1.93), [S, C, M, V, P, D].
#[must_use]
pub fn shares(x: f64) -> [f64; 6] {
    const AT: [(f64, [f64; 6]); 7] = [
        (0.70, [0.65, 0.15, 0.08, 0.08, 0.02, 0.02]),
        (0.85, [0.50, 0.25, 0.08, 0.07, 0.05, 0.05]),
        (1.00, [0.30, 0.45, 0.07, 0.03, 0.10, 0.05]),
        (1.15, [0.15, 0.55, 0.05, 0.01, 0.14, 0.10]),
        (1.30, [0.08, 0.40, 0.04, 0.00, 0.28, 0.20]),
        (1.50, [0.05, 0.15, 0.03, 0.00, 0.37, 0.40]),
        (1.90, [0.02, 0.05, 0.01, 0.00, 0.25, 0.67]),
    ];
    if x <= AT[0].0 {
        return AT[0].1;
    }
    for w in AT.windows(2) {
        let ((x0, a), (x1, b)) = (w[0], w[1]);
        if x <= x1 {
            let t = (x - x0) / (x1 - x0);
            return [0, 1, 2, 3, 4, 5].map(|k| a[k] + t * (b[k] - a[k]));
        }
    }
    AT[AT.len() - 1].1
}

/// A type drawn at `x` (distance / frost line), icy in an icy belt.
fn draw_taxon(d: &mut Draw, x: f64, icy: bool) -> Taxon {
    if icy {
        return Taxon::Icy;
    }
    let s = shares(x);
    let mut u = d.u() * s.iter().sum::<f64>();
    for (k, w) in s.iter().enumerate() {
        u -= w;
        if u <= 0.0 {
            return Taxon::ALL[k];
        }
    }
    Taxon::D
}

/// Share of the volume empty in a body `diameter_km` across: rubble piles 20–50% (a comet's
/// 40–75%), bodies of a few hundred km less, the largest and the smallest (monoliths) little.
fn porosity(d: &mut Draw, diameter_km: f64, icy: bool) -> f64 {
    match (diameter_km, icy) {
        (x, _) if x < 0.2 => d.range(0.0, 0.15),
        (x, false) if x < 100.0 => d.range(0.2, 0.5),
        (x, true) if x < 100.0 => d.range(0.4, 0.75),
        (x, false) if x < 400.0 => d.range(0.05, 0.3),
        (x, true) if x < 1000.0 => d.range(0.15, 0.5),
        _ => d.range(0.0, 0.08),
    }
}

/// The families of belt `pop` (seeded `seed`) at frost line `frost_au`.
#[must_use]
pub fn families(seed: Seed, pop: &Population, frost_au: f64, icy: bool) -> Vec<Family> {
    let mut d = Draw::new(seed.child(b"families"));
    let n = d.int(5, 30);
    (0..n)
        .map(|index| {
            let mut f = Draw::new(seed.child_u64(b"family", u64::from(index)));
            let a_au = f.range(pop.inner_au, pop.outer_au);
            let taxon = draw_taxon(&mut f, relative(a_au, frost_au), icy);
            Family { index, name: geoforge_names::name(seed.child_u64(b"family name", u64::from(index)).0), taxon, a_au, e: f.range(0.02, 0.2), i_deg: f.range(1.0, 15.0) }
        })
        .collect()
}

/// The make of member `index` of `pop` (seeded `seed`), `diameter_km` across, its orbit
/// `orbit` as drawn: its type, elements, porosity and density, and its family, whose orbit
/// it then takes (`orbit` moved onto it); `frost_au` is its star's frost line.
#[must_use]
pub fn member(seed: Seed, pop: &Population, index: u64, diameter_km: f64, frost_au: f64, orbit: &mut Elements) -> Asteroid {
    let icy = pop.composition.fraction(Species::Ice) > 0.3;
    let mut d = Draw::new(seed.child_u64(b"asteroid", index));
    let family = (d.u() < FAMILY_SHARE).then(|| {
        let all = families(seed, pop, frost_au, icy);
        all[(d.u() * all.len() as f64) as usize % all.len()].clone()
    });
    if let Some(f) = &family {
        orbit.a_au = (f.a_au + d.range(-0.01, 0.01)).clamp(pop.inner_au, pop.outer_au);
        orbit.e = (f.e + d.range(-0.02, 0.02)).max(0.0);
        orbit.i_deg = (f.i_deg + d.range(-0.5, 0.5)).max(0.0);
    }
    // A family's fragments share their parent's type, but for a tenth (a differentiated
    // parent's mantle and core, or interlopers).
    let taxon = match &family {
        Some(f) if d.u() < 0.9 => f.taxon,
        _ => draw_taxon(&mut d, relative(orbit.a_au, frost_au), icy),
    };
    let raw: Vec<(Element, f64)> = taxon.analog().iter().map(|&(e, w)| (e, w * d.range(0.85, 1.15))).collect();
    let total: f64 = raw.iter().map(|x| x.1).sum();
    let makeup = raw.into_iter().map(|(e, w)| (e, w / total)).collect();
    let porosity = porosity(&mut d, diameter_km, taxon == Taxon::Icy);
    Asteroid { taxon, makeup, porosity, density_g_cc: taxon.grain_density() * (1.0 - porosity), family }
}

/// Distance as a share of the frost line (0 round a star with none, a neutron star's).
fn relative(a_au: f64, frost_au: f64) -> f64 {
    if frost_au.is_finite() && frost_au > 0.0 {
        a_au / frost_au
    } else {
        0.0
    }
}

/// The make of the belt member `key` of `system` (as [`crate::context::belt_member`] keys
/// them), `None` for any other body.
#[must_use]
pub fn of(system: &SolarSystem, key: Key) -> Option<Asteroid> {
    let belt = key.0 >> 32;
    let index = (key.0 & 0xFFFF_FFFF).checked_sub(1)?;
    if belt == 0 {
        return None;
    }
    crate::context::belt_member_made(system, Key(belt), index).map(|(_, made)| made)
}

/// The share of belt `pop`'s members of the types `kind` admits, at frost line `frost_au`
/// (the shares by distance averaged over the belt).
#[must_use]
pub fn kind_share(pop: &Population, frost_au: f64, kind: impl Fn(Taxon) -> bool) -> f64 {
    if pop.composition.fraction(Species::Ice) > 0.3 {
        return f64::from(u8::from(kind(Taxon::Icy)));
    }
    const STEPS: usize = 16;
    (0..STEPS)
        .map(|k| {
            let a = pop.inner_au + (pop.outer_au - pop.inner_au) * (k as f64 + 0.5) / STEPS as f64;
            let s = shares(relative(a, frost_au));
            Taxon::ALL[..6].iter().zip(s).filter(|(t, _)| kind(**t)).map(|(_, w)| w).sum::<f64>()
        })
        .sum::<f64>()
        / STEPS as f64
}

/// The mix of types in belt `belt` of `system`, from its `sample` largest members past the
/// first hundred (whose few would skew it): shares, largest first.
#[must_use]
pub fn mix(system: &SolarSystem, belt: Key, sample: u64) -> Vec<(Taxon, f64)> {
    let mut count = [0usize; 7];
    let mut n = 0;
    for i in 100..100 + sample {
        let Some((_, a)) = crate::context::belt_member_made(system, belt, i) else { break };
        count[Taxon::ALL.iter().position(|&t| t == a.taxon).unwrap_or(0)] += 1;
        n += 1;
    }
    let mut out: Vec<(Taxon, f64)> = Taxon::ALL.iter().zip(count).filter(|x| x.1 > 0).map(|(&t, c)| (t, c as f64 / n as f64)).collect();
    out.sort_by(|a, b| b.1.total_cmp(&a.1));
    out
}

/// [`mix`] as a line: "C 45%, S 30%, P 12%, …" (the types over 2%).
#[must_use]
pub fn mix_line(system: &SolarSystem, belt: Key) -> Option<String> {
    let m = mix(system, belt, 400);
    (!m.is_empty()).then(|| m.iter().filter(|x| x.1 >= 0.02).map(|(t, w)| format!("{} {:.0}%", t.code(), 100.0 * w)).collect::<Vec<_>>().join(", "))
}

#[cfg(test)]
mod tests {
    use super::*;
    use geoforge_core::body::{Composition, PowerLaw};

    /// The Sun's frost line, AU.
    const SOL_FROST_AU: f64 = 2.7;

    fn belt(ice: f64) -> Population {
        Population {
            inner_au: 2.1,
            outer_au: 3.3,
            count_estimate: 1e6,
            size_law: PowerLaw { min_km: 1.0, max_km: 900.0, exponent: 3.5, knee: None },
            composition: Composition::from_weights(&[(Species::Silicate, 1.0 - ice), (Species::Ice, ice)]),
        }
    }

    /// The Sun's belt: S-types lead its inner part, C-types its outer, D-types its Trojans;
    /// shares sum to one everywhere.
    #[test]
    fn types_by_distance() {
        let s = |a: f64| shares(a / SOL_FROST_AU);
        assert!(s(2.2)[0] > s(2.2)[1] && s(3.1)[1] > s(3.1)[0] && s(5.2)[5] > 0.6);
        for x in [0.5, 0.9, 1.1, 1.4, 2.5] {
            assert!((shares(x).iter().sum::<f64>() - 1.0).abs() < 1e-9);
        }
    }

    /// Members: makeups sum to one, iron for M-types, water's hydrogen in C-types, ice in an
    /// icy belt; rubble piles are porous; a third are in families and keep near their orbit;
    /// the same seed makes the same body.
    #[test]
    fn members_are_made_by_type() {
        let pop = belt(0.0);
        let (mut in_family, mut seen) = (0, std::collections::BTreeSet::new());
        for i in 0..600 {
            let mut orbit = Elements { a_au: 2.5, e: 0.1, i_deg: 5.0, lan_deg: 0.0, argp_deg: 0.0, m0_deg: 0.0 };
            let a = member(Seed(3), &pop, i, 10.0, 2.7, &mut orbit);
            assert!((a.makeup.iter().map(|x| x.1).sum::<f64>() - 1.0).abs() < 1e-9);
            let share = |e: Element| a.makeup.iter().find(|x| x.0 == e).map_or(0.0, |x| x.1);
            match a.taxon {
                Taxon::M => assert!(share(Element::Fe) > 0.8 && share(Element::Pgm) > 1e-5),
                Taxon::C => assert!(share(Element::H) > 0.01 && share(Element::C) > 0.02),
                _ => {}
            }
            assert!((0.2..0.5).contains(&a.porosity) && a.density_g_cc < a.taxon.grain_density());
            if let Some(f) = &a.family {
                in_family += 1;
                assert!((orbit.a_au - f.a_au).abs() <= 0.0101 || orbit.a_au == pop.inner_au || orbit.a_au == pop.outer_au);
            }
            seen.insert(a.taxon);
            let mut again = Elements { a_au: 2.5, e: 0.1, i_deg: 5.0, lan_deg: 0.0, argp_deg: 0.0, m0_deg: 0.0 };
            assert_eq!(a, member(Seed(3), &pop, i, 10.0, 2.7, &mut again));
        }
        assert!((150..270).contains(&in_family), "{in_family} of 600 in families");
        assert!(seen.len() >= 5, "{seen:?}");
        let mut orbit = Elements { a_au: 40.0, e: 0.1, i_deg: 5.0, lan_deg: 0.0, argp_deg: 0.0, m0_deg: 0.0 };
        assert_eq!(member(Seed(3), &belt(0.6), 0, 300.0, 2.7, &mut orbit).taxon, Taxon::Icy);
    }

    /// The Sun's main belt as made: stony S-types lead inside 2.5 AU, carbonaceous C-types
    /// beyond 2.8; a few per cent metal; the Kuiper belt icy.
    #[test]
    fn sols_belts() {
        let sol = crate::sol::at(Seed(1), 2_451_545.0);
        let main = sol.by_name("Main Belt").unwrap().key;
        let (mut inner, mut outer, mut metal, mut n) = ([0usize; 7], [0usize; 7], 0, 0);
        for i in 200..2200 {
            let (ctx, a) = crate::context::belt_member_made(&sol, main, i).unwrap();
            let k = Taxon::ALL.iter().position(|&t| t == a.taxon).unwrap();
            if ctx.orbit.a_au < 2.5 {
                inner[k] += 1;
            } else if ctx.orbit.a_au > 2.8 {
                outer[k] += 1;
            }
            metal += usize::from(a.taxon == Taxon::M);
            n += 1;
            assert_eq!(of(&sol, ctx.key), Some(a));
        }
        assert!(inner[0] > inner[1] && outer[1] > outer[0], "inner {inner:?} outer {outer:?}");
        assert!((0.02..0.12).contains(&(metal as f64 / n as f64)), "{metal} metal of {n}");
        let kuiper = sol.by_name("Kuiper Belt").unwrap().key;
        assert_eq!(crate::context::belt_member_made(&sol, kuiper, 3).unwrap().1.taxon, Taxon::Icy);
        // About ten Kuiper belt objects over 1000 km, as the real belt has.
        let size = |i: u64| 2.0 * crate::context::belt_member(&sol, kuiper, i).unwrap().exterior.radius_km;
        assert!(size(8) > 1000.0 && size(12) < 1000.0 && size(0) > 1500.0, "{} {} {}", size(0), size(8), size(12));
    }
}
