//! Nebulae: the gas the young stars light and the dead leave, found where the galaxy's stars
//! and clouds are (so each has a star or a cloud behind it, and none is stored).
//!
//! - **Emission** (H II regions): a birth event younger than 10 Myr whose living members
//!   send out ≥ 10⁴⁷ ionising photons a second (an O or early B star; Martins et al. 2005's
//!   rates by mass). Its radius is the Strömgren sphere's in gas of 100 cm⁻³, grown as the
//!   hot gas pushes out at 10 km/s (Spitzer 1978), up to 60 pc; its light is Orion's
//!   (V 4.0 at 412 pc for 10⁴⁹ s⁻¹) scaled by the photons.
//! - **Reflection**: a birth event younger than 100 Myr with no H II region and a star of
//!   2.5 M☉ or more, still in its natal gas (under 10 Myr) or passing through a molecular
//!   cloud (the Pleiades); 0.5 pc × √(L / 10³ L☉) across, two magnitudes fainter than the star.
//! - **Dark**: a molecular cloud of 10⁴ M☉ or more ([`crate::sectors::dust`]), out to 2σ.
//! - **Planetary**: a white dwarf born within 25 kyr, its shell thrown at 20 km/s; the
//!   Helix's light (V 7.6 at 200 pc at 10 kyr), fading 2.5 magnitudes over its life.
//! - **Supernova remnants**: a neutron star or black hole born within 100 kyr: freely
//!   expanding at 5,000 km/s, then a Sedov–Taylor blast (10⁵¹ erg in 1 cm⁻³), the Crab's
//!   light at a thousand years (V 8.4 at 2 kpc), fading with age.
//!
//! Brightnesses are integrated visual magnitudes, for drawing, good to a magnitude or two.

use crate::sectors::bins::{self, NASCENT_COLLAPSED, NASCENT_WHITE_DWARFS, PLANETARY_GYR, REMNANT_GYR};
use crate::sectors::groups::{self, LEVEL};
use crate::sectors::{GalacticStar, Sectors};
use crate::stars::end_of_life_gyr;
use geoforge_core::body::StellarStage;
use geoforge_core::galactic::StarId;
use geoforge_core::seed::Seed;
use rayon::prelude::*;
use std::fmt;

/// What kind of nebula.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum NebulaKind {
    Emission,
    Reflection,
    Dark,
    Planetary,
    SupernovaRemnant,
}

impl NebulaKind {
    pub const ALL: [Self; 5] = [Self::Emission, Self::Reflection, Self::Dark, Self::Planetary, Self::SupernovaRemnant];

    /// Its name as [`geoforge_core::exotics::NEBULA_KINDS`] has it.
    #[must_use]
    pub fn code(self) -> &'static str {
        match self {
            Self::Emission => "emission",
            Self::Reflection => "reflection",
            Self::Dark => "dark",
            Self::Planetary => "planetary",
            Self::SupernovaRemnant => "supernova_remnant",
        }
    }

    #[must_use]
    pub fn label(self) -> &'static str {
        match self {
            Self::Emission => "emission nebula (H II region)",
            Self::Reflection => "reflection nebula",
            Self::Dark => "dark nebula (molecular cloud)",
            Self::Planetary => "planetary nebula",
            Self::SupernovaRemnant => "supernova remnant",
        }
    }

    /// What follows its name.
    fn suffix(self) -> &'static str {
        match self {
            Self::Emission => "Nebula",
            Self::Reflection => "Reflection Nebula",
            Self::Dark => "Cloud",
            Self::Planetary => "Planetary Nebula",
            Self::SupernovaRemnant => "Remnant",
        }
    }
}

/// What a nebula is made by, which names it: a birth event's gas, a cloud, or a star's death.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum NebulaId {
    /// Birth event `event` of group cube `cube`.
    Event { cube: [i32; 3], event: u32 },
    /// The `n`th molecular cloud of group cube `cube`.
    Cloud { cube: [i32; 3], n: u32 },
    /// The nebula a just-dead star left.
    Star(StarId),
}

impl fmt::Display for NebulaId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Event { cube: [i, j, k], event } => write!(f, "e{i},{j},{k}:{event}"),
            Self::Cloud { cube: [i, j, k], n } => write!(f, "c{i},{j},{k}:{n}"),
            Self::Star(s) => write!(f, "s{s}"),
        }
    }
}

/// A nebula.
#[derive(Debug, Clone, PartialEq)]
pub struct Nebula {
    pub id: NebulaId,
    pub kind: NebulaKind,
    pub name: String,
    /// pc from the galaxy's centre.
    pub centre_pc: [f64; 3],
    pub radius_pc: f64,
    /// Myr since its stars formed or its star died; `None` for a cloud.
    pub age_myr: Option<f64>,
    /// Gas, M☉ (a cloud's).
    pub mass_msun: Option<f64>,
    /// Ionising photons a second (an H II region's).
    pub ionising_per_s: Option<f64>,
    /// The star that lights it or left it, and that star's stage (the one that left it,
    /// not a companion).
    pub star: Option<StarId>,
    pub star_stage: Option<StellarStage>,
    /// Integrated visual magnitude at 10 pc, before dust; `None` for a dark cloud.
    pub absolute_v: Option<f64>,
    /// Its own seed (for its name and anything drawn for it).
    pub seed: Seed,
}

impl Nebula {
    /// Distance from `p`, pc.
    #[must_use]
    pub fn distance_pc(&self, p: [f64; 3]) -> f64 {
        (0..3).map(|k| (self.centre_pc[k] - p[k]).powi(2)).sum::<f64>().sqrt()
    }

    /// Its apparent magnitude from `p`, before dust.
    #[must_use]
    pub fn apparent_v(&self, p: [f64; 3]) -> Option<f64> {
        self.absolute_v.map(|m| m + 5.0 * (self.distance_pc(p).max(self.radius_pc).max(0.1) / 10.0).log10())
    }

    /// Its angular radius from `p`, degrees.
    #[must_use]
    pub fn angular_radius_deg(&self, p: [f64; 3]) -> f64 {
        let d = self.distance_pc(p);
        if d <= self.radius_pc {
            180.0
        } else {
            (self.radius_pc / d).asin().to_degrees()
        }
    }

    /// Whether `p` lies inside it.
    #[must_use]
    pub fn contains(&self, p: [f64; 3]) -> bool {
        self.distance_pc(p) <= self.radius_pc
    }
}

/// Ionising photons a second, by initial mass (M☉): log₁₀ Q for O and early B dwarfs
/// (Martins et al. 2005, Smith et al. 2002).
const IONISING: [(f64, f64); 12] =
    [(9.0, 44.5), (12.0, 45.8), (14.0, 47.4), (16.0, 47.9), (19.0, 48.3), (22.0, 48.6), (25.0, 48.9), (31.0, 49.2), (37.0, 49.4), (58.0, 49.6), (100.0, 49.9), (150.0, 50.1)];
/// Events younger than this (Gyr) are still in their natal gas.
const NATAL_GYR: f64 = 0.01;
/// Events younger than this (Gyr) may still light a reflection nebula.
const REFLECTION_GYR: f64 = 0.1;
/// The least ionising photons a second of an H II region.
const LEAST_IONISING: f64 = 1e47;
/// The least gas, M☉, of a cloud counted as a dark nebula.
const DARK_MASS: f64 = 1e4;
/// Bins whose stars can ionise or light gas: brighter than absolute V −2.
const BRIGHT_BINS: usize = 10;

/// Ionising photons a second from a living star born with `m` M☉.
#[must_use]
pub fn ionising_per_s(m: f64) -> f64 {
    if m < IONISING[0].0 {
        return 0.0;
    }
    let lm = m.ln();
    let log_q = IONISING.windows(2).find(|w| m <= w[1].0).map_or(IONISING[IONISING.len() - 1].1, |w| w[0].1 + (w[1].1 - w[0].1) * (lm - w[0].0.ln()) / (w[1].0.ln() - w[0].0.ln()));
    10f64.powf(log_q)
}

/// An H II region's radius, pc, `age_myr` after its stars lit with `q` photons a second: the
/// Strömgren radius in 100 cm⁻³ (3.15 pc at 10⁴⁹ s⁻¹), pushed out at 10 km/s.
#[must_use]
pub fn h_ii_radius_pc(q: f64, age_myr: f64) -> f64 {
    let rs = 3.15 * (q / 1e49).cbrt();
    (rs * (1.0 + 7.0 * 10.23 * age_myr / (4.0 * rs)).powf(4.0 / 7.0)).min(60.0)
}

/// A supernova remnant's radius, pc, `age_yr` after the blast.
#[must_use]
pub fn remnant_radius_pc(age_yr: f64) -> f64 {
    (5.1e-3 * age_yr).min(0.31 * age_yr.max(1.0).powf(0.4))
}

/// The nebulae within `radius_pc` of `p` (their centres), nearest first.
#[must_use]
pub fn near(sectors: &Sectors, p: [f64; 3], radius_pc: f64) -> Vec<Nebula> {
    let base = sectors.seed().child(b"nebula");
    let named = |id: NebulaId, kind: NebulaKind| {
        let seed = base.child(id.to_string().as_bytes());
        (seed, format!("{} {}", geoforge_names::name(seed.child(b"name").0), kind.suffix()))
    };
    // Birth events and clouds, by group cube.
    let cubes: Vec<[i32; 3]> = bins::cubes_near(p, radius_pc + groups::MAX_SPREAD_PC, LEVEL).collect();
    let mut out: Vec<Nebula> = cubes
        .par_iter()
        .flat_map_iter(|&cube| {
            let mut found = Vec::new();
            let clouds = sectors.clouds(cube);
            for (n, c) in clouds.iter().enumerate() {
                let id = NebulaId::Cloud { cube, n: n as u32 };
                if c.mass_msun >= DARK_MASS && dist(c.centre, p) <= radius_pc {
                    let (seed, name) = named(id, NebulaKind::Dark);
                    found.push(Nebula {
                        id,
                        kind: NebulaKind::Dark,
                        name,
                        centre_pc: c.centre,
                        radius_pc: 2.0 * c.sigma_pc,
                        age_myr: None,
                        mass_msun: Some(c.mass_msun),
                        ionising_per_s: None,
                        star: None,
                        star_stage: None,
                        absolute_v: None,
                        seed,
                    });
                }
            }
            for e in sectors.events(cube) {
                let near_p = dist(e.centre, p) <= radius_pc + e.extent_pc();
                if !near_p {
                    continue;
                }
                // The just dead.
                for bin in [NASCENT_WHITE_DWARFS, NASCENT_COLLAPSED] {
                    found.extend(sectors.members(cube, &e, bin).into_iter().filter(|s| dist(s.position_pc, p) <= radius_pc).filter_map(|s| left_by(&s, &named)));
                }
                if e.age_gyr >= REFLECTION_GYR || dist(e.centre, p) > radius_pc {
                    continue;
                }
                let lit: Vec<GalacticStar> = (0..=BRIGHT_BINS).flat_map(|b| sectors.members(cube, &e, b)).filter(|s| !s.star.stage.is_remnant()).collect();
                let id = NebulaId::Event { cube, event: e.index };
                let age_myr = e.age_gyr * 1000.0;
                let q: f64 = lit.iter().map(|s| ionising_per_s(s.kind.initial_mass)).sum();
                if e.age_gyr < NATAL_GYR && q >= LEAST_IONISING {
                    let brightest = lit.iter().max_by(|a, b| ionising_per_s(a.kind.initial_mass).total_cmp(&ionising_per_s(b.kind.initial_mass)));
                    let (seed, name) = named(id, NebulaKind::Emission);
                    found.push(Nebula {
                        id,
                        kind: NebulaKind::Emission,
                        name,
                        centre_pc: e.centre,
                        radius_pc: h_ii_radius_pc(q, age_myr),
                        age_myr: Some(age_myr),
                        mass_msun: None,
                        ionising_per_s: Some(q),
                        star: brightest.map(|s| s.id),
                        star_stage: brightest.map(|s| s.star.stage),
                        absolute_v: Some(-4.1 - 2.5 * (q / 1e49).log10()),
                        seed,
                    });
                    continue;
                }
                let in_gas = e.age_gyr < NATAL_GYR || clouds.iter().any(|c| dist(c.centre, e.centre) <= 2.0 * c.sigma_pc);
                let Some(star) = lit.iter().filter(|s| s.kind.initial_mass >= 2.5).max_by(|a, b| a.star.luminosity_solar.total_cmp(&b.star.luminosity_solar)) else { continue };
                if in_gas {
                    let (seed, name) = named(id, NebulaKind::Reflection);
                    found.push(Nebula {
                        id,
                        kind: NebulaKind::Reflection,
                        name,
                        centre_pc: star.position_pc,
                        radius_pc: (0.5 * (star.star.luminosity_solar / 1e3).sqrt()).clamp(0.2, 5.0),
                        age_myr: Some(age_myr),
                        mass_msun: None,
                        ionising_per_s: None,
                        star: Some(star.id),
                        star_stage: Some(star.star.stage),
                        absolute_v: star.absolute_v.is_finite().then_some(star.absolute_v + 2.0),
                        seed,
                    });
                }
            }
            found
        })
        .collect();
    // The field's just dead.
    let level = bins::level(NASCENT_COLLAPSED);
    let field: Vec<Nebula> = bins::cubes_near(p, radius_pc, level)
        .collect::<Vec<_>>()
        .par_iter()
        .flat_map_iter(|&cube| {
            sectors
                .field_cube(level, cube, &[NASCENT_WHITE_DWARFS, NASCENT_COLLAPSED])
                .into_iter()
                .filter(|s| dist(s.position_pc, p) <= radius_pc)
                .filter_map(|s| left_by(&s, &named))
                .collect::<Vec<_>>()
        })
        .collect();
    out.extend(field);
    // A reflection nebula sits round its star, which may lie off its group's centre.
    out.retain(|n| n.distance_pc(p) <= radius_pc);
    out.sort_by(|a, b| a.distance_pc(p).total_cmp(&b.distance_pc(p)).then_with(|| a.id.to_string().cmp(&b.id.to_string())));
    out
}

/// A nebula as seen from a place.
#[derive(Debug, Clone, PartialEq)]
pub struct SeenNebula {
    pub nebula: Nebula,
    /// Unit vector toward it in the system's frame.
    pub direction: [f64; 3],
    pub distance_pc: f64,
    /// Its integrated magnitude seen through the dust; `None` for a dark cloud.
    pub v: Option<f64>,
    pub radius_deg: f64,
}

/// The nebulae `observer` sees within `radius_pc`: the lit ones brighter than `limit`
/// through the dust, and the dark clouds wider than a degree.
#[must_use]
pub fn seen(field: &crate::starfield::Starfield, observer: &crate::sky::Observer, radius_pc: f64, limit: f64) -> Vec<SeenNebula> {
    let here = observer.position_pc;
    let dust = crate::sectors::dust::DustView::new(field.sectors(), here, radius_pc);
    near(field.sectors(), here, radius_pc)
        .into_iter()
        .filter_map(|n| {
            let d = n.distance_pc(here);
            let radius_deg = n.angular_radius_deg(here);
            let v = n.apparent_v(here).map(|m| m + if d > n.radius_pc { dust.to(n.centre_pc) } else { 0.0 });
            let shown = match v {
                Some(v) => v <= limit,
                None => radius_deg >= 1.0,
            };
            let g = [0, 1, 2].map(|k| (n.centre_pc[k] - here[k]) / d.max(1e-9));
            shown.then(|| SeenNebula { direction: observer.to_system_frame(g), distance_pc: d, v, radius_deg, nebula: n })
        })
        .collect()
}

/// The nebula a just-dead star left: a planetary nebula round a white dwarf, a remnant round
/// a neutron star or black hole.
fn left_by(s: &GalacticStar, named: &(dyn Fn(NebulaId, NebulaKind) -> (Seed, String) + Sync)) -> Option<Nebula> {
    let dead_gyr = s.kind.age_gyr - end_of_life_gyr(s.kind.initial_mass);
    let primary = crate::stars::evolve(s.kind.initial_mass, s.kind.age_gyr, s.kind.metallicity_dex);
    let id = NebulaId::Star(s.id);
    let (kind, radius_pc, absolute_v) = match primary.stage {
        StellarStage::WhiteDwarf if (0.0..=PLANETARY_GYR).contains(&dead_gyr) => (NebulaKind::Planetary, 20.4 * dead_gyr * 1000.0 + 0.01, 0.1 + 2.5 * dead_gyr / PLANETARY_GYR),
        StellarStage::NeutronStar | StellarStage::BlackHole if (0.0..=REMNANT_GYR).contains(&dead_gyr) => {
            let yr = dead_gyr * 1e9;
            (NebulaKind::SupernovaRemnant, remnant_radius_pc(yr), -3.85 + 2.5 * (1.0 + yr / 1000.0).log10())
        }
        _ => return None,
    };
    let (seed, name) = named(id, kind);
    Some(Nebula {
        id,
        kind,
        name,
        centre_pc: s.position_pc,
        radius_pc,
        age_myr: Some(dead_gyr * 1000.0),
        mass_msun: None,
        ionising_per_s: None,
        star: Some(s.id),
        star_stage: Some(primary.stage),
        absolute_v: Some(absolute_v),
        seed,
    })
}

fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|k| (a[k] - b[k]).powi(2)).sum::<f64>().sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::galaxy::Galaxy;

    /// Around the Sun's place in a Milky Way: H II regions, reflection nebulae and dark clouds
    /// within 1.5 kpc, a few tens of planetary nebulae and some supernova remnants within
    /// 1 kpc (the local rates: ~600 white dwarfs and ~25 core collapses per kpc² per Myr); each
    /// the size and age its kind allows, its star the kind that makes it; the same again.
    #[test]
    fn nebulae_round_the_sun() {
        let s = Sectors::new(&Galaxy::milky_way(), Seed(5));
        let sun = [s.galaxy().reference_kpc * 1000.0, 0.0, 20.0];
        let all = near(&s, sun, 1500.0);
        let count = |k: NebulaKind, within: f64| all.iter().filter(|n| n.kind == k && n.distance_pc(sun) <= within).count();
        let (h2, refl, dark, pn, snr) = (
            count(NebulaKind::Emission, 1500.0),
            count(NebulaKind::Reflection, 1500.0),
            count(NebulaKind::Dark, 1500.0),
            count(NebulaKind::Planetary, 1000.0),
            count(NebulaKind::SupernovaRemnant, 1000.0),
        );
        assert!((3..300).contains(&h2) && refl >= 3 && dark >= 10, "{h2} H II, {refl} reflection, {dark} dark");
        assert!((10..400).contains(&pn) && (1..80).contains(&snr), "{pn} planetary, {snr} remnants");
        for n in &all {
            assert!(n.radius_pc > 0.0 && n.radius_pc <= 110.0, "{n:?}");
            let star = n.star.and_then(|id| s.star(&id));
            // The star that made it (not a companion that may outshine it).
            let primary = star.as_ref().map(|st| crate::stars::evolve(st.kind.initial_mass, st.kind.age_gyr, st.kind.metallicity_dex).stage);
            match n.kind {
                NebulaKind::Emission => assert!(n.ionising_per_s.unwrap() >= LEAST_IONISING && n.age_myr.unwrap() < 10.0 && star.is_some()),
                NebulaKind::Planetary => assert!(n.radius_pc < 0.6 && primary == Some(StellarStage::WhiteDwarf)),
                NebulaKind::SupernovaRemnant => assert!(n.age_myr.unwrap() <= 0.1 && matches!(primary, Some(StellarStage::NeutronStar | StellarStage::BlackHole))),
                NebulaKind::Dark => assert!(n.mass_msun.unwrap() >= DARK_MASS),
                NebulaKind::Reflection => assert!(star.is_some()),
            }
        }
        assert_eq!(all, near(&s, sun, 1500.0));
    }

    /// Orion-like numbers: 10⁴⁹ photons a second light 3 pc at once and some 16 pc after
    /// 3 Myr; a remnant is 5 pc across at a thousand years, 31 pc at a hundred thousand.
    #[test]
    fn sizes() {
        assert!((h_ii_radius_pc(1e49, 0.0) - 3.15).abs() < 0.01);
        assert!((14.0..18.0).contains(&h_ii_radius_pc(1e49, 3.0)));
        assert!((ionising_per_s(37.0).log10() - 49.4).abs() < 1e-9 && ionising_per_s(5.0) == 0.0);
        assert!((4.5..5.2).contains(&remnant_radius_pc(1000.0)) && (30.0..32.0).contains(&remnant_radius_pc(1e5)));
    }
}
