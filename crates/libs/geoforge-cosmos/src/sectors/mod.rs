//! Stars that live in a galaxy, generated where they are looked at.
//!
//! Space in the galaxy's frame is cut into nested cubes ([`bins`]): each brightness bin's
//! stars live in cubes about half as wide as the bin is seen from, so a sky reads its bright
//! stars from a few big cubes and its faint ones only from small cubes nearby, and zooming
//! in generates fainter stars only where it looks. A cube's stars are a pure function of the
//! galaxy's seed and the cube's key ([`field`]): a Poisson count on the components' mass in
//! it times their stars per M☉ in the bin ([`tables`]), placed by the density inside it.
//! Every star has an identity ([`StarId`]) that seeds the system made at it.
//!
//! Stars younger than [`YOUNG_GYR`] are born in groups that spread with age ([`groups`]),
//! but for the quarter born in groups too small to tell apart, which are young field stars.
//! Dust is a smooth layer plus clouds where stars form ([`dust`]).

pub mod bins;
pub mod dust;
pub mod field;
pub mod groups;
pub mod tables;

use crate::galaxy::{ComponentKind, Galaxy, Height, Prepared, Shape};
use crate::sky::absolute_v;
use crate::stars::evolve;
use bins::BINS;
use geoforge_core::body::Star;
use geoforge_core::galactic::{StarId, StarKind};
use geoforge_core::seed::Seed;
use rayon::prelude::*;
use tables::Table;

/// Stars younger than this are born in groups, Gyr.
pub const YOUNG_GYR: f64 = 0.3;
/// Age slices of a disk whose stars thicken with age, Gyr (from [`YOUNG_GYR`]).
const DISK_SLICES: [f64; 6] = [0.6, 1.2, 2.5, 5.0, 8.0, 1e9];
/// Kinds sampled per field population.
const FIELD_KINDS: usize = 24_000;

/// A star of a galaxy, generated or edited.
#[derive(Debug, Clone, PartialEq)]
pub struct GalacticStar {
    pub id: StarId,
    /// pc from the galaxy's centre, z toward its north pole.
    pub position_pc: [f64; 3],
    pub kind: StarKind,
    /// The system's brighter star as it is now (`kind` evolved): the companion, where the
    /// primary has died and it still shines.
    pub star: Star,
    /// Absolute V of the system's stars together; infinite for a dark remnant.
    pub absolute_v: f64,
    /// The component it was born in; `None` for a star added by hand.
    pub component: Option<ComponentKind>,
}

impl GalacticStar {
    /// The star of `kind` at `position_pc`.
    #[must_use]
    pub fn new(id: StarId, position_pc: [f64; 3], kind: StarKind, component: Option<ComponentKind>) -> Self {
        let mut star = evolve(kind.initial_mass, kind.age_gyr, kind.metallicity_dex);
        let mut v = absolute_v(&star);
        if let Some(mc) = kind.companion_mass {
            let mate = evolve(mc, kind.age_gyr, kind.metallicity_dex);
            if let Some(mv) = absolute_v(&mate) {
                if v.is_none_or(|pv| mv < pv) {
                    star = mate;
                }
                v = Some(v.map_or(mv, |pv| -2.5 * (10f64.powf(-0.4 * pv) + 10f64.powf(-0.4 * mv)).log10()));
            }
        }
        Self { id, position_pc, kind, star, absolute_v: v.unwrap_or(f64::INFINITY), component }
    }

    /// Distance to `p`, pc.
    #[must_use]
    pub fn distance_pc(&self, p: [f64; 3]) -> f64 {
        (0..3).map(|k| (self.position_pc[k] - p[k]).powi(2)).sum::<f64>().sqrt()
    }
}

/// A slice of a component by age, with its stars per bin.
pub(crate) struct Population {
    pub component: usize,
    /// The age its shape and arms are taken at, Gyr.
    pub age: f64,
    /// Its stellar mass now, M☉.
    pub mass_msun: f64,
    pub shape: Prepared,
    pub table: Table,
}

/// A galaxy's stars, ready to generate anywhere.
pub struct Sectors {
    galaxy: Galaxy,
    seed: Seed,
    pub(crate) field: Vec<Population>,
    pub(crate) young: Vec<groups::YoungPopulation>,
}

impl Sectors {
    /// The stars of `galaxy`, seeded `seed`.
    #[must_use]
    pub fn new(galaxy: &Galaxy, seed: Seed) -> Self {
        let mut field = Vec::new();
        let mut young = Vec::new();
        for (c, comp) in galaxy.components.iter().enumerate() {
            let h = &comp.history;
            let (y, a) = (h.youngest_gyr, h.age_gyr);
            let first = y.max(YOUNG_GYR);
            let mut edges = vec![first];
            if matches!(comp.shape, Shape::Disk { height: Height::Aging { .. }, .. }) {
                edges.extend(DISK_SLICES.iter().copied().filter(|&e| e > first && e < a));
            }
            edges.push(a);
            for w in edges.windows(2) {
                let (a0, a1) = (w[0], w[1]);
                let share = h.share_younger(a1) - h.share_younger(a0);
                if a1 <= a0 || share <= 1e-9 {
                    continue;
                }
                let mid = 0.5 * (a0 + a1);
                let metals = |t: f64, s: f64| galaxy.metallicity(c, [galaxy.reference_kpc * 1000.0, 0.0, 0.0], t, s);
                let table = Table::field(&galaxy.imf, (a0, a1), &|t| h.lookback_density(t) / share, &metals, seed.child_u64(b"table", c as u64).child_u64(b"from", a0.to_bits()), FIELD_KINDS);
                field.push(Population { component: c, age: mid, mass_msun: comp.mass_msun * share, shape: Prepared::new(&comp.shape, mid, galaxy.edge_kpc), table });
            }
            let share = h.share_younger(YOUNG_GYR);
            if y < YOUNG_GYR && share > 1e-9 {
                let (a0, a1) = (y, a.min(YOUNG_GYR));
                let mid = 0.5 * (a0 + a1);
                let metals = |t: f64, s: f64| galaxy.metallicity(c, [galaxy.reference_kpc * 1000.0, 0.0, 0.0], t, s);
                let table = Table::field(&galaxy.imf, (a0, a1), &|t| h.lookback_density(t) / share, &metals, seed.child_u64(b"table", c as u64).child(b"young field"), FIELD_KINDS);
                let small = comp.mass_msun * share * groups::SMALL_GROUP_SHARE;
                field.push(Population { component: c, age: mid, mass_msun: small, shape: Prepared::new(&comp.shape, mid, galaxy.edge_kpc), table });
                young.push(groups::YoungPopulation::new(galaxy, c, comp.mass_msun * share - small, seed));
            }
        }
        Self { galaxy: galaxy.clone(), seed, field, young }
    }

    #[must_use]
    pub fn galaxy(&self) -> &Galaxy {
        &self.galaxy
    }

    #[must_use]
    pub fn seed(&self) -> Seed {
        self.seed
    }

    /// The mass of populations (component, age, shape) between `lo` and `hi` as cells: a
    /// disk's plane is split once for all its slices, to half an arm's width where its young
    /// stars feel the arms.
    pub(crate) fn cells(&self, pops: &[(usize, f64, &Prepared)], lo: [f64; 3], hi: [f64; 3]) -> Vec<Vec<crate::galaxy::Cell>> {
        let g = &self.galaxy;
        let mut planes: Vec<(usize, Vec<crate::galaxy::PlaneLeaf>)> = Vec::new();
        pops.iter()
            .map(|&(c, age, shape)| {
                if !shape.is_disk() {
                    return shape.cells(lo, hi, &|_, _| 1.0, f64::MAX);
                }
                if !planes.iter().any(|p| p.0 == c) {
                    let detail = match &g.arms {
                        Some(a) if g.components[c].arm_contrast > 0.0 => 500.0 * a.width_kpc,
                        _ => f64::MAX,
                    };
                    let excess = |x: f64, y: f64| g.arms.as_ref().map_or(0.0, |a| a.excess(x, y));
                    planes.push((c, shape.plane([lo[0], lo[1]], [hi[0], hi[1]], &excess, detail)));
                }
                let plane = &planes.iter().find(|p| p.0 == c).expect("just made").1;
                shape.disk_cells(plane, g.arm_contrast(c, age), lo[2], hi[2])
            })
            .collect()
    }

    /// Mean systems per pc³ at `p` in bins `bins`: the field's, and the young groups' spread
    /// out.
    #[must_use]
    pub fn density(&self, p: [f64; 3], bins: std::ops::Range<usize>) -> f64 {
        let field: f64 = self
            .field
            .iter()
            .map(|pop| {
                let per: f64 = pop.table.per_msun[bins.clone()].iter().sum();
                pop.mass_msun * pop.shape.density(p) * self.galaxy.arm_factor(pop.component, pop.age, p[0], p[1]) * per
            })
            .sum();
        field + self.young.iter().map(|y| y.density(self, p, bins.clone())).sum::<f64>()
    }

    /// Stars within `radius_pc` of `p` with absolute V brighter than `faintest`: the field's
    /// and the young groups'.
    #[must_use]
    pub fn near(&self, p: [f64; 3], radius_pc: f64, faintest: f64) -> Vec<GalacticStar> {
        // The dark and just-dead bins once the query reaches past the lit ones.
        let top = if faintest >= bins::bright_edge(bins::DARK) { BINS - 1 } else { bins::bin_of(faintest) };
        let reach = |_: usize| radius_pc;
        let mut out = self.gather(p, 0..=top, &reach);
        out.retain(|s| s.absolute_v <= faintest && s.distance_pc(p) <= radius_pc);
        out
    }

    /// Stars that could be brighter than apparent magnitude `limit` seen from `p` (before
    /// dust), out to `max_pc`.
    #[must_use]
    pub fn visible(&self, p: [f64; 3], limit: f64, max_pc: f64) -> Vec<GalacticStar> {
        let reach = |b: usize| bins::reach_pc(bins::bright_edge(b), limit).min(max_pc);
        let mut out = self.gather(p, 0..=BINS - 1, &reach);
        out.retain(|s| {
            let d = s.distance_pc(p).max(1e-3);
            s.absolute_v + 5.0 * (d / 10.0).log10() <= limit && d <= max_pc
        });
        out
    }

    /// The field and group stars of `bins` within `reach(bin)` of `p`.
    pub(crate) fn gather(&self, p: [f64; 3], bins: std::ops::RangeInclusive<usize>, reach: &(dyn Fn(usize) -> f64 + Sync)) -> Vec<GalacticStar> {
        let mut out = Vec::new();
        // Field: by level, so a cube's cells serve all its bins.
        let mut by_level: Vec<(u32, Vec<usize>)> = Vec::new();
        for b in bins.clone() {
            let l = bins::level(b);
            match by_level.iter_mut().find(|e| e.0 == l) {
                Some(e) => e.1.push(b),
                None => by_level.push((l, vec![b])),
            }
        }
        // Cubes are independent, so they are generated in parallel and kept in order.
        for (level, bs) in &by_level {
            let r = bs.iter().map(|&b| reach(b)).fold(0.0, f64::max);
            let cubes: Vec<[i32; 3]> = bins::cubes_near(p, r, *level).collect();
            let stars: Vec<Vec<GalacticStar>> = cubes
                .par_iter()
                .map(|&cube| {
                    let d = bins::distance_to(p, cube, *level);
                    let wanted: Vec<usize> = bs.iter().copied().filter(|&b| d <= reach(b)).collect();
                    if wanted.is_empty() {
                        Vec::new()
                    } else {
                        self.field_cube(*level, cube, &wanted)
                    }
                })
                .collect();
            out.extend(stars.into_iter().flatten());
        }
        // Young groups.
        let r = bins.clone().map(reach).fold(0.0, f64::max);
        let cubes: Vec<[i32; 3]> = bins::cubes_near(p, r + groups::MAX_SPREAD_PC, groups::LEVEL).collect();
        let stars: Vec<Vec<GalacticStar>> = cubes
            .par_iter()
            .map(|&cube| {
                let mut found = Vec::new();
                for event in self.events(cube) {
                    let d = (0..3).map(|k| (event.centre[k] - p[k]).powi(2)).sum::<f64>().sqrt();
                    for b in bins.clone().filter(|&b| d <= reach(b) + event.extent_pc()) {
                        found.extend(self.members(cube, &event, b));
                    }
                }
                found
            })
            .collect();
        out.extend(stars.into_iter().flatten());
        out
    }

    /// The star `id`, if the galaxy generates it.
    #[must_use]
    pub fn star(&self, id: &StarId) -> Option<GalacticStar> {
        match *id {
            StarId::Field { bin, cube, .. } => self.field_cube(bins::level(usize::from(bin)), cube, &[usize::from(bin)]).into_iter().find(|s| s.id == *id),
            StarId::Group { bin, cube, event, .. } => {
                let e = self.events(cube).into_iter().find(|e| e.index == event)?;
                self.members(cube, &e, usize::from(bin)).into_iter().find(|s| s.id == *id)
            }
            StarId::Added { .. } | StarId::Cluster { .. } => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::galaxy::GalaxyPreset;

    fn mw() -> Sectors {
        Sectors::new(&Galaxy::milky_way(), Seed(5))
    }

    /// The Sun's neighbourhood: about 0.08 systems per pc³ (the 10 pc census holds 339
    /// systems), the thick disk a few per cent of the thin disk there (4 ± 2%), and the halo
    /// a tenth of a per cent: Bland-Hawthorn & Gerhard's halo mass and slope give that, under
    /// the 0.5% of star-count fits (Jurić et al. 2008), a tension in the literature.
    #[test]
    fn densities_at_the_sun() {
        let s = mw();
        let sun = [8200.0, 0.0, 25.0];
        let all = s.density(sun, 0..bins::LIT_BINS);
        assert!((0.06..0.10).contains(&all), "{all} systems per pc³");
        // A brown dwarf to every four stars, and the dark remnants about one in a thousand.
        let dark = s.density(sun, bins::DARK..BINS);
        assert!((0.2..0.32).contains(&(dark / all)), "dark bin {dark} against {all}");
        let by = |kind: ComponentKind| s.field.iter().filter(|p| s.galaxy.components[p.component].kind == kind).map(|p| p.mass_msun * p.shape.density(sun) * p.table.total_per_msun()).sum::<f64>();
        let (thin, thick, halo) = (by(ComponentKind::ThinDisk), by(ComponentKind::ThickDisk), by(ComponentKind::Halo));
        assert!((0.01..0.08).contains(&(thick / thin)), "thick {}", thick / thin);
        assert!((0.0003..0.01).contains(&(halo / (thin + thick))), "halo {}", halo / (thin + thick));
    }

    /// The same stars from anywhere, a deeper look adds only fainter ones, and identities
    /// regenerate their stars.
    #[test]
    fn stars_are_where_they_are() {
        let s = mw();
        let (a, b) = ([8200.0, 0.0, 25.0], [8203.0, 1.0, 24.0]);
        let from_a = s.near(a, 15.0, 12.0);
        let from_b = s.near(b, 15.0, 12.0);
        assert!(from_a.len() > 200, "{} systems within 15 pc", from_a.len());
        let both: Vec<&GalacticStar> = from_a.iter().filter(|x| x.distance_pc(b) <= 15.0).collect();
        assert!(both.iter().all(|x| from_b.contains(x)));
        let deeper = s.near(a, 15.0, 16.0);
        assert!(from_a.iter().all(|x| deeper.contains(x)) && deeper.len() > from_a.len());
        for x in from_a.iter().take(20) {
            assert_eq!(s.star(&x.id).as_ref(), Some(x));
        }
        // Stars sit nearer the plane than a few scale heights.
        let high = deeper.iter().filter(|x| x.position_pc[2].abs() > 1000.0).count();
        assert_eq!(high, 0);
    }

    /// An elliptical's stars are all old: no blue ones in a big volume.
    #[test]
    fn an_elliptical_has_no_young_stars() {
        let g = GalaxyPreset::Elliptical.galaxy();
        let s = Sectors::new(&g, Seed(1));
        assert!(s.young.iter().all(|y| y.mass_msun < 1e-6 * g.mass_msun()));
        let stars = s.visible([4000.0, 0.0, 0.0], 6.5, 3000.0);
        assert!(!stars.is_empty());
        assert!(stars.iter().all(|x| x.kind.age_gyr > 2.0), "{:?}", stars.iter().map(|x| x.kind.age_gyr).fold(f64::MAX, f64::min));
    }
}
