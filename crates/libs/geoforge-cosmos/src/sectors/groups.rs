//! Young stars, born in groups that spread.
//!
//! Most stars form in groups that soon disperse: embedded clusters hold 70–90% of the stars
//! formed in giant molecular clouds, and fewer than 4–7% of them survive to the Pleiades'
//! age (Lada & Lada 2003). So stars younger than [`super::YOUNG_GYR`] come from birth events,
//! placed per 512 pc cube where young stars are (raised in the arms):
//!
//! - an event's mass formed is drawn from dN/dM ∝ M⁻² between 10 and 10⁵ M☉ (the cluster
//!   mass function), its age from the star formation history. Groups under 100 M☉ (a
//!   quarter of the mass under that law) scatter within a few Myr and are not told apart
//!   at a sky's scale: their stars are young field stars, raised in the arms like the rest
//!   ([`SMALL_GROUP_SHARE`]); the events are the groups of 100 M☉ and more;
//! - a twentieth stay bound, as open clusters with a Plummer radius of 1 pc × (M/100 M☉)^⅓;
//!   the rest spread as associations at 1.5 pc per Myr (a velocity spread of ~1.5 km/s) from
//!   2 pc, up to 150 pc, after which they are as good as field;
//! - members are drawn per bin on demand (each bin's draws seeded by the event and the bin),
//!   from the stars alive at the event's age, so a far observer sees an event's O and B
//!   stars without generating its dwarfs; none is heavier than half the event.

use super::bins;
use super::field::{cube_seed, place};
use super::tables::Table;
use super::{GalacticStar, Sectors, YOUNG_GYR};
use crate::draw::Draw;
use crate::galaxy::{Galaxy, Prepared};
use geoforge_core::galactic::{StarId, StarKind};
use geoforge_core::seed::Seed;

/// The level of group cubes: 8 pc × 2⁶ = 512 pc.
pub const LEVEL: u32 = 6;
/// The farthest a group's members are from its centre, pc (three of the largest spread).
pub const MAX_SPREAD_PC: f64 = 450.0;
/// Group masses formed, M☉, under dN/dM ∝ M⁻²: those made events, of all from 10 M☉.
const EVENT_MASS: (f64, f64) = (100.0, 1e5);
/// The share of young mass born in groups too small to be events: ln(100/10)/ln(10⁵/10).
pub const SMALL_GROUP_SHARE: f64 = 0.25;
const BOUND_SHARE: f64 = 0.05;
/// Ages the member tables are made at, Gyr.
const AGES: [f64; 10] = [0.002, 0.005, 0.01, 0.02, 0.035, 0.06, 0.1, 0.15, 0.21, 0.28];
const AGE_KINDS: usize = 5000;

/// The young stars of one component: their mass now, shape, and stars by age.
pub(crate) struct YoungPopulation {
    pub component: usize,
    pub mass_msun: f64,
    /// Ages of the young stars, Gyr.
    pub ages: (f64, f64),
    /// The age the shape and arms are taken at, Gyr.
    pub age: f64,
    pub shape: Prepared,
    pub tables: Vec<Table>,
}

impl YoungPopulation {
    pub(crate) fn new(galaxy: &Galaxy, component: usize, mass_msun: f64, seed: Seed) -> Self {
        let comp = &galaxy.components[component];
        let ages = (comp.history.youngest_gyr, comp.history.age_gyr.min(YOUNG_GYR));
        let age = 0.5 * (ages.0 + ages.1);
        let tables = AGES
            .iter()
            .enumerate()
            .map(|(i, &t)| {
                let feh = |_: f64, s: f64| galaxy.metallicity(component, [galaxy.reference_kpc * 1000.0, 0.0, 0.0], t, s);
                Table::at_age(&galaxy.imf, t, &feh, seed.child_u64(b"young", component as u64).child_u64(b"age", i as u64), AGE_KINDS)
            })
            .collect();
        Self { component, mass_msun, ages, age, shape: Prepared::new(&comp.shape, age, galaxy.edge_kpc), tables }
    }

    /// The table nearest `age_gyr`.
    fn table(&self, age_gyr: f64) -> &Table {
        let i = AGES.iter().enumerate().min_by(|a, b| (a.1.ln() - age_gyr.max(1e-4).ln()).abs().total_cmp(&(b.1.ln() - age_gyr.max(1e-4).ln()).abs())).map_or(0, |(i, _)| i);
        &self.tables[i]
    }

    /// Mean young systems per pc³ at `p` in bins `bins`.
    pub(crate) fn density(&self, s: &Sectors, p: [f64; 3], bins: std::ops::Range<usize>) -> f64 {
        let per: f64 = self.tables.iter().map(|t| t.per_msun[bins.clone()].iter().sum::<f64>()).sum::<f64>() / self.tables.len() as f64;
        self.mass_msun * self.shape.density(p) * s.galaxy().arm_factor(self.component, self.age, p[0], p[1]) * per
    }
}

/// A birth event: a cluster or an association.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Event {
    pub index: u32,
    /// Which young population it belongs to.
    pub population: usize,
    pub centre: [f64; 3],
    pub age_gyr: f64,
    /// Mass formed, M☉.
    pub mass_msun: f64,
    pub bound: bool,
    /// Plummer radius if bound, Gaussian spread if not, pc.
    pub radius_pc: f64,
    /// Its members' shared offset from their population's metallicity relation.
    pub scatter: f64,
}

impl Event {
    /// How far its members reach from its centre, pc.
    pub(crate) fn extent_pc(&self) -> f64 {
        if self.bound {
            10.0 * self.radius_pc
        } else {
            3.0 * self.radius_pc
        }
    }
}

/// The mean of dN/dM ∝ M⁻² between `lo` and `hi`.
fn mean_event_mass() -> f64 {
    let (lo, hi) = EVENT_MASS;
    (hi / lo).ln() / (1.0 / lo - 1.0 / hi)
}

impl Sectors {
    /// The young populations' cells between `lo` and `hi`.
    pub(crate) fn young_cells(&self, lo: [f64; 3], hi: [f64; 3]) -> Vec<Vec<crate::galaxy::Cell>> {
        let pops: Vec<(usize, f64, &Prepared)> = self.young.iter().map(|y| (y.component, y.age, &y.shape)).collect();
        self.cells(&pops, lo, hi)
    }

    /// The birth events of group cube `cube`.
    pub(crate) fn events(&self, cube: [i32; 3]) -> Vec<Event> {
        let (lo, hi) = bins::corners(cube, LEVEL);
        let mut out = Vec::new();
        for (yi, (yp, cells)) in self.young.iter().zip(self.young_cells(lo, hi)).enumerate() {
            let share: f64 = cells.iter().map(|c| c.mass).sum();
            let mean = yp.mass_msun * share / mean_event_mass();
            if mean <= 0.0 {
                continue;
            }
            let mut d = Draw::new(cube_seed(self.seed(), b"events", yi, cube));
            let h = &self.galaxy().components[yp.component].history;
            let peak = (0..=20).map(|i| h.lookback_density(yp.ages.0 + (yp.ages.1 - yp.ages.0) * f64::from(i) / 20.0)).fold(0.0, f64::max);
            if peak <= 0.0 {
                continue;
            }
            for _ in 0..d.poisson(mean) {
                let index = out.len() as u32;
                let centre = place(&mut d, &cells, share);
                let age_gyr = loop {
                    let t = yp.ages.0 + d.u() * (yp.ages.1 - yp.ages.0);
                    if d.u() * peak <= h.lookback_density(t) {
                        break t;
                    }
                };
                let (a, b) = EVENT_MASS;
                let mass_msun = 1.0 / (1.0 / a - d.u() * (1.0 / a - 1.0 / b));
                let bound = d.u() < BOUND_SHARE;
                let radius_pc = if bound { (mass_msun / 100.0).cbrt() } else { (2.0 + 1500.0 * age_gyr).min(MAX_SPREAD_PC / 3.0) };
                out.push(Event { index, population: yi, centre, age_gyr, mass_msun, bound, radius_pc, scatter: d.gauss() });
            }
        }
        out
    }

    /// The members of `event` (of group cube `cube`) in bin `bin`.
    pub(crate) fn members(&self, cube: [i32; 3], event: &Event, bin: usize) -> Vec<GalacticStar> {
        let yp = &self.young[event.population];
        let table = yp.table(event.age_gyr);
        let mean = event.mass_msun * table.per_msun[bin];
        if mean <= 0.0 {
            return Vec::new();
        }
        let seed = cube_seed(self.seed(), b"members", bin, cube).child_u64(b"event", u64::from(event.index));
        let mut d = Draw::new(seed);
        let kind_of = self.galaxy().components[yp.component].kind;
        let metallicity_dex = self.galaxy().metallicity(yp.component, event.centre, event.age_gyr, event.scatter);
        let mut out = Vec::new();
        for n in 0..d.poisson(mean) {
            let Some(k) = table.bins[bin].pick(d.u()) else { continue };
            let offset = if event.bound {
                let r = event.radius_pc / (d.u().max(1e-12).powf(-2.0 / 3.0) - 1.0).max(1e-12).sqrt();
                d.direction().map(|v| v * r.min(event.extent_pc()))
            } else {
                [0, 1, 2].map(|_| (event.radius_pc * d.gauss()).clamp(-event.extent_pc(), event.extent_pc()))
            };
            if k.initial_mass > 0.5 * event.mass_msun {
                continue;
            }
            let position = [0, 1, 2].map(|i| event.centre[i] + offset[i]);
            let kind = StarKind { initial_mass: k.initial_mass, age_gyr: event.age_gyr, metallicity_dex, companion_mass: k.companion };
            let id = StarId::Group { bin: bin as u8, cube, event: event.index, n: n as u32 };
            out.push(GalacticStar::new(id, position, kind, Some(kind_of)));
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::galaxy::Galaxy;

    /// Young stars clump: the brightest come in groups, their counts in 100 pc squares
    /// scattered far more than a Poisson spread's.
    #[test]
    fn young_stars_come_in_groups() {
        let s = Sectors::new(&Galaxy::milky_way(), Seed(9));
        let p = [8200.0, 0.0, 25.0];
        let young: Vec<GalacticStar> = s.near(p, 1000.0, -2.0).into_iter().filter(|x| matches!(x.id, StarId::Group { .. })).collect();
        assert!(young.len() > 30, "{} young bright stars within 1 kpc", young.len());
        assert!(young.iter().all(|x| x.kind.age_gyr < YOUNG_GYR));
        let mut counts = std::collections::HashMap::new();
        for x in &young {
            *counts.entry(((x.position_pc[0] / 100.0).floor() as i32, (x.position_pc[1] / 100.0).floor() as i32)).or_insert(0.0) += 1.0;
        }
        let cells = (std::f64::consts::PI * 100.0f64).round();
        let mean = young.len() as f64 / cells;
        let var = (counts.values().map(|c: &f64| (c - mean).powi(2)).sum::<f64>() + (cells - counts.len() as f64) * mean * mean) / cells;
        assert!(var / mean > 2.0, "dispersion {}", var / mean);
    }

    #[test]
    fn event_masses_follow_the_cluster_mass_function() {
        assert!((mean_event_mass() - 691.0).abs() < 2.0, "{}", mean_event_mass());
        let share = (EVENT_MASS.0 / 10.0f64).ln() / (EVENT_MASS.1 / 10.0f64).ln();
        assert!((share - SMALL_GROUP_SHARE).abs() < 1e-9);
    }
}
