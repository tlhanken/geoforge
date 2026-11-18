//! Crust type and age, and the inferred tectonic timeline.

use geoforge_core::grid::{Cell, Grid};
use geoforge_core::key::Key;
use geoforge_core::shapes::{Event, EventKind, Place, Table, Timeline};
use geoforge_core::surface::{BoundaryKind, BoundarySegment, CrustType};
use geoforge_core::{Seed, UnitVec3};
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::collections::{BTreeMap, BinaryHeap};

/// How crust type is generated when it is not authored.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(tag = "mode", rename_all = "snake_case")]
pub enum CrustMode {
    /// v1's rule: plates in the smallest `continental_percentile` by area are continental.
    PlateSize { continental_percentile: f64 },
    /// Continental crust covers about `target_fraction` of the surface, assembled from whole
    /// plates (smallest-first per v1, then by size) until the target is reached.
    TargetFraction { target_fraction: f64 },
}

impl Default for CrustMode {
    fn default() -> Self {
        CrustMode::TargetFraction { target_fraction: 0.40 }
    }
}

/// Dominant crust per plate.
#[must_use]
pub fn plate_kinds(areas: &BTreeMap<u32, f64>, mode: CrustMode, seed: Seed) -> BTreeMap<u32, CrustType> {
    let n = areas.len();
    match mode {
        CrustMode::PlateSize { continental_percentile } => areas
            .iter()
            .map(|(&id, &a)| {
                let smaller = areas.values().filter(|&&x| x < a).count();
                let pct = if n > 0 { smaller as f64 / n as f64 } else { 0.5 };
                (id, if pct > continental_percentile { CrustType::Oceanic } else { CrustType::Continental })
            })
            .collect(),
        CrustMode::TargetFraction { target_fraction } => {
            // Earth's largest plate (Pacific) is oceanic; continents ride on the rest. Rank the
            // other plates by a seeded draw weighted toward mid-sized plates, then take them
            // until the target is met.
            let total: f64 = areas.values().sum();
            let largest = areas.iter().max_by(|a, b| a.1.total_cmp(b.1).then(b.0.cmp(a.0))).map(|(&id, _)| id);
            let mut ranked: Vec<(f64, u32, f64)> = areas
                .iter()
                .filter(|(id, _)| Some(**id) != largest)
                .map(|(&id, &a)| {
                    let u = geoforge_core::seed::unit(seed.child_u64(b"continent", u64::from(id)));
                    let frac = a / total;
                    // Prefer plates of 3–25% of the surface; tiny and huge plates rank lower.
                    let size_pref = if frac < 0.03 { 0.5 } else if frac > 0.25 { 0.6 } else { 1.0 };
                    (u * size_pref, id, a)
                })
                .collect();
            ranked.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)));
            let mut out: BTreeMap<u32, CrustType> = areas.keys().map(|&id| (id, CrustType::Oceanic)).collect();
            let mut acc = 0.0;
            for (_, id, a) in ranked {
                if acc / total >= target_fraction {
                    break;
                }
                // Skip a plate that would overshoot badly if a smaller one could still fit.
                if (acc + a) / total > target_fraction + 0.12 && acc / total > target_fraction * 0.6 {
                    continue;
                }
                out.insert(id, CrustType::Continental);
                acc += a;
            }
            out
        }
    }
}

#[derive(PartialEq)]
struct Front(f64, Cell, f64, usize);
impl Eq for Front {}
impl Ord for Front {
    fn cmp(&self, o: &Self) -> Ordering {
        o.0.partial_cmp(&self.0).unwrap_or(Ordering::Equal).then_with(|| o.1.cmp(&self.1))
    }
}
impl PartialOrd for Front {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}

/// Oldest oceanic crust on Earth is ~180 Ma; older crust has been subducted.
pub const MAX_OCEANIC_AGE_MA: f64 = 180.0;

/// Crust age in Ma.
///
/// Oceanic crust ages away from divergent boundaries at the half-spreading rate (age =
/// distance / half rate) within its own plate. Oceanic crust with no ridge on its plate
/// ages from its plate's edges. Continental crust ages from ~400 Ma at plate edges to
/// ~3000 Ma in plate cores.
#[must_use]
pub fn crust_age(
    grid: &dyn Grid,
    owner: &[u32],
    crust: &[CrustType],
    segments: &[BoundarySegment],
    radius_km: f64,
) -> (Vec<f32>, Vec<Option<usize>>) {
    let n = grid.len();
    let mut age = vec![f64::INFINITY; n];
    let mut src = vec![None; n];
    let mut heap = BinaryHeap::new();
    for (si, s) in segments.iter().enumerate() {
        if s.kind != BoundaryKind::Divergent {
            continue;
        }
        // Half-spreading rate in km/Myr (1 cm/yr = 10 km/Myr), floored for ultra-slow ridges.
        let half = (s.rate_cm_yr * 10.0 / 2.0).max(5.0);
        for &c in &s.cells {
            if crust[c as usize] == CrustType::Oceanic && age[c as usize] > 0.0 {
                age[c as usize] = 0.0;
                src[c as usize] = Some(si);
                heap.push(Front(0.0, c, half, si));
            }
        }
    }
    while let Some(Front(a, c, half, si)) = heap.pop() {
        if a > age[c as usize] {
            continue;
        }
        let pc = grid.center(c);
        for nb in grid.neighbors(c).iter() {
            if owner[nb as usize] != owner[c as usize] || crust[nb as usize] != CrustType::Oceanic {
                continue;
            }
            let na = a + pc.angle_to(grid.center(nb)) * radius_km / half;
            if na < age[nb as usize] && na <= MAX_OCEANIC_AGE_MA {
                age[nb as usize] = na;
                src[nb as usize] = Some(si);
                heap.push(Front(na, nb, half, si));
            }
        }
    }

    // Distance from plate edges, used for ridge-less oceans and for continents.
    let edges: Vec<Cell> = grid
        .cells()
        .filter(|&c| grid.neighbors(c).iter().any(|nb| owner[nb as usize] != owner[c as usize]))
        .collect();
    let edge_dist = geoforge_core::grid::distance_field(grid, edges, std::f64::consts::PI, |_| true);
    let mut max_edge: BTreeMap<u32, f64> = BTreeMap::new();
    for c in grid.cells() {
        let d = edge_dist[c as usize];
        if d.is_finite() {
            let e = max_edge.entry(owner[c as usize]).or_insert(0.0);
            *e = e.max(d);
        }
    }

    let out = (0..n)
        .map(|i| {
            let o = owner[i];
            let rel = (edge_dist[i] / max_edge.get(&o).copied().unwrap_or(1.0).max(1e-9)).clamp(0.0, 1.0);
            match crust[i] {
                CrustType::Oceanic => {
                    if age[i].is_finite() {
                        age[i] as f32
                    } else {
                        // No ridge on this plate: old ocean, oldest in the middle.
                        (60.0 + 120.0 * rel) as f32
                    }
                }
                CrustType::Continental => (400.0 + 2600.0 * rel.powf(0.7)) as f32,
            }
        })
        .collect();
    (out, src)
}

fn mid(grid: &dyn Grid, cells: &[Cell]) -> UnitVec3 {
    UnitVec3::centroid(cells.iter().map(|&c| grid.center(c))).unwrap_or(UnitVec3::NORTH)
}

/// The tectonic timeline inferred from the present-day snapshot.
#[must_use]
pub fn infer_events(
    grid: &dyn Grid,
    owner: &[u32],
    crust: &[CrustType],
    age: &[f32],
    age_src: &[Option<usize>],
    segments: &[BoundarySegment],
) -> Timeline {
    let mut events = Table::new();
    let mut next = 0u64;
    let mut push = |e: Event, events: &mut Table<Event>| {
        next += 1;
        events.insert(Key(next), e);
    };

    // Ocean opening: the oldest crust each ridge produced.
    let mut ridge_max: BTreeMap<usize, f64> = BTreeMap::new();
    for (i, s) in age_src.iter().enumerate() {
        if let Some(si) = s {
            let e = ridge_max.entry(*si).or_insert(0.0);
            *e = e.max(f64::from(age[i]));
        }
    }

    for (si, s) in segments.iter().enumerate() {
        let at = Place::PlatePair { a: s.a, b: s.b, at: mid(grid, &s.cells) };
        let cont_frac = s.cells.iter().filter(|&&c| crust[c as usize] == CrustType::Continental).count() as f64 / s.cells.len() as f64;
        let rate_km_myr = (s.rate_cm_yr * 10.0).max(5.0);
        let (kind, t, note) = match s.kind {
            BoundaryKind::Divergent if cont_frac > 0.5 => {
                (EventKind::RiftOnset, (300.0 / rate_km_myr).clamp(2.0, 40.0), "continental rift opens")
            }
            BoundaryKind::Divergent => {
                (EventKind::OceanOpening, ridge_max.get(&si).copied().unwrap_or(10.0).max(1.0), "spreading ridge begins")
            }
            BoundaryKind::Convergent if cont_frac > 0.9 => {
                (EventKind::ContinentalCollision, (1500.0 / rate_km_myr).clamp(5.0, 80.0), "continents collide")
            }
            BoundaryKind::Convergent => {
                (EventKind::SubductionOnset, (1500.0 / rate_km_myr).clamp(5.0, 120.0), "subduction begins")
            }
            BoundaryKind::Transform => {
                (EventKind::TransformInitiation, (s.length_km * 0.3 / rate_km_myr).clamp(2.0, 60.0), "transform fault forms")
            }
        };
        push(
            Event { t_ma: (t * 10.0).round() / 10.0, kind, at, magnitude: s.rate_cm_yr, causes: vec![], note: note.into() },
            &mut events,
        );
    }

    // Craton formation: the oldest continental crust on each continental plate.
    let mut oldest: BTreeMap<u32, (f32, Cell)> = BTreeMap::new();
    for c in grid.cells() {
        if crust[c as usize] == CrustType::Continental {
            let e = oldest.entry(owner[c as usize]).or_insert((0.0, c));
            if age[c as usize] > e.0 {
                *e = (age[c as usize], c);
            }
        }
    }
    for (_, (a, c)) in oldest {
        push(
            Event {
                t_ma: f64::from(a).round(),
                kind: EventKind::CratonFormation,
                at: Place::Point { at: grid.center(c) },
                magnitude: 1.0,
                causes: vec![],
                note: "continental core stabilises".into(),
            },
            &mut events,
        );
    }
    Timeline { events }
}
