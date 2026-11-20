//! Plate boundaries, derived from the plates partition and the plates' Euler poles.
//!
//! v1 gave each plate pair one boundary type from the seed-to-seed direction. Here every
//! boundary cell is classified from the true relative velocity of the two plates at that
//! cell and the local boundary normal, using v1's 45° rule, so one plate pair can be
//! convergent along one stretch and transform along another (as Pacific–North America is).

use geoforge_core::grid::{components, Cell, Grid};
use geoforge_core::key::Key;
use geoforge_core::shapes::Partition;
use geoforge_core::surface::{BoundaryKind, BoundarySegment, Plate};
use geoforge_core::Vec3;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap, HashSet};

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct BoundaryParams {
    /// Relative motion within this angle of the boundary normal is convergent/divergent;
    /// beyond it, transform. Earth's catalogues call a boundary transform only within ~20° of
    /// pure strike-slip (PB2002: Sumatra and the Aleutians subduct at 50–70° obliquity).
    pub angle_threshold_deg: f64,
    /// Radius used to estimate the local boundary normal and to smooth the classification.
    pub normal_radius_km: f64,
    /// Stretches shorter than this merge into their neighbours, km.
    pub min_segment_km: f64,
}

impl Default for BoundaryParams {
    fn default() -> Self {
        Self { angle_threshold_deg: 70.0, normal_radius_km: 150.0, min_segment_km: 100.0 }
    }
}

/// Per-cell boundary measurements for one plate pair.
#[derive(Debug, Clone, Copy)]
struct Sample {
    kind: BoundaryKind,
    /// Normal relative speed, cm/yr, positive when converging.
    vn: f64,
    /// Relative speed, cm/yr.
    vrel: f64,
    /// Boundary normal in local (east, north) components.
    ne: f64,
    nn: f64,
}

fn bfs_rings(grid: &dyn Grid, start: Cell, rings: usize, allow: impl Fn(Cell) -> bool) -> Vec<Cell> {
    let mut seen = HashSet::new();
    seen.insert(start);
    let mut frontier = vec![start];
    let mut out = vec![start];
    for _ in 0..rings {
        let mut next = Vec::new();
        for &c in &frontier {
            for n in grid.neighbors(c).iter() {
                if allow(n) && seen.insert(n) {
                    next.push(n);
                    out.push(n);
                }
            }
        }
        if next.is_empty() {
            break;
        }
        frontier = next;
    }
    out
}

/// All boundary segments of a plates partition.
#[must_use]
pub fn analyze(grid: &dyn Grid, plates: &Partition<Plate>, radius_km: f64, p: &BoundaryParams) -> Vec<BoundarySegment> {
    let owner = &plates.owner.values;
    let spacing_km = grid.spacing_rad() * radius_km;
    let rings = ((p.normal_radius_km / spacing_km).ceil() as usize).max(2);
    // A stretch's cells lie on both sides of the boundary: two per cell of length.
    let min_cells = ((2.0 * p.min_segment_km / spacing_km).round() as usize).max(2);

    // 1. Boundary cells per unordered plate pair.
    let mut pairs: BTreeMap<(u32, u32), Vec<Cell>> = BTreeMap::new();
    for c in grid.cells() {
        let o = owner[c as usize];
        let mut seen_here: [u32; 8] = [0; 8];
        let mut k = 0;
        for n in grid.neighbors(c).iter() {
            let q = owner[n as usize];
            if q != o && q != 0 && o != 0 && !seen_here[..k].contains(&q) {
                seen_here[k] = q;
                k += 1;
                pairs.entry((o.min(q), o.max(q))).or_default().push(c);
            }
        }
    }

    // 2. Classify each boundary cell of each pair.
    let jobs: Vec<((u32, u32), Cell)> = pairs.iter().flat_map(|(&pr, cells)| cells.iter().map(move |&c| (pr, c))).collect();
    let samples: Vec<((u32, u32), Cell, Sample)> = jobs
        .par_iter()
        .map(|&((a, b), c)| {
            let own = owner[c as usize];
            let other = if own == a { b } else { a };
            let pc = grid.center(c);
            let mut s_own = Vec3::default();
            let mut s_oth = Vec3::default();
            let (mut n_own, mut n_oth) = (0usize, 0usize);
            for q in bfs_rings(grid, c, rings, |x| owner[x as usize] == own || owner[x as usize] == other) {
                let v = grid.center(q).v();
                if owner[q as usize] == own {
                    s_own = s_own + v;
                    n_own += 1;
                } else {
                    s_oth = s_oth + v;
                    n_oth += 1;
                }
            }
            let dir = s_oth * (1.0 / n_oth.max(1) as f64) - s_own * (1.0 / n_own.max(1) as f64);
            let normal = pc.tangent(dir).normalized().map_or(pc.east(), |u| u.v());
            let pa = plates.table.get(Key(u64::from(own)));
            let pb = plates.table.get(Key(u64::from(other)));
            let (va, vb) = match (pa, pb) {
                (Some(pa), Some(pb)) => (pa.velocity_km_per_myr(pc, radius_km), pb.velocity_km_per_myr(pc, radius_km)),
                _ => (Vec3::default(), Vec3::default()),
            };
            let rel = va - vb;
            let vn = rel.dot(normal) / 10.0;
            let vrel = rel.norm() / 10.0;
            let vt = (vrel * vrel - vn * vn).max(0.0).sqrt();
            let angle = vt.atan2(vn.abs()).to_degrees();
            let kind = if angle < p.angle_threshold_deg {
                if vn > 0.0 { BoundaryKind::Convergent } else { BoundaryKind::Divergent }
            } else {
                BoundaryKind::Transform
            };
            let ne = normal.dot(pc.east());
            let nn = normal.dot(pc.north());
            ((a, b), c, Sample { kind, vn, vrel, ne, nn })
        })
        .collect();

    let mut by_pair: BTreeMap<(u32, u32), HashMap<Cell, Sample>> = BTreeMap::new();
    for (pr, c, s) in samples {
        by_pair.entry(pr).or_default().insert(c, s);
    }

    // 3. Smooth, merge short stretches, and cut into segments.
    let mut segments = Vec::new();
    for ((a, b), mut cells) in by_pair {
        if cells.len() < min_cells {
            continue;
        }
        smooth_kinds(grid, &mut cells, rings);
        merge_short(grid, &mut cells, min_cells);
        let all: Vec<Cell> = {
            let mut v: Vec<Cell> = cells.keys().copied().collect();
            v.sort_unstable();
            v
        };
        for kind in [BoundaryKind::Convergent, BoundaryKind::Divergent, BoundaryKind::Transform] {
            let comps = components(grid, &all, |c| cells.get(&c).is_some_and(|s| s.kind == kind));
            for comp in comps {
                let n = comp.len() as f64;
                let rate = if kind == BoundaryKind::Transform {
                    comp.iter().map(|c| cells[c].vrel).sum::<f64>() / n
                } else {
                    comp.iter().map(|c| cells[c].vn.abs()).sum::<f64>() / n
                };
                let length_km = comp
                    .iter()
                    .map(|&c| {
                        let s = cells[&c];
                        let (ew, ns) = grid.cell_extent_rad(c);
                        radius_km / (s.ne.abs() / ew + s.nn.abs() / ns).max(1e-12)
                    })
                    .sum::<f64>()
                    / 2.0;
                segments.push(BoundarySegment {
                    a: Key(u64::from(a)),
                    b: Key(u64::from(b)),
                    kind,
                    cells: comp,
                    rate_cm_yr: rate,
                    length_km,
                });
            }
        }
    }
    segments.sort_by(|x, y| (x.a, x.b, x.kind, x.cells[0]).cmp(&(y.a, y.b, y.kind, y.cells[0])));
    segments
}

/// Majority vote of boundary kinds among nearby cells of the same pair.
fn smooth_kinds(grid: &dyn Grid, cells: &mut HashMap<Cell, Sample>, rings: usize) {
    let snapshot: HashMap<Cell, BoundaryKind> = cells.iter().map(|(&c, s)| (c, s.kind)).collect();
    let mut keys: Vec<Cell> = snapshot.keys().copied().collect();
    keys.sort_unstable();
    for c in keys {
        let mut votes = [0usize; 3];
        for q in bfs_rings(grid, c, rings, |x| snapshot.contains_key(&x)) {
            votes[snapshot[&q] as usize] += 1;
        }
        let own = snapshot[&c] as usize;
        let best = (0..3).max_by(|&i, &j| votes[i].cmp(&votes[j]).then((j == own).cmp(&(i == own)))).unwrap_or(own);
        if votes[best] > votes[own] {
            if let Some(s) = cells.get_mut(&c) {
                s.kind = [BoundaryKind::Convergent, BoundaryKind::Divergent, BoundaryKind::Transform][best];
            }
        }
    }
}

/// Stretches shorter than `min_cells` take the most common kind around them.
fn merge_short(grid: &dyn Grid, cells: &mut HashMap<Cell, Sample>, min_cells: usize) {
    for _ in 0..8 {
        let mut all: Vec<Cell> = cells.keys().copied().collect();
        all.sort_unstable();
        let mut changed = false;
        for kind in [BoundaryKind::Convergent, BoundaryKind::Divergent, BoundaryKind::Transform] {
            let comps = components(grid, &all, |c| cells.get(&c).is_some_and(|s| s.kind == kind));
            for comp in comps.into_iter().filter(|c| c.len() < min_cells) {
                let set: HashSet<Cell> = comp.iter().copied().collect();
                let mut votes = [0usize; 3];
                for &c in &comp {
                    for n in grid.neighbors(c).iter() {
                        if !set.contains(&n) {
                            if let Some(s) = cells.get(&n) {
                                votes[s.kind as usize] += 1;
                            }
                        }
                    }
                }
                if let Some(best) = (0..3).filter(|&i| votes[i] > 0).max_by(|&i, &j| votes[i].cmp(&votes[j]).then(j.cmp(&i))) {
                    let k = [BoundaryKind::Convergent, BoundaryKind::Divergent, BoundaryKind::Transform][best];
                    if k != kind {
                        for c in &comp {
                            if let Some(s) = cells.get_mut(c) {
                                s.kind = k;
                            }
                        }
                        changed = true;
                    }
                }
            }
        }
        if !changed {
            break;
        }
    }
}

/// Per-cell boundary kind for rendering: 0 none, 1 convergent, 2 divergent, 3 transform.
#[must_use]
pub fn boundary_kind_field(len: usize, segments: &[BoundarySegment]) -> Vec<u8> {
    let mut f = vec![0u8; len];
    for s in segments {
        for &c in &s.cells {
            f[c as usize] = s.kind as u8 + 1;
        }
    }
    f
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tectonics::motion::euler_pole_for;
    use geoforge_core::grid::{EquirectGrid, Grid};
    use geoforge_core::shapes::{Field, Table};
    use geoforge_core::UnitVec3;

    fn two_hemispheres(dir_north: f64, dir_south: f64) -> (EquirectGrid, Partition<Plate>) {
        let g = EquirectGrid::new(360, 180).unwrap();
        let owner: Vec<u32> = g.cells().map(|c| if g.center(c).z > 0.0 { 1 } else { 2 }).collect();
        let mut table = Table::new();
        for (k, lat, dir) in [(1u64, 45.0, dir_north), (2, -45.0, dir_south)] {
            let p = UnitVec3::from_lat_lon(lat, 0.0);
            let (axis, rate) = euler_pole_for(p, dir, 5.0, 6371.0);
            table.insert(Key(k), Plate { name: format!("P{k}"), euler_pole: axis, rate_deg_per_myr: rate, seed_point: p });
        }
        (g.clone(), Partition { owner: Field::from_values(g.spec(), owner), table })
    }

    #[test]
    fn equator_length_and_kinds() {
        // Each hemisphere rotates about an equatorial axis: near lon 0 the plates converge,
        // near lon 180 they separate (rigid rotation on a sphere cannot converge everywhere).
        let (g, plates) = two_hemispheres(270.0, 90.0);
        let segs = analyze(&g, &plates, 6371.0, &BoundaryParams::default());
        let total: f64 = segs.iter().map(|s| s.length_km).sum();
        assert!((total - 40_030.0).abs() / 40_030.0 < 0.03, "length {total}");
        let kind_at = |lat: f64, lon: f64| {
            let c = g.locate(UnitVec3::from_lat_lon(lat, lon));
            segs.iter().find(|s| s.cells.binary_search(&c).is_ok()).map(|s| s.kind)
        };
        assert_eq!(kind_at(0.5, 0.5), Some(BoundaryKind::Convergent));
        assert_eq!(kind_at(0.5, 179.5), Some(BoundaryKind::Divergent));
        // Convergence rate at lon 0 is the closing speed of two 5 cm/yr plates, ~10 cm/yr.
        let conv = segs.iter().find(|s| s.kind == BoundaryKind::Convergent).unwrap();
        assert!(conv.rate_cm_yr > 4.0 && conv.rate_cm_yr < 10.5, "{}", conv.rate_cm_yr);
    }
}
