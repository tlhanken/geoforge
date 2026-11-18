//! Plate partition on any [`Grid`]: Voronoi assignment, domain-warped boundaries, smoothing
//! and island removal (v1 Stages 1.1–1.3), plus area statistics.

use crate::noise::ValueNoise;
use geoforge_core::geom::UnitVec3;
use geoforge_core::grid::{components, Cell, Grid};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// Boundary roughening settings (v1 defaults at 1800 × 900: 80 px ≈ 16°).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct RefineParams {
    pub noise_scale: f64,
    /// Maximum displacement, degrees of arc.
    pub amplitude_deg: f64,
    pub octaves: usize,
    pub persistence: f64,
    pub smoothing_iterations: usize,
}

impl Default for RefineParams {
    fn default() -> Self {
        Self { noise_scale: 0.02, amplitude_deg: 16.0, octaves: 5, persistence: 0.5, smoothing_iterations: 1 }
    }
}

/// Assign every cell to the nearest seed (ids are 1-based seed indices).
#[must_use]
pub fn voronoi(grid: &dyn Grid, seeds: &[UnitVec3]) -> Vec<u32> {
    (0..grid.len() as Cell)
        .into_par_iter()
        .map(|c| {
            let p = grid.center(c);
            let mut best = 1u32;
            let mut min = f64::INFINITY;
            for (i, s) in seeds.iter().enumerate() {
                let d = p.angle_to(*s);
                if d < min {
                    min = d;
                    best = i as u32 + 1;
                }
            }
            best
        })
        .collect()
}

/// 3D spherical domain warping: each cell takes the plate found at a noise-displaced point.
#[must_use]
pub fn domain_warp(grid: &dyn Grid, owner: &[u32], seed: u64, p: &RefineParams) -> Vec<u32> {
    let noise = ValueNoise::new(seed, p.octaves, p.persistence);
    let amp = p.amplitude_deg.to_radians();
    let scale = p.noise_scale * 100.0;
    (0..grid.len() as Cell)
        .into_par_iter()
        .map(|c| {
            let o = grid.center(c);
            let (sx, sy, sz) = (o.x * scale, o.y * scale, o.z * scale);
            let dx = o.x + noise.fbm(sx, sy, sz, 0) * amp;
            let dy = o.y + noise.fbm(sx, sy, sz, 1) * amp;
            let dz = o.z + noise.fbm(sx, sy, sz, 2) * amp;
            let len = (dx * dx + dy * dy + dz * dz).sqrt();
            let w = UnitVec3 { x: dx / len, y: dy / len, z: dz / len };
            owner[grid.locate(w) as usize]
        })
        .collect()
}

/// Majority smoothing: a cell switches to a neighbour plate that holds at least three
/// quarters of its neighbours.
pub fn smooth(grid: &dyn Grid, owner: &mut Vec<u32>, iterations: usize) -> usize {
    let mut changed_total = 0;
    for _ in 0..iterations {
        let orig = owner.clone();
        let next: Vec<(u32, bool)> = (0..grid.len() as Cell)
            .into_par_iter()
            .map(|c| {
                let cur = orig[c as usize];
                let nb = grid.neighbors(c);
                let need = (nb.len() * 3).div_ceil(4);
                let mut counts: BTreeMap<u32, usize> = BTreeMap::new();
                for n in nb.iter() {
                    *counts.entry(orig[n as usize]).or_default() += 1;
                }
                match counts.into_iter().max_by(|a, b| a.1.cmp(&b.1).then(b.0.cmp(&a.0))) {
                    Some((p, k)) if p != cur && k >= need => (p, true),
                    _ => (cur, false),
                }
            })
            .collect();
        changed_total += next.iter().filter(|x| x.1).count();
        *owner = next.into_iter().map(|x| x.0).collect();
    }
    changed_total
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct IslandStats {
    pub plates_with_islands: usize,
    pub islands_removed: usize,
    pub cells_reassigned: usize,
}

/// Make every plate one connected piece: smaller pieces join the plate that surrounds them most.
pub fn remove_islands(grid: &dyn Grid, owner: &mut [u32]) -> IslandStats {
    let mut stats = IslandStats::default();
    let mut by_plate: BTreeMap<u32, Vec<Cell>> = BTreeMap::new();
    for (c, &o) in owner.iter().enumerate() {
        by_plate.entry(o).or_default().push(c as Cell);
    }
    for (&plate, cells) in &by_plate {
        let comps = components(grid, cells, |c| owner[c as usize] == plate);
        if comps.len() <= 1 {
            continue;
        }
        stats.plates_with_islands += 1;
        let main = comps
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.len().cmp(&b.1.len()).then(b.0.cmp(&a.0)))
            .map(|(i, _)| i)
            .unwrap_or(0);
        for (i, comp) in comps.iter().enumerate() {
            if i == main {
                continue;
            }
            let mut counts: BTreeMap<u32, usize> = BTreeMap::new();
            for &c in comp {
                for n in grid.neighbors(c).iter() {
                    let o = owner[n as usize];
                    if o != plate {
                        *counts.entry(o).or_default() += 1;
                    }
                }
            }
            let target = counts
                .into_iter()
                .max_by(|a, b| a.1.cmp(&b.1).then(b.0.cmp(&a.0)))
                .map_or(plate, |x| x.0);
            for &c in comp {
                owner[c as usize] = target;
            }
            stats.islands_removed += 1;
            stats.cells_reassigned += comp.len();
        }
    }
    stats
}

/// Area per owner in km² (cells with owner 0 are skipped).
#[must_use]
pub fn areas_km2(grid: &dyn Grid, owner: &[u32], radius_km: f64) -> BTreeMap<u32, f64> {
    let mut m: BTreeMap<u32, f64> = BTreeMap::new();
    let r2 = radius_km * radius_km;
    for (c, &o) in owner.iter().enumerate() {
        if o != 0 {
            *m.entry(o).or_default() += grid.area_sr(c as Cell) * r2;
        }
    }
    m
}
