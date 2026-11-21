//! Erosion layer: rivers cut valleys into the rough heightmap, and the result is the final
//! `elevation`.
//!
//! It runs after the climate, so rivers cut where the rain actually runs off: deep valleys
//! on wet windward ranges, little in the deserts. Without rain data (a rough heightmap with
//! no climate behind it) water is spread evenly, and a world with no liquid water keeps its
//! rough heightmap unchanged.
//!
//! Stream-power incision, `dh/dt = −K·A^m·S^n` with n = 1, solved with Braun & Willett's
//! (2013) implicit scheme: each step routes water down the steepest descent of the
//! depression-filled surface, accumulates drainage area, and lowers every cell toward its
//! receiver in one downstream-first pass. The scheme is stable for any step length, so a few
//! long steps are enough. There is no uplift: big rivers incise while ridges (small drainage
//! area) barely change, which carves dendritic valleys, breaches shallow basins and lets
//! drainage networks merge as they do on real land.
//!
//! Erodibility `K` scales with lithology (soft sediment wears fastest, cratonic basement and
//! granite slowest). Sea cells are base level and do not change.

use crate::hydrology::{drainage_with, runoff_mm, NONE};
use geoforge_core::grid::{Cell, Grid};
use geoforge_core::product::{Generator, LayerId, Precondition, ProductData, ProductId, Products, RunContext};
use geoforge_core::shapes::{Field, Partition};
use geoforge_core::surface::{GeologicProvince as G, Province};
use geoforge_core::{BodyContext, Error};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct ErosionParams {
    /// Erodibility of the softest rock, 1/yr (for m = 0.5). 0 disables erosion.
    pub k: f64,
    /// Drainage-area exponent m (n = 1).
    pub m: f64,
    /// Total erosion time, Myr.
    pub duration_myr: f64,
    /// Implicit steps (routing is recomputed each step).
    pub steps: usize,
    /// Hillslope diffusion per step (fraction toward the neighbour mean, land only).
    pub hillslope: f64,
    /// Rock uplift in active mountain belts during the erosion time, as a fraction of their
    /// tectonic height: ridges rise while rivers cut, so ranges stay high and rugged.
    pub orogen_uplift: f64,
    /// Runoff at which a cell counts as much as its area, mm/yr (about Earth's land mean).
    pub reference_runoff_mm: f64,
    /// How strongly a cell's runoff scales its share of drainage area: weight =
    /// (runoff / reference)^exponent. Below 1 because erosion depends only weakly on mean
    /// runoff: deserts erode in rare floods, and wet slopes are held by forest.
    pub runoff_exponent: f64,
    /// Fu's Budyko parameter for runoff from rain and potential evaporation (as hydrology).
    pub budyko_w: f64,
    /// Extra sediment fill of closed land depressions after erosion, m (rifts excepted).
    /// Erosion already breaches sills where rivers are big, so the default is 0.
    pub sediment_fill_m: f64,
}

impl Default for ErosionParams {
    fn default() -> Self {
        Self { k: 3e-6, m: 0.5, duration_myr: 1.0, steps: 8, hillslope: 0.02, orogen_uplift: 0.5, reference_runoff_mm: 300.0, runoff_exponent: 0.5, budyko_w: 3.0, sediment_fill_m: 0.0 }
    }
}

/// Erode `h` (m above sea level; cells below 0 are sea) in place. `erodibility` scales `k`
/// per cell (1 = softest); `uplift` is total rock uplift over the run per cell (m); `water`
/// weights each cell's share of drainage area (1 = the reference runoff; `None` = even
/// rain); `jitter` breaks ties on flats (see `hydrology::flat_jitter`).
#[allow(clippy::too_many_arguments)]
pub fn erode(grid: &dyn Grid, h: &mut [f32], erodibility: &[f64], uplift: &[f64], water: Option<&[f64]>, jitter: &[f32], radius_km: f64, p: &ErosionParams) {
    if p.k <= 0.0 || p.steps == 0 || p.duration_myr <= 0.0 {
        return;
    }
    let n = grid.len();
    let r_m = radius_km * 1000.0;
    let area_m2: Vec<f64> = (0..n as Cell).map(|c| grid.area_sr(c) * r_m * r_m * water.map_or(1.0, |w| w[c as usize])).collect();
    let dt = p.duration_myr * 1e6 / p.steps as f64;
    let mut receiver = vec![NONE; n];
    let mut dist = vec![1.0f64; n];
    for _ in 0..p.steps {
        for (c, u) in uplift.iter().enumerate() {
            if h[c] >= 0.0 && *u > 0.0 {
                h[c] += (*u / p.steps as f64) as f32;
            }
        }
        // Route on the depression-filled surface so water crosses basins.
        let d = drainage_with(grid, h, Some(jitter));
        for c in 0..n {
            receiver[c] = NONE;
            if h[c] < 0.0 {
                continue;
            }
            let pc = grid.center(c as Cell);
            let mut best = 0.0;
            for nb in grid.neighbors(c as Cell).iter() {
                let dd = pc.angle_to(grid.center(nb)) * r_m;
                let s = (d.filled[c] - d.filled[nb as usize]) / dd;
                if s > best {
                    best = s;
                    receiver[c] = nb;
                    dist[c] = dd;
                }
            }
            if receiver[c] == NONE {
                receiver[c] = d.parent[c];
                if receiver[c] != NONE {
                    dist[c] = pc.angle_to(grid.center(receiver[c])) * r_m;
                }
            }
        }
        // Drainage area, upstream first (reverse of the flood order, which is by filled level).
        let mut acc = area_m2.clone();
        for &c in d.order.iter().rev() {
            let r = receiver[c as usize];
            if r != NONE {
                acc[r as usize] += acc[c as usize];
            }
        }
        // Implicit incision, downstream first.
        for &c in &d.order {
            let ci = c as usize;
            let r = receiver[ci];
            if r == NONE || h[ci] < 0.0 {
                continue;
            }
            let hr = f64::from(h[r as usize]);
            let hi = f64::from(h[ci]);
            if hi <= hr {
                continue; // in a depression: left for sediment fill
            }
            let f = p.k * erodibility[ci] * acc[ci].powf(p.m) * dt / dist[ci];
            h[ci] = ((hi + f * hr) / (1.0 + f)).max(0.0) as f32;
        }
        // Hillslope diffusion (land only), so valley sides are slopes rather than steps.
        if p.hillslope > 0.0 {
            let prev = h.to_vec();
            for c in 0..n {
                if prev[c] <= 0.0 {
                    continue;
                }
                let nbs = grid.neighbors(c as Cell);
                let (mut s, mut k) = (0.0f64, 0.0f64);
                for nb in nbs.iter() {
                    s += f64::from(prev[nb as usize].max(0.0));
                    k += 1.0;
                }
                if k > 0.0 {
                    let v = f64::from(prev[c]) + p.hillslope * (s / k - f64::from(prev[c]));
                    h[c] = v.max(0.0) as f32;
                }
            }
        }
    }
}

pub struct Erosion;

impl Generator for Erosion {
    type Params = ErosionParams;
    const NAME: &'static str = "erosion";
    const VERSION: u32 = 1;
    const LAYER: LayerId = LayerId::Erosion;

    fn needs() -> &'static [ProductId] {
        &[ProductId::RoughElevation]
    }
    fn wants() -> &'static [ProductId] {
        &[ProductId::Precipitation, ProductId::Temperature, ProductId::Lithology, ProductId::Provinces]
    }
    fn applies(body: &BodyContext) -> Result<(), Precondition> {
        crate::preconditions::rocky_surface(body)
    }

    fn run(ctx: &RunContext<'_>, inputs: &Products, p: &ErosionParams) -> Result<Products, Error> {
        let mut out = Products::new();
        if inputs.elevation().is_some() {
            return Ok(out);
        }
        let rough = inputs.rough_elevation().ok_or(Error::MissingInput(ProductId::RoughElevation))?;
        let mut field = rough.values.clone();
        if ctx.body.interior.has_surface_liquid_water() && p.k > 0.0 {
            let grid = ctx.grid;
            let n = grid.len();
            let mut erodibility = inputs.lithology().map_or_else(|| vec![0.6; n], crate::geology::lithology::erodibility);
            let mut uplift = vec![0.0; n];
            if let Some(provinces) = inputs.provinces() {
                let interior = crate::elevation::plateau_interior_km(grid, provinces, ctx.radius_km());
                for (i, (k, u)) in erodibility.iter_mut().zip(uplift.iter_mut()).enumerate() {
                    let c = i as Cell;
                    let kind = provinces.row_at(c).map(|pr| pr.kind);
                    // Plateau interiors sit in their front range's rain shadow and drain into
                    // their own basins: no river crosses Tibet's interior to the sea.
                    if kind == Some(G::OrogenicPlateau) && interior.get(&c).is_some_and(|&d| d >= crate::elevation::PLATEAU_INTERIOR_KM) {
                        *k = 0.0;
                    }
                    // Active belts keep rising while they erode.
                    let active = match kind {
                        Some(G::CollisionOrogen) => 1.0,
                        Some(G::VolcanicArc) => 0.8,
                        Some(G::ContinentalRift | G::ContinentalHotspotTrack) => 0.4,
                        _ => 0.0,
                    };
                    *u = active * p.orogen_uplift * f64::from(field[i]).max(0.0);
                }
            }
            // Runoff from the climate's water balance, as a weight on drainage area.
            let water: Option<Vec<f64>> = inputs.precipitation().zip(inputs.temperature()).map(|(pr, t)| {
                pr.values
                    .iter()
                    .zip(&t.values)
                    .map(|(&rain, &temp)| (runoff_mm(f64::from(rain), crate::climate::pet_mm(f64::from(temp)), p.budyko_w) / p.reference_runoff_mm).powf(p.runoff_exponent).min(10.0))
                    .collect()
            });
            let jitter = crate::hydrology::flat_jitter(grid, ctx.seed(ProductId::Elevation).value() ^ 0xE205);
            erode(grid, &mut field, &erodibility, &uplift, water.as_deref(), &jitter, ctx.radius_km(), p);
            if p.sediment_fill_m > 0.0 {
                sediment_fill(grid, &mut field, inputs.provinces(), p.sediment_fill_m);
            }
        }
        out.insert(ProductData::Elevation(Field::from_values(rough.grid, field)));
        Ok(out)
    }
}

/// Closed land basins fill up (or their sills are cut) by up to `fill` m, except in active
/// rifts.
fn sediment_fill(grid: &dyn Grid, field: &mut [f32], provinces: Option<&Partition<Province>>, fill: f64) {
    let d = crate::hydrology::drainage(grid, field);
    for (i, h) in field.iter_mut().enumerate() {
        if *h < 0.0 {
            continue;
        }
        let depth = d.filled[i] - f64::from(*h);
        if depth <= 0.0 {
            continue;
        }
        let rift = provinces.and_then(|pv| pv.row_at(i as Cell)).is_some_and(|pr| pr.kind == G::ContinentalRift);
        if !rift {
            *h = (f64::from(*h) + depth.min(fill)) as f32;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use geoforge_core::grid::EquirectGrid;

    /// A tilted plateau draining to a sea strip: erosion cuts valleys but never below sea
    /// level, never raises land, and leaves the sea alone.
    #[test]
    fn incision_is_bounded() {
        let g = EquirectGrid::new(128, 64).unwrap();
        let n = g.len();
        let mut h: Vec<f32> = (0..n as u32)
            .map(|c| {
                let (x, y) = g.xy(c);
                if x < 8 { -500.0 } else { 3000.0 * (f64::from(x) / 128.0) as f32 + 50.0 * ((x * 7 + y * 13) % 5) as f32 }
            })
            .collect();
        let before = h.clone();
        let jitter = vec![0.5f32; n];
        erode(&g, &mut h, &vec![1.0; n], &vec![0.0; n], None, &jitter, 6371.0, &ErosionParams { hillslope: 0.0, ..ErosionParams::default() });
        let mut lowered = 0;
        for c in 0..n {
            if before[c] < 0.0 {
                assert_eq!(h[c], before[c]);
            } else {
                assert!(h[c] >= 0.0);
                assert!(h[c] <= before[c] + 1e-3, "raised at {c}");
                if h[c] < before[c] - 10.0 {
                    lowered += 1;
                }
            }
        }
        assert!(lowered > n / 20, "erosion should lower a good share of land: {lowered}");
    }
}
