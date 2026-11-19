//! Elevation layer: metres above sea level.
//!
//! - Ocean floor depth follows crust age (GDH1, Stein & Stein 1992): 2600 m at the ridge
//!   crest, deepening as the plate cools to ~5650 m.
//! - Continents sit at a freeboard of a few hundred metres; provinces add their relief
//!   (plateaus behind collisions, narrow trenches, rift valleys with shoulders, seamount
//!   chains that sink as they age).
//! - Continental crust near the continent–ocean edge is submerged as shelf and slope, so
//!   land is less than the continental-crust area, as on Earth.
//! - Sea level is set so the oceans hold exactly the body's water (ocean global-equivalent
//!   layer × surface area); the stored field is then shifted so sea level is 0.
//! - Rivers then cut valleys (stream-power erosion, see `erosion`), while active mountain
//!   belts keep rising, so ranges stay high and drainage networks branch and merge.

use crate::noise::TerrainNoise;
use crate::tectonics::boundaries::{self, BoundaryParams};
use geoforge_core::grid::{distance_field, Cell, Grid};
use geoforge_core::product::{Generator, LayerId, Precondition, ProductData, ProductId, Products, RunContext};
use geoforge_core::shapes::Field;
use geoforge_core::surface::{CrustType, GeologicProvince as G};
use geoforge_core::{BodyContext, Error};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct ElevationParams {
    /// Continental freeboard: base elevation of stable continental crust, m.
    pub continental_base_m: f64,
    /// Scale of terrain noise (1 = Earth-like relief).
    pub roughness_scale: f64,
    /// Height of broad continental plateaus (dynamic topography, old uplifts), m.
    pub plateau_m: f64,
    /// Width of the continental slope, km.
    pub slope_width_km: f64,
    /// Extra sediment fill of closed land depressions, m, after erosion (rifts excepted).
    /// Erosion already breaches sills where rivers are big, so the default is 0.
    pub sediment_fill_m: f64,
    /// Fluvial erosion (stream power), run before the sediment fill.
    pub erosion: crate::erosion::ErosionParams,
    pub boundaries: BoundaryParams,
}

impl Default for ElevationParams {
    fn default() -> Self {
        Self { continental_base_m: 250.0, roughness_scale: 1.0, plateau_m: 600.0, slope_width_km: 130.0, sediment_fill_m: 0.0, erosion: crate::erosion::ErosionParams::default(), boundaries: BoundaryParams::default() }
    }
}

/// Ocean depth (positive, m) for crust age `t` Ma, GDH1.
#[must_use]
pub fn gdh1_depth_m(t: f64) -> f64 {
    let t = t.max(0.0);
    if t < 20.0 { 2600.0 + 365.0 * t.sqrt() } else { 5651.0 - 2473.0 * (-0.0278 * t).exp() }
}

fn smoothstep(e0: f64, e1: f64, x: f64) -> f64 {
    let t = ((x - e0) / (e1 - e0)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

pub struct Elevation;

impl Generator for Elevation {
    type Params = ElevationParams;
    const NAME: &'static str = "elevation";
    const VERSION: u32 = 1;
    const LAYER: LayerId = LayerId::Elevation;

    fn needs() -> &'static [ProductId] {
        &[ProductId::CrustType, ProductId::Provinces]
    }
    fn wants() -> &'static [ProductId] {
        &[ProductId::CrustAge, ProductId::Plates, ProductId::Lithology]
    }
    fn applies(body: &BodyContext) -> Result<(), Precondition> {
        crate::preconditions::plate_tectonics(body)
    }

    fn run(ctx: &RunContext<'_>, inputs: &Products, p: &ElevationParams) -> Result<Products, Error> {
        let mut out = Products::new();
        if inputs.elevation().is_some() {
            return Ok(out);
        }
        let crust = inputs.crust_type().ok_or(Error::MissingInput(ProductId::CrustType))?;
        let provinces = inputs.provinces().ok_or(Error::MissingInput(ProductId::Provinces))?;
        let grid = ctx.grid;
        let r = ctx.radius_km();
        let n = grid.len();
        let crust = &crust.values;
        let age: Vec<f32> = inputs.crust_age().map_or_else(|| vec![80.0; n], |f| f.values.clone());
        let noise = TerrainNoise::new(ctx.seed(ProductId::Elevation).value());
        let plateau_noise = TerrainNoise::new(ctx.seed(ProductId::Elevation).value() ^ 0x0071_A7EA);

        // Distance to the continent–ocean crust edge, both sides (km).
        let coast: Vec<Cell> = grid
            .cells()
            .filter(|&c| grid.neighbors(c).iter().any(|nb| crust[nb as usize] != crust[c as usize]))
            .collect();
        let edge_raw: Vec<f64> = distance_field(grid, coast, 1.0, |_| true).into_iter().map(|d| d * r).collect();
        // Coastlines are fractal: perturb distance-to-edge so margins form bays, capes and
        // offshore banks instead of following the smooth plate boundary.
        let coast_noise = TerrainNoise::new(ctx.seed(ProductId::Elevation).value() ^ 0xC0A57);
        let edge_km: Vec<f64> = grid
            .cells()
            .map(|c| {
                let pt = grid.center(c);
                // fbm is roughly ±0.3; these give ~±300 km capes and bays at ~800 km scale
                // plus ~±90 km detail at ~250 km scale.
                let n1 = coast_noise.fbm(pt, 8.0, 3);
                let n2 = coast_noise.fbm(pt, 25.0, 4);
                (edge_raw[c as usize] + 1000.0 * n1 + 300.0 * n2).max(0.0)
            })
            .collect();

        // Distance from each boundary-built province's source boundary (km).
        let segments = inputs.plates().map(|pl| boundaries::analyze(grid, pl, r, &p.boundaries)).unwrap_or_default();
        let by_owner = provinces.cells_by_owner();
        let mut need: std::collections::BTreeMap<usize, f64> = std::collections::BTreeMap::new();
        for key in by_owner.keys() {
            let Some(prov) = provinces.table.get(*key) else { continue };
            let Some(si) = prov.boundary else { continue };
            let reach_km = match prov.kind {
                G::CollisionOrogen => prov.width_km / 2.0 + 50.0,
                G::ContinentalRift => 450.0,
                G::MidOceanRidge => 200.0,
                G::OceanTrench => 80.0,
                _ => continue,
            };
            let e = need.entry(si).or_insert(0.0);
            *e = e.max(reach_km);
        }
        let mut seg_dist: HashMap<usize, HashMap<Cell, f64>> = HashMap::new();
        for (&si, &km) in &need {
            if let Some(seg) = segments.get(si) {
                seg_dist.insert(si, sparse_dist(grid, &seg.cells, km / r));
            }
        }
        let mut bdist: HashMap<Cell, f64> = HashMap::new();
        for (key, cells) in &by_owner {
            let Some(prov) = provinces.table.get(*key) else { continue };
            let Some(df) = prov.boundary.and_then(|si| seg_dist.get(&si)) else { continue };
            for c in cells {
                bdist.insert(*c, df.get(c).map_or(f64::INFINITY, |v| v * r));
            }
        }

        let spacing_km = grid.spacing_rad() * r;
        let rs = p.roughness_scale;
        // Depth inside each paleo-orogen (distance from the belt's edge), for a rounded profile.
        let mut interior_km: HashMap<Cell, f64> = HashMap::new();
        for (key, cells) in &by_owner {
            let Some(prov) = provinces.table.get(*key) else { continue };
            if prov.kind != G::PaleoOrogen {
                continue;
            }
            let set: std::collections::HashSet<Cell> = cells.iter().copied().collect();
            let rim: Vec<Cell> = cells.iter().copied().filter(|&c| grid.neighbors(c).iter().any(|n| !set.contains(&n))).collect();
            let d = sparse_dist(grid, &rim, prov.width_km / r);
            for &c in cells {
                interior_km.insert(c, d.get(&c).map_or(prov.width_km, |v| v * r));
            }
        }
        let base_cont = p.continental_base_m;

        // 1. Base + province relief (smooth part).
        let mut e: Vec<f64> = (0..n as Cell)
            .into_par_iter()
            .map(|c| {
                let i = c as usize;
                let ocean = crust[i] == CrustType::Oceanic;
                // Broad plateaus lift every continental province alike (Africa, Iran, Tibet).
                let base_cont = if ocean { base_cont } else { base_cont + p.plateau_m * smoothstep(0.0, 0.3, plateau_noise.fbm(grid.center(c), 1.1, 2)) };
                let mut z = if ocean {
                    // Sediment fill (isostatically compensated) thickens with age and toward
                    // continents (the continental rise).
                    let a = f64::from(age[i]);
                    let sediment = 380.0 * (a / 100.0).sqrt().min(1.4) + 1400.0 * (-edge_km[i] / 280.0).exp();
                    -gdh1_depth_m(a) + sediment
                } else {
                    let pt = grid.center(c);
                    base_cont + 650.0 * noise.fbm(pt, 1.6, 3) + 280.0 * noise.fbm(pt, 5.0, 3)
                };
                let Some(key) = provinces.owner_of(c) else { return z };
                let Some(prov) = provinces.table.get(key) else { return z };
                let d = bdist.get(&c).copied().unwrap_or(0.0);
                let half = (prov.width_km / 2.0).max(spacing_km);
                let pa = prov.anchor.map(|a| a.angle_to(grid.center(c)) * r);
                match prov.kind {
                    G::CollisionOrogen => {
                        let s = smoothstep(0.0, 1.0, 1.0 - d / half);
                        let h = 2800.0 + 2200.0 * prov.intensity;
                        z = base_cont + h * (1.6 * s).min(1.0) * (0.55 + 0.45 * s);
                    }
                    G::VolcanicArc => {
                        if ocean {
                            z = z.max(-900.0 + 1300.0 * prov.intensity);
                        } else {
                            z = base_cont + 1600.0 + 1400.0 * prov.intensity;
                        }
                    }
                    G::AccretionaryWedge => {
                        z = if ocean { z * 0.45 } else { base_cont + 300.0 };
                    }
                    G::ForearcBasin => {
                        z = if ocean { z.max(-3200.0) } else { base_cont - 250.0 };
                    }
                    G::BackarcBasin => {
                        z = if ocean { z.min(-2800.0) } else { base_cont - 300.0 };
                    }
                    G::OceanTrench => {}
                    G::MidOceanRidge => {
                        if prov.convergence_rate_cm_yr.abs() < 5.0 {
                            z -= 900.0 * (-(d / 18.0).powi(2)).exp();
                        }
                    }
                    G::OceanicPlateau => {
                        let rad = (prov.width_km / std::f64::consts::PI.sqrt()).max(spacing_km);
                        z += 2400.0 * smoothstep(0.0, 1.0, 1.0 - pa.unwrap_or(0.0) / rad).max(0.25);
                    }
                    G::OceanicHotspotTrack => {
                        // Seamounts sink as they age (they ride away from the hotspot).
                        let along = pa.unwrap_or(0.0);
                        let age_ma = along / 80.0;
                        z += 4200.0 * (-age_ma / 22.0).exp() + 1200.0;
                    }
                    G::ContinentalHotspotTrack => z = base_cont + 1100.0,
                    G::ContinentalFloodBasalt => z = base_cont + 600.0,
                    G::Craton => z = base_cont + 120.0,
                    G::Platform => z = base_cont - 60.0,
                    G::IntracratonicBasin => z = base_cont - 230.0,
                    G::ContinentalRift => {
                        z = base_cont + 900.0 * (1.0 - (-(d / 90.0).powi(2)).exp()) - 1300.0 * (-(d / 30.0).powi(2)).exp() + 250.0;
                    }
                    G::ExtendedCrust => {}
                    G::PaleoOrogen => {
                        let inside = interior_km.get(&c).copied().unwrap_or(0.0);
                        let s = smoothstep(0.0, (prov.width_km * 0.5).max(spacing_km), inside);
                        z = base_cont + (250.0 + 1400.0 * (-prov.age_ma / 450.0).exp()) * s;
                    }
                    G::AbyssalPlain | G::OceanicFractureZone | G::HotspotTrack => {}
                    // Rigid-lid kinds: drawn by `lid::LidElevation`.
                    G::CrateredHighlands
                    | G::CrateredIce
                    | G::VolcanicPlains
                    | G::ImpactBasin
                    | G::VolcanicRise
                    | G::ShieldVolcano
                    | G::Corona
                    | G::Tessera => {}
                }
                z
            })
            .collect();

        // 2. Soften province edges (two diffusion passes), then cut trenches.
        for _ in 0..2 {
            let prev = e.clone();
            e = (0..n as Cell)
                .into_par_iter()
                .map(|c| {
                    let nb = grid.neighbors(c);
                    let m: f64 = nb.iter().map(|x| prev[x as usize]).sum::<f64>() / nb.len() as f64;
                    0.5 * prev[c as usize] + 0.5 * m
                })
                .collect();
        }
        for (key, cells) in &by_owner {
            let Some(prov) = provinces.table.get(*key) else { continue };
            if prov.kind != G::OceanTrench {
                continue;
            }
            let depth = 7000.0 + 2500.0 * prov.intensity;
            for &c in cells {
                let d = bdist.get(&c).copied().unwrap_or(0.0);
                let z = -depth * (-(d / 35.0).powi(2)).exp().max(0.55);
                e[c as usize] = e[c as usize].min(z);
            }
        }

        // Continental margin zone: continental crust within slope + shelf width of the
        // continent–ocean edge. Shelves are built relative to sea level (below).
        let shelf_w = |c: Cell| -> f64 {
            match provinces.owner_of(c).and_then(|k| provinces.table.get(k)).map(|p| p.kind) {
                Some(G::ExtendedCrust) => 380.0,
                Some(G::VolcanicArc | G::AccretionaryWedge | G::CollisionOrogen | G::ForearcBasin) => 50.0,
                _ => 200.0,
            }
        };
        let sw = p.slope_width_km;
        let margin: Vec<Option<(f64, f64)>> = grid
            .cells()
            .map(|c| {
                let i = c as usize;
                let shelf = shelf_w(c);
                (crust[i] == CrustType::Continental && edge_km[i] < sw + shelf).then_some((edge_km[i], shelf))
            })
            .collect();

        // 3. Terrain detail scaled by province roughness (quiet on shelves).
        let freq = 6.0;
        e = (0..n as Cell)
            .into_par_iter()
            .map(|c| {
                let i = c as usize;
                let pt = grid.center(c);
                let prov = provinces.owner_of(c).and_then(|k| provinces.table.get(k));
                let rough = prov.map_or(0.2, |p| p.roughness);
                let kind = prov.map(|p| p.kind);
                let base = noise.fbm(pt, freq, 7);
                let land = crust[i] == CrustType::Continental;
                let amp = if margin[i].is_some() { 60.0 } else if land { 120.0 + 900.0 * rough } else { 120.0 + 350.0 * rough };
                let mut z = e[i] + rs * amp * base;
                if matches!(kind, Some(G::CollisionOrogen | G::VolcanicArc | G::PaleoOrogen | G::ContinentalRift)) {
                    let ridges = noise.ridged(pt, freq * 3.0, 6);
                    let mountain_amp = match kind {
                        Some(G::CollisionOrogen) => 2600.0 * prov.map_or(0.0, |p| 0.5 + 0.5 * p.intensity),
                        Some(G::VolcanicArc) => if crust[i] == CrustType::Oceanic { 1200.0 } else { 1800.0 },
                        Some(G::PaleoOrogen) => 700.0,
                        _ => 500.0,
                    };
                    z += rs * mountain_amp * (ridges - 0.35);
                }
                if matches!(kind, Some(G::OceanicHotspotTrack)) {
                    z += rs * 900.0 * (noise.ridged(pt, freq * 8.0, 3) - 0.5);
                }
                z
            })
            .collect();

        // 4. Sea level from ocean volume, with shelves planed to sea level (wave base):
        //    iterate shelves ↔ sea level to a fixed point.
        let gel_m = ctx.body.interior.surface_ocean_gel_km() * 1000.0;
        let areas: Vec<f64> = grid.cells().map(|c| grid.area_sr(c)).collect();
        let interior = e.clone();
        let detail: Vec<f64> = grid.cells().map(|c| 60.0 * rs * noise.fbm(grid.center(c), freq * 2.0, 4)).collect();
        let mut sea = if gel_m > 0.0 { sea_level_for_volume(&e, &areas, gel_m) } else { 0.0 };
        for _ in 0..6 {
            for (i, m) in margin.iter().enumerate() {
                let Some((d, shelf)) = *m else { continue };
                e[i] = if d < sw {
                    let t = smoothstep(0.0, sw, d);
                    sea - 3000.0 + (3000.0 - 130.0) * t + detail[i]
                } else {
                    let t = (d - sw) / shelf;
                    let shelf_z = sea - 130.0 + 140.0 * t + detail[i];
                    let w = smoothstep(0.75, 1.0, t);
                    shelf_z * (1.0 - w) + interior[i].max(shelf_z) * w
                };
            }
            if gel_m <= 0.0 {
                break;
            }
            let next = sea_level_for_volume(&e, &areas, gel_m);
            let done = (next - sea).abs() < 1.0;
            sea = next;
            if done {
                break;
            }
        }
        let mut field: Vec<f32> = e.iter().map(|z| (z - sea) as f32).collect();

        // 5. Rivers cut valleys (stream-power erosion; hard rock resists).
        if p.erosion.k > 0.0 {
            let erodibility = inputs.lithology().map_or_else(|| vec![0.6; n], crate::geology::lithology::erodibility);
            // Active belts keep rising while they erode.
            let uplift: Vec<f64> = (0..n)
                .map(|i| {
                    let kind = provinces.row_at(i as Cell).map(|pr| pr.kind);
                    let active = match kind {
                        Some(G::CollisionOrogen) => 1.0,
                        Some(G::VolcanicArc) => 0.8,
                        Some(G::ContinentalRift | G::ContinentalHotspotTrack) => 0.4,
                        _ => 0.0,
                    };
                    active * p.erosion.orogen_uplift * f64::from(field[i]).max(0.0)
                })
                .collect();
            let jitter = crate::hydrology::flat_jitter(grid, ctx.seed(ProductId::Elevation).value() ^ 0xE205);
            crate::erosion::erode(grid, &mut field, &erodibility, &uplift, &jitter, r, &p.erosion);
        }

        // 6. Long-term sediment fill: closed land basins fill up (or their sills
        //    are cut) by up to `sediment_fill_m`, except in active rifts.
        if p.sediment_fill_m > 0.0 {
            let d = crate::hydrology::drainage(grid, &field);
            let fill = p.sediment_fill_m;
            for i in 0..n {
                if field[i] < 0.0 {
                    continue;
                }
                let depth = d.filled[i] - f64::from(field[i]);
                if depth <= 0.0 {
                    continue;
                }
                let rift = provinces.owner_of(i as Cell).and_then(|k| provinces.table.get(k)).is_some_and(|pr| pr.kind == G::ContinentalRift);
                if !rift {
                    field[i] = (f64::from(field[i]) + depth.min(fill)) as f32;
                }
            }
        }
        out.insert(ProductData::Elevation(Field::from_values(grid.spec(), field)));
        Ok(out)
    }
}

fn sparse_dist(grid: &dyn Grid, sources: &[Cell], max_rad: f64) -> HashMap<Cell, f64> {
    use std::collections::BinaryHeap;
    #[derive(PartialEq)]
    struct It(f64, Cell);
    impl Eq for It {}
    impl Ord for It {
        fn cmp(&self, o: &Self) -> std::cmp::Ordering {
            o.0.partial_cmp(&self.0).unwrap_or(std::cmp::Ordering::Equal).then(o.1.cmp(&self.1))
        }
    }
    impl PartialOrd for It {
        fn partial_cmp(&self, o: &Self) -> Option<std::cmp::Ordering> {
            Some(self.cmp(o))
        }
    }
    let mut best: HashMap<Cell, f64> = HashMap::new();
    let mut heap = BinaryHeap::new();
    for &s in sources {
        if best.insert(s, 0.0).is_none() {
            heap.push(It(0.0, s));
        }
    }
    while let Some(It(d, c)) = heap.pop() {
        if best.get(&c).is_some_and(|&b| d > b) {
            continue;
        }
        let pc = grid.center(c);
        for nb in grid.neighbors(c).iter() {
            let nd = d + pc.angle_to(grid.center(nb));
            if nd <= max_rad && best.get(&nb).is_none_or(|&b| nd < b) {
                best.insert(nb, nd);
                heap.push(It(nd, nb));
            }
        }
    }
    best
}

/// The sea level (same datum as `e`) at which water of global-equivalent depth `gel_m`
/// fills exactly the cells below it.
#[must_use]
pub fn sea_level_for_volume(e: &[f64], area_sr: &[f64], gel_m: f64) -> f64 {
    let total: f64 = area_sr.iter().sum();
    let target = gel_m * total;
    let vol = |s: f64| e.iter().zip(area_sr).map(|(&z, &a)| if z < s { (s - z) * a } else { 0.0 }).sum::<f64>();
    let (mut lo, mut hi) = (e.iter().cloned().fold(f64::INFINITY, f64::min), e.iter().cloned().fold(f64::NEG_INFINITY, f64::max) + gel_m * 2.0);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if vol(mid) < target { lo = mid } else { hi = mid }
    }
    0.5 * (lo + hi)
}

/// Derived: is this cell under the sea? (Elevation is stored relative to sea level.)
#[must_use]
pub fn is_ocean(elevation_m: f32) -> bool {
    elevation_m < 0.0
}

/// Derived: slope magnitude (rise over run, dimensionless) at a cell.
#[must_use]
pub fn slope(grid: &dyn Grid, elev: &[f32], c: Cell, radius_km: f64) -> f64 {
    let pc = grid.center(c);
    grid.neighbors(c)
        .iter()
        .map(|n| (f64::from(elev[n as usize]) - f64::from(elev[c as usize])).abs() / (pc.angle_to(grid.center(n)) * radius_km * 1000.0))
        .fold(0.0, f64::max)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gdh1_matches_published_values() {
        assert!((gdh1_depth_m(0.0) - 2600.0).abs() < 1.0);
        assert!((gdh1_depth_m(10.0) - 3754.0).abs() < 5.0);
        assert!((gdh1_depth_m(100.0) - 5498.0).abs() < 10.0);
        assert!(gdh1_depth_m(180.0) < 5652.0);
    }

    #[test]
    fn sea_level_holds_volume() {
        let e = vec![-4000.0, -2000.0, 0.0, 1000.0];
        let a = vec![1.0; 4];
        let s = sea_level_for_volume(&e, &a, 1000.0);
        let v: f64 = e.iter().map(|z| (s - z).max(0.0)).sum();
        assert!((v - 4000.0).abs() < 1e-6, "{s} {v}");
    }
}
