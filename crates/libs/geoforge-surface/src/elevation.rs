//! Elevation layer: the rough heightmap, metres above sea level, before rivers cut it.
//! The erosion layer (`erosion`) cuts valleys into it once the climate says where the rain
//! falls, and writes the final `elevation`.
//!
//! - Ocean floor depth follows crust age (GDH1, Stein & Stein 1992): 2600 m at the ridge
//!   crest, deepening as the plate cools to ~5650 m.
//! - Continents sit at a freeboard of a few hundred metres; provinces add their relief
//!   (plateaus behind collisions, narrow trenches, rift valleys with shoulders, seamount
//!   chains that sink as they age).
//! - Continental collisions are asymmetric: a foreland basin on the plate going under, a
//!   steep front range across the boundary, and a plateau behind it whose interior barely
//!   erodes (see `geology::provinces::Collision`). Their heights scale with gravity
//!   (`relief`).
//! - Ocean–continent margins rise from the trench to a coastal range, a volcanic cordillera
//!   (higher where convergence is fast) and, behind an old slab, an Andean plateau; their
//!   coasts run straight along the trench, with a narrow shelf.
//! - Continental crust near the continent–ocean edge is submerged as shelf and slope, so
//!   land is less than the continental-crust area, as on Earth.
//! - Sea level is set so the oceans hold exactly the body's water (ocean global-equivalent
//!   layer × surface area); the stored field is then shifted so sea level is 0.

use crate::geology::provinces as pc;
use crate::noise::TerrainNoise;
use crate::tectonics::boundaries::{self, BoundaryParams};
use geoforge_core::grid::{distance_field, Cell, Grid};
use geoforge_core::product::{Generator, LayerId, Precondition, ProductData, ProductId, Products, RunContext};
use geoforge_core::shapes::{Field, Partition};
use geoforge_core::surface::{CrustType, GeologicProvince, GeologicProvince as G, Province};
use geoforge_core::{BodyContext, Error, UnitVec3};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};

/// Island arcs: the submarine ridge, the tallest volcanoes on it (m), and their spacing (km).
const ARC_RIDGE_M: f64 = -1800.0;
/// Block uplifts over a flat slab: the regional rise, the ranges' height above it, the
/// ranges' spacing and the taper at the province's edge.
const BASEMENT_RISE_M: f64 = 1200.0;
const BASEMENT_RANGE_M: f64 = 2200.0;
const BASEMENT_BLOCK_KM: f64 = 180.0;
const BASEMENT_EDGE_KM: f64 = 120.0;
/// Peaks of an island arc's strongest stretches: the Aleutians' reach 2.9 km, the Lesser
/// Antilles' 1.5 km, the Marianas' 1 km.
const ARC_PEAK_M: f64 = 2400.0;
/// A plume under a continent lifts a dome this high and this wide (1/e radius): Ethiopia's
/// and East Africa's are 1.5–2 km over ~1,000 km, the Hoggar's and Tibesti's ~1 km.
const PLUME_DOME_M: f64 = 1500.0;
const PLUME_DOME_RADIUS_KM: f64 = 650.0;
const ARC_VOLCANO_SPACING_KM: f64 = 70.0;
/// Length of the stretches an arc's height varies over.
const ARC_SEGMENT_KM: f64 = 900.0;
/// Flood basalts younger than this still sit on their plume's dome (Ethiopia's are 30 Ma).
const PLUME_HEAD_MAX_AGE_MA: f64 = 45.0;
/// Dynamic topography: continents are lowest within the first distance behind an active
/// margin's arc or orogen and highest beyond the second, km; `DYNAMIC_LOW_SHARE` of the swell is below
/// the base; and smooth noise adds or takes this much, m.
const DYNAMIC_NEAR_KM: f64 = 1000.0;
const DYNAMIC_FAR_KM: f64 = 3500.0;
const DYNAMIC_LOW_SHARE: f64 = 0.0;
const DYNAMIC_NOISE_M: f64 = 300.0;
/// The broad swell under orogens and arcs, m.
const OROGEN_SWELL_M: f64 = 600.0;
/// Plates faster than the first speed lose swell, up to `DYNAMIC_FAST_CUT` of it at the
/// second, cm/yr (Australia moves 7 cm/yr, Africa 2).
const DYNAMIC_SLOW_CM_YR: f64 = 3.0;
const DYNAMIC_FAST_CM_YR: f64 = 7.0;
const DYNAMIC_FAST_CUT: f64 = 0.7;
/// An orogenic plateau rises to full height over this far in from its outer edges, km.
const PLATEAU_EDGE_KM: f64 = 120.0;
/// Within this distance of a convergent boundary coastlines are not roughened: an active
/// margin's coast runs straight along its trench.
const STRAIGHT_COAST_KM: f64 = 800.0;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct ElevationParams {
    /// Continental freeboard: base elevation of stable continental crust, m.
    pub continental_base_m: f64,
    /// Scale of terrain noise (1 = Earth-like relief).
    pub roughness_scale: f64,
    /// Dynamic topography: how much higher a continent stands far from every convergent
    /// boundary than right behind one, m (mantle flow drags continents down above sinking
    /// slabs and lifts them far from any: Africa and Brazil's highlands stand high, the
    /// Amazon and Mesopotamia low).
    pub plateau_m: f64,
    /// Width of the continental slope, km.
    pub slope_width_km: f64,
    pub boundaries: BoundaryParams,
}

impl Default for ElevationParams {
    fn default() -> Self {
        Self { continental_base_m: 250.0, roughness_scale: 1.0, plateau_m: 1800.0, slope_width_km: 130.0, boundaries: BoundaryParams::default() }
    }
}

/// Ocean depth (positive, m) for crust age `t` Ma, GDH1.
#[must_use]
pub fn gdh1_depth_m(t: f64) -> f64 {
    let t = t.max(0.0);
    if t < 20.0 {
        2600.0 + 365.0 * t.sqrt()
    } else {
        5651.0 - 2473.0 * (-0.0278 * t).exp()
    }
}

pub(crate) fn smoothstep(e0: f64, e1: f64, x: f64) -> f64 {
    let t = ((x - e0) / (e1 - e0)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

pub struct Elevation;

impl Generator for Elevation {
    type Params = ElevationParams;
    const NAME: &'static str = "elevation";
    const VERSION: u32 = 9;
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
        if inputs.rough_elevation().is_some() {
            return Ok(out);
        }
        // Randomness is drawn from the `elevation` seed, as before the heightmap was split
        // into rough and eroded layers, so the split changed no existing world.
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
        let coast: Vec<Cell> = grid.cells().filter(|&c| grid.neighbors(c).iter().any(|nb| crust[nb as usize] != crust[c as usize])).collect();
        let edge_raw: Vec<f64> = distance_field(grid, coast, 1.0, |_| true).into_iter().map(|d| d * r).collect();
        // Coastlines are fractal: perturb distance-to-edge so margins form bays, capes and
        // offshore banks instead of following the smooth plate boundary. Not at active
        // margins: there the coast runs straight along the trench (Chile, Peru, Sumatra).
        let segments = inputs.plates().map(|pl| boundaries::analyze(grid, pl, r, &p.boundaries)).unwrap_or_default();
        let convergent: Vec<Cell> = segments.iter().filter(|s| s.kind == geoforge_core::surface::BoundaryKind::Convergent).flat_map(|s| s.cells.iter().copied()).collect();
        let near_trench = sparse_dist(grid, &convergent, STRAIGHT_COAST_KM / r);
        // Distance behind the nearest active margin's upper plate (its arc or orogen), where
        // the sinking slab drags the surface down, and how fast each cell's plate moves.
        let upper: Vec<Cell> = grid.cells().filter(|&c| provinces.row_at(c).is_some_and(|p| matches!(p.kind, G::VolcanicArc | G::CollisionOrogen | G::OrogenicPlateau))).collect();
        let behind_margin_km: Vec<f64> =
            if upper.is_empty() { vec![DYNAMIC_FAR_KM; n] } else { distance_field(grid, upper, 1.05 * DYNAMIC_FAR_KM / r, |_| true).into_iter().map(|d| (d * r).min(DYNAMIC_FAR_KM)).collect() };
        let plate_speed_cm_yr: Vec<f64> = match inputs.plates() {
            Some(pl) => {
                grid.cells().map(|c| pl.table.get(geoforge_core::key::Key(u64::from(pl.owner.values[c as usize]))).map_or(0.0, |p| p.velocity_km_per_myr(grid.center(c), r).norm() / 10.0)).collect()
            }
            None => vec![0.0; n],
        };
        let coast_noise = TerrainNoise::new(ctx.seed(ProductId::Elevation).value() ^ 0xC0A57);
        let edge_km: Vec<f64> = grid
            .cells()
            .map(|c| {
                let pt = grid.center(c);
                // fbm is roughly ±0.3; these give ~±300 km capes and bays at ~800 km scale
                // plus ~±90 km detail at ~250 km scale.
                let n1 = coast_noise.fbm(pt, 8.0, 3);
                let n2 = coast_noise.fbm(pt, 25.0, 4);
                let open = near_trench.get(&c).map_or(1.0, |d| smoothstep(0.5 * STRAIGHT_COAST_KM, STRAIGHT_COAST_KM, d * r));
                (edge_raw[c as usize] + open * 1000.0 * n1 + (0.3 + 0.7 * open) * 300.0 * n2).max(0.0)
            })
            .collect();

        // Distance (km) from the nearest source boundary of the province's own kind. Measured
        // per kind rather than per province so that where rifts or belts meet at a junction
        // their valleys and crests merge, instead of one branch's flank crossing another's axis.
        let by_owner = provinces.cells_by_owner();
        let mut axes: BTreeMap<G, (Vec<Cell>, f64)> = BTreeMap::new();
        for key in by_owner.keys() {
            let Some(prov) = provinces.table.get(*key) else { continue };
            let Some(seg) = prov.boundary.and_then(|si| segments.get(si)) else { continue };
            let reach_km = match prov.kind {
                G::CollisionOrogen | G::ForelandBasin => prov.width_km * 1.2 + 50.0,
                G::ContinentalRift => 450.0,
                G::MidOceanRidge => 200.0,
                G::OceanTrench => 80.0,
                _ => continue,
            };
            let axis = axes.entry(prov.kind).or_default();
            axis.0.extend_from_slice(&seg.cells);
            axis.1 = axis.1.max(reach_km);
        }
        let axis_km: BTreeMap<G, HashMap<Cell, f64>> = axes.into_iter().map(|(kind, (cells, reach_km))| (kind, sparse_dist(grid, &cells, reach_km / r))).collect();
        let bdist = |c: Cell| -> Option<f64> {
            let kind = provinces.row_at(c)?.kind;
            Some(axis_km.get(&kind)?.get(&c).map_or(f64::INFINITY, |v| v * r))
        };

        let spacing_km = grid.spacing_rad() * r;
        let rs = p.roughness_scale;
        // Depth inside each paleo-orogen (distance from the belt's edge), for a rounded profile.
        let mut interior_km: HashMap<Cell, f64> = HashMap::new();
        for (key, cells) in &by_owner {
            let Some(prov) = provinces.table.get(*key) else { continue };
            if !matches!(prov.kind, G::PaleoOrogen | G::BasementUplift) {
                continue;
            }
            let set: std::collections::HashSet<Cell> = cells.iter().copied().collect();
            let rim: Vec<Cell> = cells.iter().copied().filter(|&c| grid.neighbors(c).iter().any(|n| !set.contains(&n))).collect();
            let d = sparse_dist(grid, &rim, prov.width_km / r);
            for &c in cells {
                interior_km.insert(c, d.get(&c).map_or(prov.width_km, |v| v * r));
            }
        }
        interior_km.extend(plateau_interior_km(grid, provinces, r));
        let base_cont = p.continental_base_m;
        let relief = crate::relief::relief_scale(ctx.body);
        let collision = |prov: &geoforge_core::surface::Province| {
            let c = pc::Collision::at_rate(prov.convergence_rate_cm_yr);
            pc::Collision { range_m: c.range_m * relief, plateau_m: c.plateau_m * relief, foreland_m: c.foreland_m * relief, ..c }
        };

        // Arcs rise and fall along their length, in stretches of several hundred km: the
        // central Andes stand 2 km above the northern; volcano-free gaps open where the
        // factor is lowest.
        let along_arc = |pt: UnitVec3| (0.85 + 1.4 * plateau_noise.fbm(pt, r / ARC_SEGMENT_KM, 2)).clamp(0.35, 1.3);
        // Plume domes: under every active continental hotspot and young flood basalt.
        let plumes: Vec<UnitVec3> = provinces
            .table
            .rows
            .values()
            .filter(|p| p.kind == G::ContinentalHotspotTrack || (p.kind == G::ContinentalFloodBasalt && p.age_ma < PLUME_HEAD_MAX_AGE_MA))
            .filter_map(|p| p.anchor)
            .collect();
        let dome = |pt: UnitVec3| plumes.iter().map(|a| (-(a.angle_to(pt) * r / PLUME_DOME_RADIUS_KM).powi(2)).exp()).fold(0.0, f64::max);

        // 1. Base + province relief (smooth part).
        let mut e: Vec<f64> = (0..n as Cell)
            .into_par_iter()
            .map(|c| {
                let i = c as usize;
                let ocean = crust[i] == CrustType::Oceanic;
                // Broad plateaus lift every continental province alike (Africa, Iran, Tibet).
                let kind = provinces.row_at(c).map(|p| p.kind);
                let base_cont = if ocean {
                    base_cont
                } else if matches!(kind, Some(G::CollisionOrogen | G::OrogenicPlateau | G::VolcanicArc)) {
                    // Orogens and arcs stand on their own crustal roots; only a broad, smooth
                    // swell lifts some of them further (calibrated with their heights).
                    base_cont + OROGEN_SWELL_M * smoothstep(0.0, 0.3, plateau_noise.fbm(grid.center(c), 1.1, 2))
                } else {
                    // Far from any slab and slow over the mantle (Africa) stands high; behind an
                    // active margin (the Amazon) or racing over old slabs (Australia, India), low.
                    let far = smoothstep(DYNAMIC_NEAR_KM, DYNAMIC_FAR_KM, behind_margin_km[i]);
                    let slow = 1.0 - DYNAMIC_FAST_CUT * smoothstep(DYNAMIC_SLOW_CM_YR, DYNAMIC_FAST_CM_YR, plate_speed_cm_yr[i]);
                    // Upwellings are patchy: the swell covers part of a far, slow interior
                    // (southern and eastern Africa, not the Sahara).
                    let up = (0.65 + 0.45 * plateau_noise.fbm(grid.center(c), 1.1, 3) / 0.3).clamp(0.3, 1.0);
                    let swell = (plateau_noise.fbm(grid.center(c), 2.3, 3) / 0.3).clamp(-0.2, 1.5);
                    base_cont + p.plateau_m * (far * slow * up - DYNAMIC_LOW_SHARE) + DYNAMIC_NOISE_M * swell
                };
                let mut z = if ocean {
                    // Sediment fill (isostatically compensated) thickens with age and toward
                    // continents (the continental rise).
                    let a = f64::from(age[i]);
                    let sediment = 380.0 * (a / 100.0).sqrt().min(1.4) + 1400.0 * (-edge_km[i] / 280.0).exp();
                    -gdh1_depth_m(a) + sediment
                } else {
                    let pt = grid.center(c);
                    base_cont + 650.0 * noise.fbm(pt, 1.6, 3) + 280.0 * noise.fbm(pt, 5.0, 3) + relief * PLUME_DOME_M * dome(pt)
                };
                let Some(key) = provinces.owner_of(c) else { return z };
                let Some(prov) = provinces.table.get(key) else { return z };
                let d = bdist(c).unwrap_or(0.0);
                let pa = prov.anchor.map(|a| a.angle_to(grid.center(c)) * r);
                match prov.kind {
                    // Steep front from the foreland up to the crest, then down the back to
                    // the plateau (or to the plain behind, where there is none).
                    G::CollisionOrogen => {
                        let shape = collision(prov);
                        let crest = 0.25 * shape.range_km;
                        z = if d < crest {
                            base_cont + shape.range_m * smoothstep(0.0, crest, d)
                        } else {
                            base_cont + shape.range_m + (shape.plateau_m - shape.range_m) * smoothstep(crest, shape.range_km, d)
                        };
                    }
                    G::OrogenicPlateau => {
                        // Ragged margins: the edge ramp starts early or late along the edge.
                        let inside = interior_km.get(&c).copied().unwrap_or(2.0 * PLATEAU_EDGE_KM) + PLATEAU_EDGE_KM * noise.fbm(grid.center(c), 12.0, 3);
                        z = base_cont + collision(prov).plateau_m * prov.elevation_intensity * smoothstep(0.0, PLATEAU_EDGE_KM, inside);
                    }
                    // Bent down under the range's load and filled with its debris.
                    G::ForelandBasin => {
                        let shape = collision(prov);
                        z = base_cont - shape.foreland_m * (1.0 - smoothstep(0.0, shape.foreland_km, d));
                    }
                    G::VolcanicArc => {
                        let seg = along_arc(grid.center(c));
                        if ocean {
                            // The submarine ridge an island arc's volcanoes stand on (step 3);
                            // its strongest stretches have thickened into larger islands.
                            z = z.max(ARC_RIDGE_M + 600.0 * prov.intensity + 1600.0 * smoothstep(0.9, 1.3, seg));
                        } else {
                            // A continental cordillera: ~1.5 km on a slow margin (the
                            // Cascades), ~4 km where convergence is fast (the Andes).
                            // Cordilleras vary less along strike than island arcs (±20%).
                            let fast = smoothstep(2.0, 8.0, prov.convergence_rate_cm_yr);
                            z = base_cont + relief * (1000.0 + 3000.0 * fast) * (1.0 + 0.6 * (seg - 0.85)).clamp(0.75, 1.15);
                        }
                    }
                    G::AccretionaryWedge => {
                        z = if ocean { z * 0.45 } else { base_cont + 300.0 };
                    }
                    G::ForearcBasin => {
                        // Onshore at an Andean margin: the coastal range and central valley.
                        z = if ocean { z.max(-3200.0) } else { base_cont + 300.0 * relief };
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
                        // Higher and lower stretches along the belt (the Appalachians reach 2 km in
                        // the Smokies and 600 m in Pennsylvania's ridges).
                        let along = (0.6 + 1.6 * plateau_noise.fbm(grid.center(c), 9.0, 2)).clamp(0.25, 1.4);
                        z = base_cont + (250.0 + 1400.0 * (-prov.age_ma / 450.0).exp()) * s * along;
                    }
                    G::BasementUplift => {
                        // Block ranges of old basement ~100–200 km long between broad basins
                        // (the Wind River and Bighorn ranges), on a gentle regional rise.
                        let inside = interior_km.get(&c).copied().unwrap_or(0.0);
                        let s = smoothstep(0.0, BASEMENT_EDGE_KM.max(spacing_km), inside);
                        let blocks = smoothstep(0.45, 0.8, plateau_noise.ridged(grid.center(c), r / BASEMENT_BLOCK_KM, 2));
                        z = base_cont + relief * (BASEMENT_RISE_M + BASEMENT_RANGE_M * blocks) * s;
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
                    | G::Tessera
                    | G::RidgedPlains
                    | G::ChaosTerrain
                    | G::DilationalBand
                    | G::RidgeLinea
                    | G::GroovedTerrain => {}
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
                let d = bdist(c).unwrap_or(0.0);
                let z = -depth * (-(d / 35.0).powi(2)).exp().max(0.55);
                e[c as usize] = e[c as usize].min(z);
            }
        }

        // Continental margin zone: continental crust within slope + shelf width of the
        // continent–ocean edge. Shelves are built relative to sea level (below).
        // Slope and shelf widths (km): at an active margin the slope drops straight from a
        // shelf only a few tens of km wide into the trench (Chile's is 5–30 km, and the land
        // ~50 km from the 2 km isobath already stands 1–2 km high). Thickened
        // crust (ranges, plateaus) is not drowned into bays: its margin follows the true edge.
        let sw = p.slope_width_km;
        let margin_w = |c: Cell| -> (f64, f64, bool) {
            match provinces.owner_of(c).and_then(|k| provinces.table.get(k)).map(|p| p.kind) {
                Some(G::ExtendedCrust) => (sw, 380.0, true),
                Some(G::VolcanicArc | G::AccretionaryWedge | G::ForearcBasin) => (0.3 * sw, 15.0, true),
                Some(G::CollisionOrogen | G::OrogenicPlateau) => (0.5 * sw, 50.0, false),
                _ => (sw, 200.0, true),
            }
        };
        let margin: Vec<Option<(f64, f64, f64)>> = grid
            .cells()
            .map(|c| {
                let i = c as usize;
                let (slope, shelf, bays) = margin_w(c);
                let d = if bays { edge_km[i] } else { edge_raw[i] };
                (crust[i] == CrustType::Continental && d < slope + shelf).then_some((d, slope, shelf))
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
                let amp = if margin[i].is_some() {
                    60.0
                } else if land {
                    120.0 + 900.0 * rough
                } else {
                    120.0 + 350.0 * rough
                };
                let mut z = e[i] + rs * amp * base;
                if matches!(kind, Some(G::CollisionOrogen | G::OrogenicPlateau | G::VolcanicArc | G::PaleoOrogen | G::ContinentalRift)) {
                    let ridges = noise.ridged(pt, freq * 3.0, 6);
                    let mountain_amp = match kind {
                        Some(G::CollisionOrogen) => 2600.0 * prov.map_or(0.0, |p| 0.5 + 0.5 * p.intensity),
                        // Ranges and basins across the plateau (Tibet's ~1 km internal relief).
                        Some(G::OrogenicPlateau) => 1100.0,
                        Some(G::VolcanicArc) => {
                            if crust[i] == CrustType::Oceanic {
                                300.0
                            } else {
                                1800.0
                            }
                        }
                        Some(G::PaleoOrogen) => 700.0,
                        _ => 500.0,
                    };
                    z += rs * mountain_amp * (ridges - 0.35);
                }
                // An island arc is a chain of volcanoes, ~70 km apart, on its ridge: a few reach
                // the surface (the Marianas, Tonga, the Aleutians), not one smooth ridge.
                if kind == Some(G::VolcanicArc) && crust[i] == CrustType::Oceanic {
                    let cone = smoothstep(0.1, 0.55, noise.fbm(pt, r / ARC_VOLCANO_SPACING_KM, 2));
                    // Stretches with few or no volcanoes where the arc is weakest.
                    let seg = smoothstep(0.45, 1.2, along_arc(pt));
                    z += rs * (ARC_PEAK_M - ARC_RIDGE_M) * cone * seg * prov.map_or(0.6, |p| 0.6 + 0.4 * p.intensity);
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
                let Some((d, sw, shelf)) = *m else { continue };
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
        let field: Vec<f32> = e.iter().map(|z| (z - sea) as f32).collect();

        out.insert(ProductData::RoughElevation(Field::from_values(grid.spec(), field)));
        Ok(out)
    }
}

/// Depth (km) inside each orogenic plateau, from its outer edges (not the edge against its
/// front range or arc), so the plateau rises from its back and sides and meets the range at
/// full height. Cells more than twice the edge width inside read as that width.
pub(crate) fn plateau_interior_km(grid: &dyn Grid, provinces: &Partition<Province>, r: f64) -> HashMap<Cell, f64> {
    use GeologicProvince as G;
    let mut out = HashMap::new();
    for (key, cells) in &provinces.cells_by_owner() {
        let Some(prov) = provinces.table.get(*key) else { continue };
        if prov.kind != G::OrogenicPlateau {
            continue;
        }
        let set: std::collections::HashSet<Cell> = cells.iter().copied().collect();
        let outside = |n: Cell| !set.contains(&n) && provinces.row_at(n).is_none_or(|pr| !matches!(pr.kind, G::CollisionOrogen | G::VolcanicArc));
        let rim: Vec<Cell> = cells.iter().copied().filter(|&c| grid.neighbors(c).iter().any(outside)).collect();
        let d = sparse_dist(grid, &rim, 2.0 * PLATEAU_EDGE_KM / r);
        for &c in cells {
            out.insert(c, d.get(&c).map_or(2.0 * PLATEAU_EDGE_KM, |v| v * r));
        }
    }
    out
}

/// Plateau interiors deeper than this sit in their front range's rain shadow and drain into
/// their own basins: rivers cut only their rims.
pub(crate) const PLATEAU_INTERIOR_KM: f64 = PLATEAU_EDGE_KM;

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
        if vol(mid) < target {
            lo = mid
        } else {
            hi = mid
        }
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
    grid.neighbors(c).iter().map(|n| (f64::from(elev[n as usize]) - f64::from(elev[c as usize])).abs() / (pc.angle_to(grid.center(n)) * radius_km * 1000.0)).fold(0.0, f64::max)
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
