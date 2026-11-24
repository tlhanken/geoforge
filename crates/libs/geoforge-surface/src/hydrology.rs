//! Hydrology layer: river network, lakes and drainage basins.
//!
//! 1. **Runoff** from the climate's water balance: actual evaporation follows Budyko's curve
//!    (Fu's form), so arid land yields almost nothing and wet land most of its rain.
//! 2. **Drainage**: priority-flood from the sea (Barnes et al. 2014). Every land cell gets a
//!    downstream neighbour, depressions included, in one pass.
//! 3. **Lakes**: each closed depression is tested with its water balance. If inflow plus
//!    rain on the lake beats evaporation from the full lake, it fills to its sill and
//!    overflows. If not, it is endorheic: the lake shrinks until evaporation matches inflow,
//!    and nothing flows on (Caspian, Great Salt Lake, Chad); a basin with too little water is
//!    left as a dry pan. Overflowing lakes only last where something keeps the basin open
//!    against the river cutting its sill and sediment filling it: active rifting (Baikal,
//!    Tanganyika), great depth, or recent glaciation (cold regions). Elsewhere the
//!    depression is treated as an alluvial plain the river runs straight through.
//!    Depressions under ice (very cold) hold no lakes.
//! 4. **Rivers**: discharge accumulated downstream; reaches above a threshold form the
//!    network, ordered by Strahler number.
//! 5. **Watersheds**: each cell's terminal outlet, the sea or a closed lake.

use crate::climate::pet_mm;
use geoforge_core::grid::{Cell, Grid};
use geoforge_core::key::Key;
use geoforge_core::product::{Generator, LayerId, Precondition, ProductData, ProductId, Products, RunContext};
use geoforge_core::shapes::{Edge, Field, Network, Partition, Table};
use geoforge_core::surface::{Lake, RiverNode, RiverReach, Watershed};
use geoforge_core::{BodyContext, Error};
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::collections::BinaryHeap;

/// Share of the river threshold a cell must receive from upstream to carry a river.
const RIVER_INFLOW_SHARE: f64 = 0.5;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct HydrologyParams {
    /// Shallowest depression kept as a lake basin, m.
    pub min_lake_depth_m: f64,
    /// Overflowing lakes at least this deep survive anywhere, m.
    pub deep_lake_m: f64,
    /// Overflowing lakes survive where the mean temperature is below this (glacial terrain), °C.
    pub glacial_lake_max_c: f64,
    /// No lakes where the mean temperature is below this (ice-covered), °C.
    pub ice_c: f64,
    /// Smallest mean discharge drawn as a river, m³/s.
    pub min_river_m3_s: f64,
    /// Smallest drainage basin given its own watershed key, km².
    pub min_basin_km2: f64,
    /// Fu's Budyko parameter (higher = more evaporation for the same aridity).
    pub budyko_w: f64,
    /// Open-water evaporation relative to land potential evapotranspiration.
    pub lake_evaporation_factor: f64,
    /// Fraction of a river's flow lost per 1000 km through fully arid land (seepage and
    /// evaporation from the channel, as on the Nile or Niger).
    pub arid_channel_loss_per_1000km: f64,
    /// Largest area an overflowing lake keeps; bigger basins are taken as partly filled with
    /// sediment and keep only their deepest part, km².
    pub max_open_lake_km2: f64,
}

impl Default for HydrologyParams {
    fn default() -> Self {
        Self {
            min_lake_depth_m: 15.0,
            deep_lake_m: 300.0,
            glacial_lake_max_c: 3.0,
            ice_c: -10.0,
            min_river_m3_s: 50.0,
            min_basin_km2: 20_000.0,
            budyko_w: 3.0,
            lake_evaporation_factor: 1.4,
            arid_channel_loss_per_1000km: 0.35,
            max_open_lake_km2: 80_000.0,
        }
    }
}

/// Mean annual runoff (mm/yr) from precipitation and potential evaporation (Fu 1981).
#[must_use]
pub fn runoff_mm(p: f64, pet: f64, w: f64) -> f64 {
    if p <= 0.0 {
        return 0.0;
    }
    let phi = pet / p;
    let aet = p * (1.0 + phi - (1.0 + phi.powf(w)).powf(1.0 / w));
    (p - aet).max(0.0)
}

#[derive(PartialEq)]
struct Item(f64, u64, Cell);
impl Eq for Item {}
impl Ord for Item {
    fn cmp(&self, o: &Self) -> Ordering {
        // Min-heap on level, FIFO among equals.
        o.0.partial_cmp(&self.0).unwrap_or(Ordering::Equal).then_with(|| o.1.cmp(&self.1))
    }
}
impl PartialOrd for Item {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}

/// Drainage of a surface: downstream neighbour per land cell, filled levels and the order
/// cells were reached (downstream first).
pub struct Drainage {
    pub parent: Vec<u32>,
    pub filled: Vec<f64>,
    pub order: Vec<Cell>,
}

pub const NONE: u32 = u32::MAX;

/// Priority-flood from the sea. Cells with `elev < 0` are sea; `parent` of sea cells is NONE.
#[must_use]
pub fn drainage(grid: &dyn Grid, elev: &[f32]) -> Drainage {
    drainage_with(grid, elev, None)
}

/// Coherent micro-relief in [0, 1] for routing water across flats (sediment plains, filled
/// basins): without it, flow on a flat follows the grid in straight parallel lines.
#[must_use]
pub fn flat_jitter(grid: &dyn Grid, seed: u64) -> Vec<f32> {
    let noise = crate::noise::TerrainNoise::new(seed);
    (0..grid.len() as Cell).map(|c| (0.5 + noise.fbm(grid.center(c), 30.0, 4)).clamp(0.0, 1.0) as f32).collect()
}

/// [`drainage`] with optional flat-routing jitter (see [`flat_jitter`]).
#[must_use]
pub fn drainage_with(grid: &dyn Grid, elev: &[f32], jitter: Option<&[f32]>) -> Drainage {
    let n = grid.len();
    let mut parent = vec![NONE; n];
    let mut filled: Vec<f64> = elev.iter().map(|&e| f64::from(e)).collect();
    let mut seen = vec![false; n];
    let mut heap = BinaryHeap::new();
    let mut tick = 0u64;
    for c in grid.cells() {
        if elev[c as usize] < 0.0 {
            seen[c as usize] = true;
            if grid.neighbors(c).iter().any(|nb| elev[nb as usize] >= 0.0) {
                heap.push(Item(filled[c as usize], tick, c));
                tick += 1;
            }
        }
    }
    let mut order = Vec::new();
    // No land touches the sea (a world of one ocean-free continent): drain to the lowest cell.
    if heap.is_empty() {
        if let Some(c) = grid.cells().min_by(|&a, &b| elev[a as usize].total_cmp(&elev[b as usize])) {
            seen[c as usize] = true;
            heap.push(Item(filled[c as usize], 0, c));
            order.push(c);
        }
    }
    while let Some(Item(level, _, c)) = heap.pop() {
        for nb in grid.neighbors(c).iter() {
            let i = nb as usize;
            if seen[i] {
                continue;
            }
            seen[i] = true;
            parent[i] = c;
            let eps = 1e-3 * (1.0 + 20.0 * jitter.map_or(0.0, |j| f64::from(j[i])));
            filled[i] = f64::from(elev[i]).max(level + eps);
            order.push(nb);
            heap.push(Item(filled[i], tick, nb));
            tick += 1;
        }
    }
    Drainage { parent, filled, order }
}

pub struct Hydrology;

struct Basin {
    cells: Vec<Cell>,
    remaining: usize,
    ext_inflow: f64,
    out: f64,
    exit: u32,
}

impl Generator for Hydrology {
    type Params = HydrologyParams;
    const NAME: &'static str = "hydrology";
    const VERSION: u32 = 3;
    const LAYER: LayerId = LayerId::Hydrology;

    fn needs() -> &'static [ProductId] {
        &[ProductId::Elevation, ProductId::Precipitation, ProductId::Temperature]
    }
    fn wants() -> &'static [ProductId] {
        &[ProductId::Provinces, ProductId::RoughElevation]
    }
    fn applies(body: &BodyContext) -> Result<(), Precondition> {
        crate::preconditions::surface_water(body)
    }

    #[allow(clippy::too_many_lines)]
    fn run(ctx: &RunContext<'_>, inputs: &Products, p: &HydrologyParams) -> Result<Products, Error> {
        if inputs.rivers().is_some() && inputs.lakes().is_some() && inputs.watersheds().is_some() {
            return Ok(Products::new());
        }
        let elev = &inputs.elevation().ok_or(Error::MissingInput(ProductId::Elevation))?.values;
        let precip = &inputs.precipitation().ok_or(Error::MissingInput(ProductId::Precipitation))?.values;
        let temp = &inputs.temperature().ok_or(Error::MissingInput(ProductId::Temperature))?.values;
        let temp = crate::climate::on_ground(temp, inputs.rough_elevation().map(|f| f.values.as_slice()), elev, ctx.body);
        let pet: Vec<f64> = temp.iter().map(|&t| pet_mm(f64::from(t))).collect();
        let precip: Vec<f64> = precip.iter().map(|&v| f64::from(v)).collect();
        let water = Water { elev, precip: &precip, pet: &pet, lake_temp: Some(&temp), sea: true, inventory_km3: None };
        Ok(drain(ctx, inputs, &water, p))
    }
}

/// What a world's liquid does on its ground: the heights, the rain and the evaporative demand
/// (mm/yr of the liquid), and the mean temperature that freezes or keeps lakes (`None`: a
/// liquid whose lakes neither freeze nor fill, every closed depression keeps its lake). `sea`:
/// ground below datum is sea; without one the world drains to its lowest point. A world with
/// little liquid fills only as many lakes as its `inventory_km3` holds.
pub(crate) struct Water<'a> {
    pub elev: &'a [f32],
    pub precip: &'a [f64],
    pub pet: &'a [f64],
    pub lake_temp: Option<&'a [f32]>,
    pub sea: bool,
    /// All the liquid the lakes may hold, km³ (`None`: no limit).
    pub inventory_km3: Option<f64>,
}

/// Rivers, lakes and watersheds of `w`, those `inputs` lack.
#[allow(clippy::too_many_lines)]
pub(crate) fn drain(ctx: &RunContext<'_>, inputs: &Products, w: &Water, p: &HydrologyParams) -> Products {
    let mut out = Products::new();
    let grid = ctx.grid;
    let r = ctx.radius_km();
    let n = grid.len();
    // Without a sea, heights from just under the lowest point, which takes what reaches it.
    let datum = if w.sea { 0.0 } else { w.elev.iter().copied().fold(f32::INFINITY, f32::min) + 1e-3 };
    let shifted: Vec<f32>;
    let elev: &[f32] = if w.sea {
        w.elev
    } else {
        shifted = w.elev.iter().map(|&e| e - datum).collect();
        &shifted
    };
    let (precip, pet) = (w.precip, w.pet);
    let area_km2: Vec<f64> = (0..n).map(|c| grid.area_sr(c as Cell) * r * r).collect();
    let runoff: Vec<f64> = (0..n).map(|c| runoff_mm(precip[c], pet[c], p.budyko_w)).collect();
    let e_lake: Vec<f64> = pet.iter().map(|v| v * p.lake_evaporation_factor).collect();
    // Fraction of flow kept crossing each cell (channel losses where the land is arid).
    let keep: Vec<f64> = (0..n)
        .map(|c| {
            let aridity = if pet[c] > 0.0 { (1.0 - precip[c] / pet[c]).clamp(0.0, 1.0) } else { 0.0 };
            let len_km = (area_km2[c]).sqrt();
            (1.0 - p.arid_channel_loss_per_1000km).max(1e-6).powf(aridity * len_km / 1000.0)
        })
        .collect();
    let jitter = flat_jitter(grid, ctx.seed(ProductId::Rivers).value());
    let d = drainage_with(grid, elev, Some(&jitter));
    let land = |c: usize| elev[c] >= 0.0;

    // Depressions deep enough to hold a lake.
    let raised: Vec<bool> = (0..n).map(|c| land(c) && d.filled[c] - f64::from(elev[c]) > 0.5).collect();
    let mut comp = vec![u32::MAX; n];
    let mut basins: Vec<Basin> = Vec::new();
    let mut open_ok: Vec<bool> = Vec::new();
    let rift_cell = |c: Cell| inputs.provinces().is_some_and(|pv| pv.owner_of(c).and_then(|k| pv.table.get(k)).is_some_and(|pr| pr.kind == geoforge_core::surface::GeologicProvince::ContinentalRift));
    for s in 0..n {
        if !raised[s] || comp[s] != u32::MAX {
            continue;
        }
        let id = basins.len() as u32;
        let mut stack = vec![s as Cell];
        comp[s] = id;
        let mut cells = Vec::new();
        while let Some(c) = stack.pop() {
            cells.push(c);
            for nb in grid.neighbors(c).iter() {
                if raised[nb as usize] && comp[nb as usize] == u32::MAX {
                    comp[nb as usize] = id;
                    stack.push(nb);
                }
            }
        }
        let depth = cells.iter().map(|&c| d.filled[c as usize] - f64::from(elev[c as usize])).fold(0.0, f64::max);
        let t_mean = w.lake_temp.map(|t| cells.iter().map(|&c| f64::from(t[c as usize])).sum::<f64>() / cells.len() as f64);
        let rift = cells.iter().any(|&c| rift_cell(c));
        let persistent = rift || depth >= p.deep_lake_m || t_mean.is_none_or(|t| t < p.glacial_lake_max_c);
        if depth < p.min_lake_depth_m || t_mean.is_some_and(|t| t < p.ice_c) {
            for &c in &cells {
                comp[c as usize] = u32::MAX - 1; // shallow or frozen: ordinary land
            }
        }
        open_ok.push(persistent);
        basins.push(Basin { remaining: cells.len(), cells, ext_inflow: 0.0, out: 0.0, exit: NONE });
    }
    let in_basin = |c: usize| -> Option<usize> {
        if comp[c] < u32::MAX - 1 {
            Some(comp[c] as usize)
        } else {
            None
        }
    };

    // Accumulate upstream to downstream, resolving each lake when its last cell is done.
    // Units: mm/yr × km² = 1000 m³/yr.
    let mut acc = vec![0.0f64; n];
    let mut lake_cells: Vec<(usize, Vec<Cell>, bool, f64)> = Vec::new(); // (basin, cells, endorheic, level)
    let mut endorheic = vec![false; basins.len()];
    for &c in d.order.iter().rev() {
        let ci = c as usize;
        let b = in_basin(ci);
        acc[ci] += match b {
            Some(_) => (precip[ci] - e_lake[ci]) * area_km2[ci],
            None => runoff[ci] * area_km2[ci],
        };
        let par = d.parent[ci];
        if par == NONE {
            continue;
        }
        let pb = in_basin(par as usize);
        match b {
            Some(k) if pb != Some(k) => {
                basins[k].out += acc[ci];
                basins[k].exit = par;
            }
            _ => {
                if b.is_none() && acc[ci] > 0.0 {
                    acc[ci] *= keep[ci];
                }
                if let (Some(k), true) = (pb, b.is_none()) {
                    basins[k].ext_inflow += acc[ci].max(0.0);
                }
                acc[par as usize] += acc[ci];
            }
        }
        if let Some(k) = b {
            basins[k].remaining -= 1;
            if basins[k].remaining == 0 {
                let bk = &basins[k];
                let level = bk.cells.iter().map(|&x| d.filled[x as usize]).fold(f64::NEG_INFINITY, f64::max);
                if bk.out >= 0.0 {
                    let (exit, outflow) = (bk.exit, bk.out);
                    if exit != NONE {
                        acc[exit as usize] += outflow;
                        if let Some(k2) = in_basin(exit as usize) {
                            if k2 != k {
                                basins[k2].ext_inflow += outflow;
                            }
                        }
                    }
                    if open_ok[k] {
                        // Large basins keep only their deepest part as open water.
                        let mut cells = basins[k].cells.clone();
                        cells.sort_by(|a, b| elev[*a as usize].total_cmp(&elev[*b as usize]));
                        let mut a = 0.0;
                        cells.retain(|&x| {
                            a += area_km2[x as usize];
                            a <= p.max_open_lake_km2 || a - area_km2[x as usize] == 0.0
                        });
                        lake_cells.push((k, cells, false, level));
                    }
                } else {
                    endorheic[k] = true;
                    // Grow the lake from its lowest cell while inflow still covers evaporation.
                    let mut cells = bk.cells.clone();
                    cells.sort_by(|a, b| elev[*a as usize].total_cmp(&elev[*b as usize]));
                    let dry: f64 = cells.iter().map(|&x| runoff[x as usize] * area_km2[x as usize]).sum();
                    let mut bal = bk.ext_inflow + dry;
                    let mut j = 0;
                    for &x in &cells {
                        let xi = x as usize;
                        let next = bal - runoff[xi] * area_km2[xi] + (precip[xi] - e_lake[xi]) * area_km2[xi];
                        if next < 0.0 {
                            break;
                        }
                        bal = next;
                        j += 1;
                    }
                    cells.truncate(j);
                    if !cells.is_empty() {
                        let lvl = f64::from(elev[*cells.last().expect("non-empty") as usize]);
                        lake_cells.push((k, cells, true, lvl));
                    }
                }
            }
        }
    }
    let q = |c: usize| acc[c].max(0.0) * 1000.0 / 3.155_76e7;

    // Lakes partition, largest first.
    let mut lakes_sorted: Vec<_> = lake_cells.iter().map(|(k, cs, endo, lvl)| (cs.iter().map(|&c| area_km2[c as usize]).sum::<f64>(), *k, cs, *endo, *lvl)).collect();
    // A world with only so much liquid fills its lakes where evaporation asks least, the
    // coldest first, until it runs out; the rest are dry pans.
    if let Some(inventory_km3) = w.inventory_km3 {
        let demand = |cs: &[Cell]| cs.iter().map(|&c| pet[c as usize]).sum::<f64>() / cs.len().max(1) as f64;
        let volume = |cs: &[Cell], lvl: f64| cs.iter().map(|&c| (lvl - f64::from(elev[c as usize])).max(0.0) / 1000.0 * area_km2[c as usize]).sum::<f64>();
        lakes_sorted.sort_by(|a, b| demand(a.2).total_cmp(&demand(b.2)).then(a.1.cmp(&b.1)));
        let mut held = 0.0;
        lakes_sorted.retain(|l| {
            held += volume(l.2, l.4);
            held <= inventory_km3
        });
    }
    lakes_sorted.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)));
    let mut lake_owner = vec![0u32; n];
    let mut lake_table = Table::new();
    for (i, (_, _, cs, endo, lvl)) in lakes_sorted.iter().enumerate() {
        let key = i as u32 + 1;
        for &c in cs.iter() {
            lake_owner[c as usize] = key;
        }
        lake_table.insert(Key(u64::from(key)), Lake { surface_elevation_m: *lvl + f64::from(datum), endorheic: *endo });
    }

    // Rivers: land cells (not lake) above the threshold, plus their mouths. A cell's own
    // runoff runs off in many small streams, so a river also needs water from upstream:
    // in a very wet land every cell would otherwise count as one.
    let own = |c: usize| runoff[c].max(0.0) * area_km2[c] * 1000.0 / 3.155_76e7;
    let is_river = |c: usize| land(c) && lake_owner[c] == 0 && q(c) >= p.min_river_m3_s && q(c) - own(c) >= RIVER_INFLOW_SHARE * p.min_river_m3_s;
    let mut strahler = vec![0u8; n];
    let mut best = vec![(0u8, 0u8); n]; // (max upstream order, how many have it)
    for &c in d.order.iter().rev() {
        let ci = c as usize;
        if is_river(ci) || lake_owner[ci] != 0 {
            let (m, k) = best[ci];
            strahler[ci] = if lake_owner[ci] != 0 {
                m
            } else if m == 0 {
                1
            } else if k >= 2 {
                m.saturating_add(1)
            } else {
                m
            };
            let par = d.parent[ci];
            if par != NONE && strahler[ci] > 0 {
                let e = &mut best[par as usize];
                match strahler[ci].cmp(&e.0) {
                    Ordering::Greater => *e = (strahler[ci], 1),
                    Ordering::Equal => e.1 = e.1.saturating_add(1),
                    Ordering::Less => {}
                }
            }
        }
    }
    let mut nodes = Table::new();
    let mut edges = Vec::new();
    for (c, &order) in strahler.iter().enumerate() {
        if !is_river(c) {
            continue;
        }
        nodes.insert(Key(c as u64 + 1), RiverNode { cell: c as Cell, discharge_m3_s: q(c) });
        let par = d.parent[c];
        if par == NONE {
            continue;
        }
        let pi = par as usize;
        if !is_river(pi) {
            nodes.insert(Key(pi as u64 + 1), RiverNode { cell: par, discharge_m3_s: q(c) });
        }
        edges.push(Edge { from: Key(c as u64 + 1), to: Key(pi as u64 + 1), attrs: RiverReach { strahler_order: order.max(1), discharge_m3_s: q(c) } });
    }

    // Watersheds: terminal outlet per cell, downstream first.
    let mut root = vec![NONE; n];
    let mut roots: Vec<(u32, bool)> = Vec::new(); // (outlet cell, reaches sea)
    for &c in &d.order {
        let ci = c as usize;
        if !land(ci) {
            continue;
        }
        let par = d.parent[ci];
        let b = in_basin(ci);
        let new_root = match (par, b) {
            (NONE, _) => Some((c, false)),
            (pp, _) if !land(pp as usize) => Some((pp, true)),
            (pp, Some(k)) if endorheic[k] && in_basin(pp as usize) != Some(k) => {
                let low = basins[k].cells.iter().copied().min_by(|a, b| elev[*a as usize].total_cmp(&elev[*b as usize])).unwrap_or(c);
                Some((low, false))
            }
            _ => None,
        };
        root[ci] = match new_root {
            Some(rt) => {
                roots.push(rt);
                (roots.len() - 1) as u32
            }
            None => root[par as usize],
        };
    }
    let mut basin_area = vec![0.0; roots.len()];
    for c in 0..n {
        if root[c] != NONE {
            basin_area[root[c] as usize] += area_km2[c];
        }
    }
    let mut big: Vec<usize> = (0..roots.len()).filter(|&i| basin_area[i] >= p.min_basin_km2).collect();
    big.sort_by(|a, b| basin_area[*b].total_cmp(&basin_area[*a]).then(a.cmp(b)));
    let mut ws_key = vec![0u32; roots.len()];
    let mut ws_table = Table::new();
    for (i, &b) in big.iter().enumerate() {
        ws_key[b] = i as u32 + 1;
        ws_table.insert(Key(i as u64 + 1), Watershed { outlet: roots[b].0, reaches_sea: w.sea && roots[b].1 });
    }
    let ws_owner: Vec<u32> = (0..n).map(|c| if root[c] == NONE { 0 } else { ws_key[root[c] as usize] }).collect();

    let spec = grid.spec();
    if inputs.rivers().is_none() {
        out.insert(ProductData::Rivers(Network { nodes, edges }));
    }
    if inputs.lakes().is_none() {
        out.insert(ProductData::Lakes(Partition { owner: Field::from_values(spec, lake_owner), table: lake_table }));
    }
    if inputs.watersheds().is_none() {
        out.insert(ProductData::Watersheds(Partition { owner: Field::from_values(spec, ws_owner), table: ws_table }));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn budyko_limits() {
        // Energy-limited (wet, cold): most rain runs off. Water-limited (hot, dry): almost none.
        assert!(runoff_mm(2000.0, 300.0, 2.6) > 1650.0);
        assert!(runoff_mm(200.0, 2000.0, 2.6) < 5.0);
        // Earth-like humid temperate: ~30–45% runoff.
        let r = runoff_mm(900.0, 750.0, 2.6) / 900.0;
        assert!((0.25..0.5).contains(&r), "{r}");
    }
}
