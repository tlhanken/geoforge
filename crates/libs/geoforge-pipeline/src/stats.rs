//! Summary statistics of a generated world, with Earth reference values alongside, for
//! checking realism.

use geoforge_core::grid::{AnyGrid, Cell, Grid};
use geoforge_core::product::Products;
use geoforge_core::shapes::Partition;
use geoforge_core::surface::Category;
use serde::Serialize;
use std::collections::BTreeMap;

#[derive(Debug, Clone, Serialize)]
pub struct Stat {
    pub value: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub earth: Option<f64>,
}

fn st(value: f64, earth: Option<f64>) -> Stat {
    Stat { value: (value * 1000.0).round() / 1000.0, earth }
}

/// Grouped statistics: section → name → value (with Earth reference where known).
pub type Summary = BTreeMap<String, BTreeMap<String, Stat>>;

#[allow(clippy::too_many_lines)]
#[must_use]
/// `sea`: the world has surface water, so ground below 0 m is sea; on a dry world every cell
/// is land.
pub fn summarize(pr: &Products, grid: &AnyGrid, radius_km: f64, sea: bool) -> Summary {
    let n = grid.len();
    let area: Vec<f64> = (0..n).map(|c| grid.area_sr(c as Cell)).collect();
    let total: f64 = area.iter().sum();
    let lat: Vec<f64> = (0..n).map(|c| grid.center(c as Cell).to_lat_lon().lat_deg).collect();
    let mut out: Summary = BTreeMap::new();
    let wmean = |vals: &dyn Fn(usize) -> Option<f64>| -> f64 {
        let (mut s, mut w) = (0.0, 0.0);
        for (c, a) in area.iter().enumerate() {
            if let Some(v) = vals(c) {
                s += v * a;
                w += a;
            }
        }
        if w > 0.0 {
            s / w
        } else {
            f64::NAN
        }
    };
    if let Some(pv) = pr.provinces() {
        // Share of the surface and number of regions per province kind.
        let mut area_by: BTreeMap<String, (f64, usize)> = BTreeMap::new();
        for (k, cells) in pv.cells_by_owner() {
            if let Some(row) = pv.table.get(k) {
                let e = area_by.entry(row.kind.name().to_lowercase().replace(' ', "_")).or_default();
                e.0 += cells.iter().map(|&c| area[c as usize]).sum::<f64>() / total;
                e.1 += 1;
            }
        }
        let mut s = BTreeMap::new();
        for (kind, (share, count)) in area_by {
            s.insert(format!("{kind}_of_surface"), st(share, None));
            s.insert(format!("{kind}_count"), st(count as f64, None));
        }
        out.insert("provinces".into(), s);
    }

    let Some(elev) = pr.elevation().map(|f| &f.values) else { return out };
    let land = |c: usize| !sea || elev[c] >= 0.0;
    let land_area: f64 = (0..n).filter(|&c| land(c)).map(|c| area[c]).sum();
    let km2 = |sr: f64| sr * radius_km * radius_km;

    let mut s = BTreeMap::new();
    s.insert("land_fraction".into(), st(land_area / total, Some(0.292)));
    s.insert("mean_land_elevation_m".into(), st(wmean(&|c| land(c).then(|| f64::from(elev[c]))), Some(840.0)));
    s.insert("mean_ocean_depth_m".into(), st(-wmean(&|c| (!land(c)).then(|| f64::from(elev[c]))), Some(3680.0)));
    s.insert("max_elevation_m".into(), st(elev.iter().copied().fold(f32::MIN, f32::max).into(), Some(8849.0)));
    s.insert("min_elevation_m".into(), st(elev.iter().copied().fold(f32::MAX, f32::min).into(), Some(-10935.0)));
    let above = |h: f32| (0..n).filter(|&c| elev[c] >= h).map(|c| area[c]).sum::<f64>() / land_area.max(1e-12);
    s.insert("land_above_1000m_fraction".into(), st(above(1000.0), Some(0.27)));
    s.insert("land_above_3000m_fraction".into(), st(above(3000.0), Some(0.05)));
    // Hypsometry about the mean radius: plains-dominated worlds (Venus) are narrow.
    let mean = wmean(&|c| Some(f64::from(elev[c])));
    let near = (0..n).filter(|&c| (f64::from(elev[c]) - mean).abs() <= 1000.0).map(|c| area[c]).sum::<f64>() / total;
    s.insert("within_1km_of_mean_fraction".into(), st(near, None));
    out.insert("elevation".into(), s);

    if let Some(crust) = pr.crust_type() {
        // Oceanic on plate worlds; the lowlands (thin crust) on rigid-lid worlds.
        let thin = (0..n).filter(|&c| crust.values[c] == geoforge_core::surface::CrustType::Oceanic).map(|c| area[c]).sum::<f64>() / total;
        out.insert("crust".into(), BTreeMap::from([("thin_crust_fraction".to_string(), st(thin, Some(0.6)))]));
    }

    if let Some(pl) = pr.plates() {
        let mut s = BTreeMap::new();
        let by = pl.cells_by_owner();
        let mut fr: Vec<f64> = by.values().map(|cs| cs.iter().map(|&c| area[c as usize]).sum::<f64>() / total).collect();
        fr.sort_by(|a, b| b.total_cmp(a));
        s.insert("count".into(), st(fr.len() as f64, None));
        s.insert("largest_fraction".into(), st(fr.first().copied().unwrap_or(0.0), Some(0.203)));
        s.insert("seven_largest_fraction".into(), st(fr.iter().take(7).sum(), Some(0.81)));
        out.insert("plates".into(), s);
    }

    if let Some(l) = pr.lithology() {
        let mut s = BTreeMap::new();
        for (slug, a) in area_by_kind(l, |r| &r.slug, &area, &land) {
            s.insert(slug + "_of_land", st(a / land_area.max(1e-12), None));
        }
        if let Some(age) = pr.rock_age().map(|f| &f.values) {
            let frac = |lo: f32, hi: f32| (0..n).filter(|&c| land(c) && age[c] >= lo && age[c] < hi).map(|c| area[c]).sum::<f64>() / land_area.max(1e-12);
            s.insert("land_rock_cenozoic_fraction".into(), st(frac(0.0, 66.0), Some(0.28)));
            s.insert("land_rock_mesozoic_fraction".into(), st(frac(66.0, 252.0), Some(0.22)));
            s.insert("land_rock_paleozoic_fraction".into(), st(frac(252.0, 541.0), Some(0.18)));
            s.insert("land_rock_precambrian_fraction".into(), st(frac(541.0, 5000.0), Some(0.32)));
        }
        out.insert("lithology".into(), s);
    }

    if let Some(t) = pr.temperature().map(|f| &f.values) {
        let mut s = BTreeMap::new();
        s.insert("global_mean_c".into(), st(wmean(&|c| Some(f64::from(t[c]))), Some(14.5)));
        s.insert("land_mean_c".into(), st(wmean(&|c| land(c).then(|| f64::from(t[c]))), Some(8.5)));
        for (lo, hi, e) in [(-5.0, 5.0, 26.5), (25.0, 35.0, 20.0), (40.0, 50.0, 10.0), (55.0, 65.0, 0.5), (75.0, 90.0, -16.0)] {
            s.insert(format!("zonal_{lo}_{hi}N_c"), st(wmean(&|c| (lat[c] >= lo && lat[c] < hi).then(|| f64::from(t[c]))), Some(e)));
        }
        s.insert("zonal_75_90S_c".into(), st(wmean(&|c| (lat[c] <= -75.0).then(|| f64::from(t[c]))), Some(-40.0)));
        s.insert("coldest_c".into(), st(t.iter().copied().fold(f32::MAX, f32::min).into(), Some(-58.0)));
        s.insert("warmest_c".into(), st(t.iter().copied().fold(f32::MIN, f32::max).into(), Some(35.0)));
        if let Some(r) = pr.temperature_range().map(|f| &f.values) {
            s.insert("range_land_40_60_c".into(), st(wmean(&|c| (land(c) && lat[c].abs() >= 40.0 && lat[c].abs() < 60.0).then(|| f64::from(r[c]))), Some(28.0)));
            s.insert("range_ocean_40_60_c".into(), st(wmean(&|c| (!land(c) && lat[c].abs() >= 40.0 && lat[c].abs() < 60.0).then(|| f64::from(r[c]))), Some(7.0)));
            s.insert("range_tropics_land_c".into(), st(wmean(&|c| (land(c) && lat[c].abs() < 10.0).then(|| f64::from(r[c]))), Some(3.0)));
            s.insert("range_max_c".into(), st(r.iter().copied().fold(f32::MIN, f32::max).into(), Some(62.0)));
        }
        out.insert("temperature".into(), s);
    }

    if let Some(p) = pr.precipitation().map(|f| &f.values) {
        let mut s = BTreeMap::new();
        s.insert("global_mean_mm".into(), st(wmean(&|c| Some(f64::from(p[c]))), Some(1000.0)));
        s.insert("land_mean_mm".into(), st(wmean(&|c| land(c).then(|| f64::from(p[c]))), Some(780.0)));
        s.insert("ocean_mean_mm".into(), st(wmean(&|c| (!land(c)).then(|| f64::from(p[c]))), Some(1100.0)));
        let frac = |f: &dyn Fn(f32) -> bool| (0..n).filter(|&c| land(c) && f(p[c])).map(|c| area[c]).sum::<f64>() / land_area.max(1e-12);
        s.insert("land_below_250mm_fraction".into(), st(frac(&|v| v < 250.0), Some(0.2)));
        s.insert("land_above_2000mm_fraction".into(), st(frac(&|v| v >= 2000.0), Some(0.11)));
        for (lo, hi, e) in [(-5.0, 10.0, 2000.0), (20.0, 30.0, 800.0), (40.0, 55.0, 1000.0), (70.0, 90.0, 250.0)] {
            s.insert(format!("zonal_{lo}_{hi}N_mm"), st(wmean(&|c| (lat[c] >= lo && lat[c] < hi).then(|| f64::from(p[c]))), Some(e)));
        }
        out.insert("precipitation".into(), s);
    }

    if let (Some(rv), Some(lk)) = (pr.rivers(), pr.lakes()) {
        let mut s = BTreeMap::new();
        // Discharge into the sea: edges whose downstream node is a sea cell (on a world without
        // one, a lake: the rivers' ends).
        let mut to_sea = 0.0;
        let mut mouths: Vec<f64> = Vec::new();
        for e in &rv.edges {
            if let Some(to) = rv.nodes.get(e.to) {
                if if sea { elev[to.cell as usize] < 0.0 } else { lk.owner.values[to.cell as usize] != 0 } {
                    to_sea += e.attrs.discharge_m3_s;
                    mouths.push(e.attrs.discharge_m3_s);
                }
            }
        }
        mouths.sort_by(|a, b| b.total_cmp(a));
        s.insert("river_discharge_to_sea_km3_yr".into(), st(to_sea * 3.155_76e7 / 1e9, Some(37_000.0)));
        s.insert("largest_river_m3_s".into(), st(mouths.first().copied().unwrap_or(0.0), Some(209_000.0)));
        s.insert("second_river_m3_s".into(), st(mouths.get(1).copied().unwrap_or(0.0), Some(41_000.0)));
        s.insert("rivers_over_10000_m3_s".into(), st(mouths.iter().filter(|&&q| q >= 10_000.0).count() as f64, Some(8.0)));
        s.insert("rivers_over_1000_m3_s".into(), st(mouths.iter().filter(|&&q| q >= 1000.0).count() as f64, Some(80.0)));
        s.insert("max_strahler".into(), st(rv.edges.iter().map(|e| e.attrs.strahler_order).max().unwrap_or(0).into(), None));
        let by = lk.cells_by_owner();
        let mut la: Vec<(f64, bool)> = by.iter().map(|(k, cs)| (km2(cs.iter().map(|&c| area[c as usize]).sum()), lk.table.get(*k).is_some_and(|l| l.endorheic))).collect();
        la.sort_by(|a, b| b.0.total_cmp(&a.0));
        s.insert("lakes".into(), st(la.len() as f64, None));
        s.insert("lakes_over_10000_km2".into(), st(la.iter().filter(|l| l.0 >= 10_000.0).count() as f64, Some(16.0)));
        s.insert("largest_lake_km2".into(), st(la.first().map_or(0.0, |l| l.0), Some(371_000.0)));
        s.insert("lake_fraction_of_land".into(), st(km2(by.values().flatten().map(|&c| area[c as usize]).sum::<f64>()) / km2(land_area).max(1e-9), Some(0.018)));
        s.insert("endorheic_lakes".into(), st(la.iter().filter(|l| l.1).count() as f64, None));
        if let Some(ws) = pr.watersheds() {
            let by = ws.cells_by_owner();
            let mut wa: Vec<(f64, bool)> = by.iter().map(|(k, cs)| (km2(cs.iter().map(|&c| area[c as usize]).sum()), ws.table.get(*k).is_some_and(|w| w.reaches_sea))).collect();
            wa.sort_by(|a, b| b.0.total_cmp(&a.0));
            s.insert("largest_basin_km2".into(), st(wa.first().map_or(0.0, |w| w.0), Some(7_000_000.0)));
            let endo: f64 = wa.iter().filter(|w| !w.1).map(|w| w.0).sum();
            s.insert("endorheic_fraction_of_land".into(), st(endo / km2(land_area).max(1e-9), Some(0.18)));
        }
        out.insert("hydrology".into(), s);
    }

    if let Some(b) = pr.biome() {
        let mut s = BTreeMap::new();
        // Earth's land cover by biome (share of land).
        let earth: BTreeMap<&str, f64> = [
            ("ice_sheet", 0.105),
            ("tundra", 0.06),
            ("boreal_forest", 0.11),
            ("temperate_deciduous_forest", 0.08),
            ("temperate_grassland", 0.08),
            ("mediterranean_woodland", 0.02),
            ("hot_desert", 0.14),
            ("cold_desert", 0.05),
            ("savanna", 0.15),
            ("tropical_rainforest", 0.07),
            ("tropical_seasonal_forest", 0.05),
        ]
        .into_iter()
        .collect();
        for (slug, a) in area_by_kind(b, |r| &r.slug, &area, &|_| true) {
            if slug != "ocean" {
                let reference = earth.get(slug.as_str()).copied();
                s.insert(slug + "_of_land", st(a / land_area.max(1e-12), reference));
            }
        }
        out.insert("biomes".into(), s);
    }
    out
}

/// Area of each catalog kind (by slug) over the cells `keep` accepts; every row is listed,
/// used or not.
fn area_by_kind<E>(pt: &Partition<E>, slug: impl Fn(&E) -> &String, area: &[f64], keep: &dyn Fn(usize) -> bool) -> BTreeMap<String, f64> {
    let mut by_key: BTreeMap<u32, f64> = BTreeMap::new();
    for (c, &o) in pt.owner.values.iter().enumerate() {
        if keep(c) {
            *by_key.entry(o).or_default() += area[c];
        }
    }
    pt.table.rows.iter().map(|(k, r)| (slug(r).clone(), by_key.get(&(k.0 as u32)).copied().unwrap_or(0.0))).collect()
}
