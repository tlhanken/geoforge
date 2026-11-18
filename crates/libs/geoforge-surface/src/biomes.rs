//! Biomes layer: a vegetation/land-cover class per cell from climate, with lakes and
//! wetlands from hydrology.
//!
//! The classification follows Köppen's thresholds where they are physical (the 10 °C
//! warmest-month tree line, the 18 °C coldest-month tropics, the aridity threshold that
//! scales with temperature) and Whittaker's precipitation steps within each thermal zone.
//! Monthly extremes are taken as mean ± range/2.

use crate::climate::pet_mm;
use geoforge_core::grid::Cell;
use geoforge_core::product::{Generator, LayerId, Precondition, ProductData, ProductId, Products, RunContext};
use geoforge_core::shapes::Field;
use geoforge_core::surface::Biome;
use geoforge_core::{BodyContext, Error};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct BiomesParams {
    /// Elevation above which cold climates count as alpine rather than polar, m.
    pub alpine_min_m: f64,
    /// River discharge above which flat lowland floodplains become wetland, m³/s.
    pub wetland_discharge_m3_s: f64,
    /// Steepest slope (rise/run) still counted as a floodplain.
    pub wetland_max_slope: f64,
}

impl Default for BiomesParams {
    fn default() -> Self {
        Self { alpine_min_m: 1800.0, wetland_discharge_m3_s: 2500.0, wetland_max_slope: 4e-4 }
    }
}

/// Climate of one land cell.
#[derive(Debug, Clone, Copy)]
pub struct CellClimate {
    pub t_mean: f64,
    pub t_range: f64,
    pub precip_mm: f64,
    pub elevation_m: f64,
    pub lat_deg: f64,
    /// Share of the year's precipitation that falls in the warm half-year (0–1).
    pub warm_share: f64,
}

/// Köppen's aridity threshold (mm/yr) for mean temperature `t`: higher when rain falls in
/// the warm season (more of it evaporates), lower when it falls in winter.
#[must_use]
pub fn aridity_threshold(t: f64, warm_share: f64) -> f64 {
    (20.0 * t + 140.0 + 140.0 * ((warm_share - 0.5) / 0.2).clamp(-1.0, 1.0)).max(60.0)
}

/// Classify a land cell.
#[must_use]
pub fn classify(c: &CellClimate, p: &BiomesParams) -> Biome {
    let warm = c.t_mean + c.t_range / 2.0;
    let cold = c.t_mean - c.t_range / 2.0;
    // Alpine: cold because of altitude. At sea level the same place would be warm enough
    // for trees, so it is a mountain climate rather than a polar one.
    let warm_at_sea_level = c.t_mean + c.t_range / 2.0 + 6.5 * c.elevation_m.max(0.0) / 1000.0;
    let high = c.elevation_m >= p.alpine_min_m && warm_at_sea_level >= 10.0;
    if warm < 0.0 {
        return if high { Biome::Alpine } else { Biome::IceSheet };
    }
    if warm < 10.0 {
        return if high { Biome::Alpine } else { Biome::Tundra };
    }
    let th = aridity_threshold(c.t_mean, c.warm_share);
    let pr = c.precip_mm;
    if pr < 0.5 * th {
        return if c.t_mean >= 18.0 { Biome::HotDesert } else { Biome::ColdDesert };
    }
    let semi_arid = pr < th;
    if cold >= 18.0 {
        // Tropical.
        return if semi_arid {
            Biome::Savanna
        } else if pr >= 1800.0 && c.t_range < 6.0 || pr >= 2300.0 {
            Biome::TropicalRainforest
        } else if pr >= 1300.0 {
            Biome::TropicalSeasonalForest
        } else {
            Biome::Savanna
        };
    }
    if cold < -3.0 && warm < 22.0 {
        // Cold-winter, cool-summer: boreal, or steppe where too dry for trees.
        return if semi_arid { Biome::TemperateGrassland } else { Biome::BorealForest };
    }
    // Temperate (and hot-summer continental).
    if semi_arid {
        return if c.t_mean >= 18.0 { Biome::Savanna } else { Biome::TemperateGrassland };
    }
    // Dry summers, wet mild winters (Köppen Cs).
    if c.warm_share < 0.35 && cold > -3.0 && pr < 1200.0 {
        return Biome::Mediterranean;
    }
    if pr >= 2000.0 || (pr >= 1500.0 && c.t_range < 14.0 && c.t_mean < 16.0) {
        return Biome::TemperateRainforest;
    }
    if pr < 1.6 * th {
        return Biome::TemperateGrassland;
    }
    if c.t_mean >= 20.0 && cold >= 10.0 {
        return if pr >= 1300.0 { Biome::TropicalSeasonalForest } else { Biome::Savanna };
    }
    Biome::TemperateForest
}

pub struct Biomes;

impl Generator for Biomes {
    type Params = BiomesParams;
    const LAYER: LayerId = LayerId::Biomes;

    fn needs() -> &'static [ProductId] {
        &[ProductId::Elevation, ProductId::Temperature, ProductId::TemperatureRange, ProductId::Precipitation]
    }
    fn wants() -> &'static [ProductId] {
        &[ProductId::Lakes, ProductId::Rivers, ProductId::PrecipitationWarmShare]
    }
    fn applies(body: &BodyContext) -> Result<(), Precondition> {
        crate::climate::Climate::applies(body)
    }

    fn run(ctx: &RunContext<'_>, inputs: &Products, p: &BiomesParams) -> Result<Products, Error> {
        let mut out = Products::new();
        if inputs.biome().is_some() {
            return Ok(out);
        }
        let grid = ctx.grid;
        let n = grid.len();
        let elev = &inputs.elevation().ok_or(Error::MissingInput(ProductId::Elevation))?.values;
        let t = &inputs.temperature().ok_or(Error::MissingInput(ProductId::Temperature))?.values;
        let tr = &inputs.temperature_range().ok_or(Error::MissingInput(ProductId::TemperatureRange))?.values;
        let pr = &inputs.precipitation().ok_or(Error::MissingInput(ProductId::Precipitation))?.values;
        let lakes = inputs.lakes();
        let share = inputs.precipitation_warm_share().map(|f| &f.values);
        let mut wet = vec![false; n];
        if let Some(rv) = inputs.rivers() {
            for node in rv.nodes.rows.values() {
                let c = node.cell as usize;
                if node.discharge_m3_s >= p.wetland_discharge_m3_s && elev[c] >= 0.0 && elev[c] < 300.0 {
                    let s = crate::elevation::slope(grid, elev, node.cell, ctx.radius_km());
                    if s < p.wetland_max_slope {
                        wet[c] = true;
                        // Floodplains are wider than the channel.
                        for nb in grid.neighbors(node.cell).iter() {
                            if elev[nb as usize] >= 0.0 {
                                wet[nb as usize] = true;
                            }
                        }
                    }
                }
            }
        }
        let values: Vec<Biome> = (0..n)
            .map(|i| {
                if lakes.is_some_and(|l| l.owner.values[i] != 0) {
                    return Biome::Lake;
                }
                if elev[i] < 0.0 {
                    return Biome::Ocean;
                }
                let lat = grid.center(i as Cell).to_lat_lon().lat_deg;
                let warm_share = share.map_or(0.5, |s| f64::from(s[i]));
                let c = CellClimate { t_mean: f64::from(t[i]), t_range: f64::from(tr[i]), precip_mm: f64::from(pr[i]), elevation_m: f64::from(elev[i]), lat_deg: lat, warm_share };
                let b = classify(&c, p);
                if wet[i] && !matches!(b, Biome::IceSheet | Biome::HotDesert | Biome::ColdDesert | Biome::Alpine) {
                    // Wet floodplain, unless the climate leaves it bone dry.
                    if f64::from(pr[i]) > 0.4 * pet_mm(c.t_mean) {
                        return Biome::Wetland;
                    }
                }
                b
            })
            .collect();
        out.insert(ProductData::Biome(Field::from_values(grid.spec(), values)));
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn at(t_mean: f64, t_range: f64, precip_mm: f64, lat: f64) -> Biome {
        at_s(t_mean, t_range, precip_mm, lat, 0.5)
    }

    fn at_s(t_mean: f64, t_range: f64, precip_mm: f64, lat: f64, warm_share: f64) -> Biome {
        classify(&CellClimate { t_mean, t_range, precip_mm, elevation_m: 200.0, lat_deg: lat, warm_share }, &BiomesParams::default())
    }

    #[test]
    fn earth_stations() {
        assert_eq!(at(26.5, 2.0, 2300.0, -3.0), Biome::TropicalRainforest, "Manaus");
        assert_eq!(at(30.0, 18.0, 25.0, 25.0), Biome::HotDesert, "Sahara");
        assert_eq!(at(-10.0, 62.0, 240.0, 62.0), Biome::BorealForest, "Yakutsk (dry boreal)");
        assert_eq!(at(-12.0, 33.0, 200.0, 72.0), Biome::Tundra, "Taymyr");
        assert_eq!(at(-30.0, 30.0, 100.0, 75.0), Biome::IceSheet, "Greenland");
        assert_eq!(at_s(17.5, 15.0, 600.0, 38.0, 0.2), Biome::Mediterranean, "Athens-like");
        assert_eq!(at_s(15.0, 10.0, 450.0, 37.0, 0.15), Biome::Mediterranean, "San Francisco-like");
        assert_ne!(at_s(17.5, 15.0, 600.0, 38.0, 0.6), Biome::Mediterranean, "summer rain is not Mediterranean");
        assert_eq!(at(10.5, 12.0, 2200.0, 47.5), Biome::TemperateRainforest, "Olympic Peninsula");
        assert_eq!(at(11.0, 27.0, 900.0, 41.0), Biome::TemperateForest, "Ohio");
        assert_eq!(at(9.0, 30.0, 380.0, 47.0), Biome::TemperateGrassland, "steppe");
        assert_eq!(at(26.0, 6.0, 900.0, 10.0), Biome::Savanna, "Sahel/savanna");
        assert_eq!(at(9.0, 18.0, 150.0, 40.0), Biome::ColdDesert, "Great Basin");
    }
}
