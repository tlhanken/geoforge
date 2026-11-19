//! Biomes layer: a vegetation/land-cover class per cell from climate, with lakes and
//! wetlands from hydrology.
//!
//! The classification follows Köppen's thresholds where they are physical (the 10 °C
//! warmest-month tree line, the 18 °C coldest-month tropics, the aridity threshold that
//! scales with temperature) and Whittaker's precipitation steps within each thermal zone.
//! Monthly extremes are taken as mean ± range/2.
//!
//! Biomes come from a [`Catalog`]: Earth's sixteen are built in and classified as above;
//! `[params.biomes] extra_biomes` adds more (fantasy or otherwise) as climate envelopes,
//! checked first, in order. The biome partition's table is the catalog, so a painted edit can
//! use any of its biomes.

use crate::climate::pet_mm;
use crate::noise::{uniform_ranks, TerrainNoise};
use geoforge_core::grid::Cell;
use geoforge_core::key::Key;
use geoforge_core::product::{Generator, LayerId, Precondition, ProductData, ProductId, Products, RunContext};
use geoforge_core::shapes::{Field, Partition, Table};
use geoforge_core::surface::BiomeClass;
use geoforge_core::{BodyContext, Error};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct BiomesParams {
    /// Elevation above which cold climates count as alpine rather than polar, m.
    pub alpine_min_m: f64,
    /// River discharge above which flat lowland floodplains become wetland, m³/s.
    pub wetland_discharge_m3_s: f64,
    /// Steepest slope (rise/run) still counted as a floodplain.
    pub wetland_max_slope: f64,
    /// Biomes added to Earth's (see [`BiomeDef`]).
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub extra_biomes: Vec<BiomeDef>,
}

impl Default for BiomesParams {
    fn default() -> Self {
        Self { alpine_min_m: 1800.0, wetland_discharge_m3_s: 2500.0, wetland_max_slope: 4e-4, extra_biomes: Vec::new() }
    }
}

/// Earth's biomes. Keys are the variant number + 1, so stored biomes stay stable.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Earth {
    Ocean,
    Lake,
    IceSheet,
    Tundra,
    BorealForest,
    TemperateRainforest,
    TemperateForest,
    TemperateGrassland,
    Mediterranean,
    ColdDesert,
    HotDesert,
    Savanna,
    TropicalSeasonalForest,
    TropicalRainforest,
    Alpine,
    Wetland,
}

impl Earth {
    pub const ALL: [Self; 16] = [
        Self::Ocean,
        Self::Lake,
        Self::IceSheet,
        Self::Tundra,
        Self::BorealForest,
        Self::TemperateRainforest,
        Self::TemperateForest,
        Self::TemperateGrassland,
        Self::Mediterranean,
        Self::ColdDesert,
        Self::HotDesert,
        Self::Savanna,
        Self::TropicalSeasonalForest,
        Self::TropicalRainforest,
        Self::Alpine,
        Self::Wetland,
    ];

    #[must_use]
    pub fn key(self) -> Key {
        Key(self as u64 + 1)
    }

    #[must_use]
    pub fn class(self) -> BiomeClass {
        let (slug, name, color) = match self {
            Self::Ocean => ("ocean", "Ocean", [28, 62, 120]),
            Self::Lake => ("lake", "Lake", [60, 110, 190]),
            Self::IceSheet => ("ice_sheet", "Ice Sheet", [240, 248, 255]),
            Self::Tundra => ("tundra", "Tundra", [176, 184, 160]),
            Self::BorealForest => ("boreal_forest", "Boreal Forest", [58, 98, 72]),
            Self::TemperateRainforest => ("temperate_rainforest", "Temperate Rainforest", [38, 120, 80]),
            Self::TemperateForest => ("temperate_deciduous_forest", "Temperate Deciduous Forest", [92, 148, 70]),
            Self::TemperateGrassland => ("temperate_grassland", "Temperate Grassland", [186, 196, 110]),
            Self::Mediterranean => ("mediterranean_woodland", "Mediterranean Woodland", [160, 150, 70]),
            Self::ColdDesert => ("cold_desert", "Cold Desert", [196, 182, 150]),
            Self::HotDesert => ("hot_desert", "Hot Desert", [234, 208, 140]),
            Self::Savanna => ("savanna", "Savanna", [200, 180, 80]),
            Self::TropicalSeasonalForest => ("tropical_seasonal_forest", "Tropical Seasonal Forest", [120, 160, 50]),
            Self::TropicalRainforest => ("tropical_rainforest", "Tropical Rainforest", [20, 110, 30]),
            Self::Alpine => ("alpine", "Alpine", [150, 140, 130]),
            Self::Wetland => ("wetland", "Wetland", [70, 130, 120]),
        };
        BiomeClass { slug: slug.into(), name: name.into(), color, tags: Vec::new() }
    }

    #[must_use]
    pub fn from_slug(slug: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|e| e.class().slug == slug)
    }
}

/// Where an extra biome grows. Every bound is optional; ranges are `[min, max]`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Envelope {
    /// Mean annual temperature, °C.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub temperature_c: Option<[f64; 2]>,
    /// Annual precipitation, mm.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub precipitation_mm: Option<[f64; 2]>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub elevation_m: Option<[f64; 2]>,
    /// Earth biomes (slugs) it may replace; empty = any land.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub within: Vec<String>,
    /// Share of the matching ground it takes, in map-scale patches (1 = all of it).
    #[serde(default = "full_cover")]
    pub cover: f64,
}

fn full_cover() -> f64 {
    1.0
}

/// A biome a world adds to Earth's.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BiomeDef {
    #[serde(flatten)]
    pub class: BiomeClass,
    #[serde(flatten)]
    pub envelope: Envelope,
}

/// Keys of extra biomes start here, clear of the built-ins.
const EXTRA_KEY_BASE: u64 = 101;

struct Extra {
    key: Key,
    def: BiomeDef,
    /// Earth biomes it may replace; empty = any land.
    within: Vec<Earth>,
}

impl Extra {
    fn admits(&self, base: Earth, c: &CellClimate) -> bool {
        let inside = |r: Option<[f64; 2]>, v: f64| r.is_none_or(|[lo, hi]| (lo..=hi).contains(&v));
        let e = &self.def.envelope;
        let ground = if self.within.is_empty() { !matches!(base, Earth::Ocean | Earth::Lake) } else { self.within.contains(&base) };
        ground && inside(e.temperature_c, c.t_mean) && inside(e.precipitation_mm, c.precip_mm) && inside(e.elevation_m, c.elevation_m)
    }
}

/// Earth's biomes plus a world's extras, keyed as the biome partition stores them.
pub struct Catalog {
    extras: Vec<Extra>,
}

impl Catalog {
    /// # Errors
    /// If a slug repeats or `within` names a biome Earth does not have.
    pub fn new(extra: &[BiomeDef]) -> Result<Self, Error> {
        let mut slugs: BTreeSet<String> = Earth::ALL.iter().map(|e| e.class().slug).collect();
        let mut extras = Vec::new();
        for (i, def) in extra.iter().enumerate() {
            if !slugs.insert(def.class.slug.clone()) {
                return Err(Error::Invalid(format!("biome '{}' is defined twice", def.class.slug)));
            }
            let within = def
                .envelope
                .within
                .iter()
                .map(|s| Earth::from_slug(s).ok_or_else(|| Error::Invalid(format!("biome '{}': within names unknown biome '{s}'", def.class.slug))))
                .collect::<Result<_, _>>()?;
            extras.push(Extra { key: Key(EXTRA_KEY_BASE + i as u64), def: def.clone(), within });
        }
        Ok(Self { extras })
    }

    /// The biome partition's table.
    #[must_use]
    pub fn table(&self) -> Table<BiomeClass> {
        let mut t = Table::new();
        for e in Earth::ALL {
            t.insert(e.key(), e.class());
        }
        for x in &self.extras {
            t.insert(x.key, x.def.class.clone());
        }
        t
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

/// Classify a land cell among Earth's biomes.
#[must_use]
pub fn classify(c: &CellClimate, p: &BiomesParams) -> Earth {
    let warm = c.t_mean + c.t_range / 2.0;
    let cold = c.t_mean - c.t_range / 2.0;
    // Alpine: cold because of altitude. At sea level the same place would be warm enough
    // for trees, so it is a mountain climate rather than a polar one.
    let warm_at_sea_level = c.t_mean + c.t_range / 2.0 + 6.5 * c.elevation_m.max(0.0) / 1000.0;
    let high = c.elevation_m >= p.alpine_min_m && warm_at_sea_level >= 10.0;
    if warm < 0.0 {
        return if high { Earth::Alpine } else { Earth::IceSheet };
    }
    if warm < 10.0 {
        return if high { Earth::Alpine } else { Earth::Tundra };
    }
    let th = aridity_threshold(c.t_mean, c.warm_share);
    let pr = c.precip_mm;
    if pr < 0.5 * th {
        return if c.t_mean >= 18.0 { Earth::HotDesert } else { Earth::ColdDesert };
    }
    let semi_arid = pr < th;
    if cold >= 18.0 {
        // Tropical.
        return if semi_arid {
            Earth::Savanna
        } else if pr >= 1800.0 && c.t_range < 6.0 || pr >= 2300.0 {
            Earth::TropicalRainforest
        } else if pr >= 1300.0 {
            Earth::TropicalSeasonalForest
        } else {
            Earth::Savanna
        };
    }
    if cold < -3.0 && warm < 22.0 {
        // Cold-winter, cool-summer: boreal, or steppe where too dry for trees.
        return if semi_arid { Earth::TemperateGrassland } else { Earth::BorealForest };
    }
    // Temperate (and hot-summer continental).
    if semi_arid {
        return if c.t_mean >= 18.0 { Earth::Savanna } else { Earth::TemperateGrassland };
    }
    // Dry summers, wet mild winters (Köppen Cs).
    if c.warm_share < 0.35 && cold > -3.0 && pr < 1200.0 {
        return Earth::Mediterranean;
    }
    if pr >= 2000.0 || (pr >= 1500.0 && c.t_range < 14.0 && c.t_mean < 16.0) {
        return Earth::TemperateRainforest;
    }
    if pr < 1.6 * th {
        return Earth::TemperateGrassland;
    }
    if c.t_mean >= 20.0 && cold >= 10.0 {
        return if pr >= 1300.0 { Earth::TropicalSeasonalForest } else { Earth::Savanna };
    }
    Earth::TemperateForest
}

pub struct Biomes;

impl Generator for Biomes {
    type Params = BiomesParams;
    const NAME: &'static str = "biomes";
    const VERSION: u32 = 2;
    const LAYER: LayerId = LayerId::Biomes;

    fn needs() -> &'static [ProductId] {
        &[ProductId::Elevation, ProductId::Temperature, ProductId::TemperatureRange, ProductId::Precipitation]
    }
    fn wants() -> &'static [ProductId] {
        &[ProductId::Lakes, ProductId::Rivers, ProductId::PrecipitationWarmShare]
    }
    fn applies(body: &BodyContext) -> Result<(), Precondition> {
        crate::preconditions::surface_water(body)
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
        let catalog = Catalog::new(&p.extra_biomes)?;
        // Patches for extras that take only part of their ground: one smooth uniform field each.
        let patches: Vec<Option<Vec<f64>>> = catalog
            .extras
            .iter()
            .map(|x| {
                (x.def.envelope.cover < 1.0).then(|| {
                    let noise = TerrainNoise::new(ctx.seed(ProductId::Biome).child(x.def.class.slug.as_bytes()).value());
                    uniform_ranks(&(0..n).into_par_iter().map(|i| noise.fbm(grid.center(i as Cell), 6.0, 3)).collect::<Vec<_>>())
                })
            })
            .collect();
        let owner: Vec<u32> = (0..n)
            .map(|i| {
                let lat = grid.center(i as Cell).to_lat_lon().lat_deg;
                let warm_share = share.map_or(0.5, |s| f64::from(s[i]));
                let c = CellClimate { t_mean: f64::from(t[i]), t_range: f64::from(tr[i]), precip_mm: f64::from(pr[i]), elevation_m: f64::from(elev[i]), lat_deg: lat, warm_share };
                let base = if lakes.is_some_and(|l| l.owner.values[i] != 0) {
                    Earth::Lake
                } else if elev[i] < 0.0 {
                    Earth::Ocean
                } else {
                    let b = classify(&c, p);
                    // Wet floodplain, unless the climate leaves it bone dry.
                    let floods = wet[i] && !matches!(b, Earth::IceSheet | Earth::HotDesert | Earth::ColdDesert | Earth::Alpine);
                    if floods && f64::from(pr[i]) > 0.4 * pet_mm(c.t_mean) { Earth::Wetland } else { b }
                };
                let extra = catalog.extras.iter().zip(&patches).find(|(x, patch)| x.admits(base, &c) && patch.as_ref().is_none_or(|u| u[i] < x.def.envelope.cover));
                extra.map_or(base.key(), |(x, _)| x.key).0 as u32
            })
            .collect();
        out.insert(ProductData::Biome(Partition { owner: Field::from_values(grid.spec(), owner), table: catalog.table() }));
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn at(t_mean: f64, t_range: f64, precip_mm: f64, lat: f64) -> Earth {
        at_s(t_mean, t_range, precip_mm, lat, 0.5)
    }

    fn at_s(t_mean: f64, t_range: f64, precip_mm: f64, lat: f64, warm_share: f64) -> Earth {
        classify(&CellClimate { t_mean, t_range, precip_mm, elevation_m: 200.0, lat_deg: lat, warm_share }, &BiomesParams::default())
    }

    #[test]
    fn earth_stations() {
        assert_eq!(at(26.5, 2.0, 2300.0, -3.0), Earth::TropicalRainforest, "Manaus");
        assert_eq!(at(30.0, 18.0, 25.0, 25.0), Earth::HotDesert, "Sahara");
        assert_eq!(at(-10.0, 62.0, 240.0, 62.0), Earth::BorealForest, "Yakutsk (dry boreal)");
        assert_eq!(at(-12.0, 33.0, 200.0, 72.0), Earth::Tundra, "Taymyr");
        assert_eq!(at(-30.0, 30.0, 100.0, 75.0), Earth::IceSheet, "Greenland");
        assert_eq!(at_s(17.5, 15.0, 600.0, 38.0, 0.2), Earth::Mediterranean, "Athens-like");
        assert_eq!(at_s(15.0, 10.0, 450.0, 37.0, 0.15), Earth::Mediterranean, "San Francisco-like");
        assert_ne!(at_s(17.5, 15.0, 600.0, 38.0, 0.6), Earth::Mediterranean, "summer rain is not Mediterranean");
        assert_eq!(at(10.5, 12.0, 2200.0, 47.5), Earth::TemperateRainforest, "Olympic Peninsula");
        assert_eq!(at(11.0, 27.0, 900.0, 41.0), Earth::TemperateForest, "Ohio");
        assert_eq!(at(9.0, 30.0, 380.0, 47.0), Earth::TemperateGrassland, "steppe");
        assert_eq!(at(26.0, 6.0, 900.0, 10.0), Earth::Savanna, "Sahel/savanna");
        assert_eq!(at(9.0, 18.0, 150.0, 40.0), Earth::ColdDesert, "Great Basin");
    }

    fn glowcap(cover: f64, within: &[&str]) -> BiomeDef {
        BiomeDef {
            class: BiomeClass { slug: "glowcap_forest".into(), name: "Glowcap Forest".into(), color: [90, 58, 138], tags: vec!["magical".into()] },
            envelope: Envelope { temperature_c: Some([0.0, 15.0]), precipitation_mm: None, elevation_m: None, within: within.iter().map(|s| (*s).to_string()).collect(), cover },
        }
    }

    #[test]
    fn extras_extend_the_catalog() {
        let cat = Catalog::new(&[glowcap(0.3, &["boreal_forest", "temperate_rainforest"])]).unwrap();
        let t = cat.table();
        assert_eq!(t.rows.len(), 17);
        assert_eq!(t.rows[&Key(EXTRA_KEY_BASE)].slug, "glowcap_forest");
        assert_eq!(t.rows[&Earth::Ocean.key()].slug, "ocean");
        let x = &cat.extras[0];
        let cool = CellClimate { t_mean: 5.0, t_range: 20.0, precip_mm: 800.0, elevation_m: 300.0, lat_deg: 55.0, warm_share: 0.6 };
        assert!(x.admits(Earth::BorealForest, &cool));
        assert!(!x.admits(Earth::Tundra, &cool), "only within the listed biomes");
        assert!(!x.admits(Earth::BorealForest, &CellClimate { t_mean: 20.0, ..cool }), "outside the temperature range");
        assert!(Catalog::new(&[glowcap(1.0, &["nowhere"])]).is_err());
        assert!(Catalog::new(&[glowcap(1.0, &[]), glowcap(1.0, &[])]).is_err());
        // Any land when `within` is empty, never sea.
        let any = Catalog::new(&[glowcap(1.0, &[])]).unwrap();
        assert!(any.extras[0].admits(Earth::Tundra, &cool) && !any.extras[0].admits(Earth::Ocean, &cool));
    }

    #[test]
    fn extras_parse_from_toml() {
        let d: BiomeDef = toml::from_str(r##"
            slug = "ashlands"
            name = "Ashlands"
            color = "#4a4040"
            within = ["hot_desert", "savanna"]
            temperature_c = [20.0, 40.0]
            cover = 0.25
        "##)
        .unwrap();
        assert_eq!(d.envelope.within.len(), 2);
        assert_eq!(d.envelope.temperature_c, Some([20.0, 40.0]));
        assert!((d.envelope.cover - 0.25).abs() < 1e-12);
    }
}
