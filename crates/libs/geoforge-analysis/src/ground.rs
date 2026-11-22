//! What the rule-based views (resources, soils, hazards) read from a generated world, per
//! cell, gathered once from the stored products.

use crate::koppen::{self, Koppen, KoppenRules};
use geoforge_core::grid::{Cell, Grid};
use geoforge_core::product::{ProductId, Products};
use geoforge_core::surface::GeologicProvince;

/// A world's ground, cell by cell.
pub struct Ground<'a> {
    pub grid: &'a dyn Grid,
    pub radius_km: f64,
    /// Below sea level on a world with seas.
    pub sea: Vec<bool>,
    pub elevation: Vec<f64>,
    pub lat: Vec<f64>,
    pub province: Vec<Option<GeologicProvince>>,
    /// Rock slug, "" where none.
    pub rock: Vec<&'a str>,
    /// Rock age, Ma.
    pub rock_age: &'a [f32],
    pub temperature: &'a [f32],
    pub temperature_range: &'a [f32],
    pub precipitation: &'a [f32],
    pub koppen: Vec<Option<Koppen>>,
    /// Highest Strahler order of a river through the cell, 0 for none.
    pub river_order: Vec<u8>,
    /// A closed (salt) lake.
    pub salt_lake: Vec<bool>,
    /// In a basin with no outlet to the sea.
    pub closed_basin: Vec<bool>,
    /// Biome slug, "" where none.
    pub biome: Vec<&'a str>,
    /// Steepest drop to a neighbour, m per km.
    pub slope: Vec<f64>,
}

/// The products a [`Ground`] reads.
pub const READS: [ProductId; 12] = [
    ProductId::Provinces,
    ProductId::Lithology,
    ProductId::RockAge,
    ProductId::Elevation,
    ProductId::Temperature,
    ProductId::TemperatureRange,
    ProductId::Precipitation,
    ProductId::PrecipitationWarmShare,
    ProductId::Rivers,
    ProductId::Lakes,
    ProductId::Watersheds,
    ProductId::Biome,
];

impl<'a> Ground<'a> {
    /// `None` until every product in [`READS`] is stored.
    #[must_use]
    pub fn read(p: &'a Products, grid: &'a dyn Grid, radius_km: f64, sea_world: bool) -> Option<Self> {
        let elev: Vec<f64> = p.elevation()?.values.iter().map(|&v| f64::from(v)).collect();
        let n = grid.len();
        let provinces = p.provinces()?;
        let lithology = p.lithology()?;
        let biome = p.biome()?;
        let lakes = p.lakes()?;
        let rivers = p.rivers()?;
        let watersheds = p.watersheds()?;
        let mut river_order = vec![0u8; n];
        for e in &rivers.edges {
            for k in [e.from, e.to] {
                if let Some(node) = rivers.nodes.get(k) {
                    let o = &mut river_order[node.cell as usize];
                    *o = (*o).max(e.attrs.strahler_order);
                }
            }
        }
        let spacing_km = grid.spacing_rad() * radius_km;
        let slope = (0..n).map(|c| grid.neighbors(c as Cell).iter().map(|nb| (elev[c] - elev[nb as usize]).max(0.0)).fold(0.0, f64::max) / spacing_km).collect();
        Some(Self {
            grid,
            radius_km,
            sea: elev.iter().map(|&e| sea_world && e < 0.0).collect(),
            lat: (0..n).map(|c| grid.center(c as Cell).to_lat_lon().lat_deg).collect(),
            province: (0..n).map(|c| provinces.row_at(c as Cell).map(|r| r.kind)).collect(),
            rock: (0..n).map(|c| lithology.row_at(c as Cell).map_or("", |r| r.slug.as_str())).collect(),
            rock_age: &p.rock_age()?.values,
            temperature: &p.temperature()?.values,
            temperature_range: &p.temperature_range()?.values,
            precipitation: &p.precipitation()?.values,
            koppen: koppen::map(p, sea_world, &KoppenRules::default())?,
            river_order,
            salt_lake: (0..n).map(|c| lakes.row_at(c as Cell).is_some_and(|l| l.endorheic)).collect(),
            closed_basin: (0..n).map(|c| watersheds.row_at(c as Cell).is_some_and(|w| !w.reaches_sea)).collect(),
            biome: (0..n).map(|c| biome.row_at(c as Cell).map_or("", |r| r.slug.as_str())).collect(),
            elevation: elev,
            slope,
        })
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.sea.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.sea.is_empty()
    }

    /// The cell's Köppen group letter ('B' for arid), or ' ' at sea.
    #[must_use]
    pub fn group(&self, c: usize) -> char {
        self.koppen[c].map_or(' ', |k| k.code().as_bytes()[0] as char)
    }

    #[must_use]
    pub fn code(&self, c: usize) -> &'static str {
        self.koppen[c].map_or("", Koppen::code)
    }

    #[must_use]
    pub fn province_is(&self, c: usize, kinds: &[GeologicProvince]) -> bool {
        self.province[c].is_some_and(|p| kinds.contains(&p))
    }

    /// Coldest month, °C.
    #[must_use]
    pub fn coldest(&self, c: usize) -> f64 {
        f64::from(self.temperature[c]) - 0.5 * f64::from(self.temperature_range[c])
    }

    /// Warmest month, °C.
    #[must_use]
    pub fn warmest(&self, c: usize) -> f64 {
        f64::from(self.temperature[c]) + 0.5 * f64::from(self.temperature_range[c])
    }

    /// The cell or a neighbour has rock `slug`.
    #[must_use]
    pub fn touches_rock(&self, c: usize, slug: &str) -> bool {
        self.rock[c] == slug || self.grid.neighbors(c as Cell).iter().any(|n| self.rock[n as usize] == slug)
    }

    /// Land next to the sea.
    #[must_use]
    pub fn coast(&self, c: usize) -> bool {
        !self.sea[c] && self.grid.neighbors(c as Cell).iter().any(|n| self.sea[n as usize])
    }

    /// Cell area, km².
    #[must_use]
    pub fn area_km2(&self, c: usize) -> f64 {
        self.grid.area_sr(c as Cell) * self.radius_km * self.radius_km
    }

    /// Geodesic distance from `sources` over land, km, up to `max_km` (beyond: infinity).
    #[must_use]
    pub fn distance_over_land_km(&self, sources: impl IntoIterator<Item = Cell>, max_km: f64) -> Vec<f64> {
        geoforge_core::grid::distance_field(self.grid, sources, max_km / self.radius_km, |c| !self.sea[c as usize]).into_iter().map(|d| d * self.radius_km).collect()
    }

    /// Geodesic distance from `sources` over land and sea, km, up to `max_km`.
    #[must_use]
    pub fn distance_over_land_or_sea_km(&self, sources: impl IntoIterator<Item = Cell>, max_km: f64) -> Vec<f64> {
        geoforge_core::grid::distance_field(self.grid, sources, max_km / self.radius_km, |_| true).into_iter().map(|d| d * self.radius_km).collect()
    }
}
