//! Derived views: maps that are pure functions of stored products (and the body). They are
//! worked out when rendered or exported and never stored, so they cannot go stale and cost
//! nothing to keep.
//!
//! Three kinds: **class maps** (Köppen, soils, lightning tiers), exported as exact-colour PNGs
//! with a legend; **chance maps** (each deposit kind, each hazard, fertility, timber),
//! exported as 8-bit grayscale PNGs, 0–255 for 0–1, with a note of what 1 means; and
//! **quantity maps** (sea-level pressure), exported as 16-bit grayscale PNGs with the range
//! and units.

use geoforge_core::grid::{AnyGrid, Cell, Grid};
use geoforge_core::magnetism::MagneticField;
use geoforge_core::product::{ProductId, Products};
use geoforge_core::BodyContext;
use geoforge_surface::ground::{self, Ground};
use geoforge_surface::hazards::{self, Hazard, Lightning, Measure};
use geoforge_surface::koppen::{self, Koppen, KoppenRules};
use geoforge_surface::oceans::{Reef, SeaIce};
use geoforge_surface::resources::{self, Deposit, Soil};
use serde::Serialize;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum View {
    /// Köppen–Geiger climate classes.
    Koppen,
    /// Soil classes.
    Soil,
    /// Soil fertility for farming.
    Fertility,
    /// Standing timber.
    Timber,
    /// Chance a 10,000 km² region holds a deposit of this kind.
    Deposit(Deposit),
    /// Chance of this hazard in a 10,000 km² region in a year.
    Hazard(Hazard),
    /// Lightning tiers.
    Lightning,
    /// Sea-level pressure in January (`false`) or July (`true`).
    Pressure { july: bool },
    /// Sea ice: none, seasonal or perennial.
    SeaIce,
    /// Coral reefs and atolls.
    Reefs,
    /// The sea's productivity.
    Productivity,
}

impl View {
    /// Every view.
    #[must_use]
    pub fn all() -> Vec<Self> {
        let mut v = vec![Self::Koppen, Self::Soil, Self::Fertility, Self::Timber];
        v.extend(Deposit::ALL.map(Self::Deposit));
        v.extend(Hazard::ALL.map(Self::Hazard));
        v.push(Self::Lightning);
        v.extend([Self::Pressure { july: false }, Self::Pressure { july: true }]);
        v.extend([Self::SeaIce, Self::Reefs, Self::Productivity]);
        v
    }

    /// Views by group name: `climate_maps` (Köppen, pressure), `oceans` (sea ice, reefs,
    /// productivity), `resources` (deposits, soils, fertility, timber) or `hazards`.
    #[must_use]
    pub fn group(name: &str) -> Option<Vec<Self>> {
        match name {
            "climate_maps" => Some(Self::all().into_iter().filter(|v| matches!(v, Self::Koppen | Self::Pressure { .. })).collect()),
            "oceans" => Some(vec![Self::SeaIce, Self::Reefs, Self::Productivity]),
            "resources" => Some(Self::all().into_iter().filter(|v| matches!(v, Self::Soil | Self::Fertility | Self::Timber | Self::Deposit(_))).collect()),
            "hazards" => Some(Self::all().into_iter().filter(|v| matches!(v, Self::Hazard(_) | Self::Lightning)).collect()),
            _ => None,
        }
    }

    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Self::Koppen => "koppen",
            Self::Soil => "soil",
            Self::Fertility => "fertility",
            Self::Timber => "timber",
            Self::Deposit(d) => d.slug(),
            Self::Hazard(h) => h.slug(),
            Self::Lightning => "lightning",
            Self::Pressure { july: false } => "pressure_jan",
            Self::Pressure { july: true } => "pressure_jul",
            Self::SeaIce => "sea_ice",
            Self::Reefs => "reefs",
            Self::Productivity => "productivity",
        }
    }

    #[must_use]
    pub fn from_name(name: &str) -> Option<Self> {
        Self::all().into_iter().find(|v| v.name() == name)
    }

    /// The products the view is worked out from.
    #[must_use]
    pub fn reads(self) -> &'static [ProductId] {
        match self {
            Self::Koppen => &[ProductId::Temperature, ProductId::TemperatureRange, ProductId::Precipitation, ProductId::PrecipitationWarmShare],
            Self::Pressure { july: false } => &[ProductId::WindJan, ProductId::Temperature, ProductId::TemperatureRange, ProductId::Elevation],
            Self::Pressure { july: true } => &[ProductId::WindJul, ProductId::Temperature, ProductId::TemperatureRange, ProductId::Elevation],
            Self::SeaIce | Self::Reefs | Self::Productivity => &[ProductId::Elevation, ProductId::Temperature, ProductId::TemperatureRange, ProductId::SstAnomaly],
            _ => &ground::READS,
        }
    }
}

/// One class of a class map, with its share of the mapped area (the land, for land maps).
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Class {
    pub code: String,
    pub name: String,
    pub color: [u8; 3],
    pub area_share: f64,
}

/// A class map: each cell's class (an index into `legend`), `None` at sea.
#[derive(Debug, Clone, PartialEq)]
pub struct Classes {
    pub cells: Vec<Option<usize>>,
    pub legend: Vec<Class>,
}

impl Classes {
    /// Each cell's colour, `None` at sea.
    #[must_use]
    pub fn colors(&self) -> Vec<Option<[u8; 3]>> {
        self.cells.iter().map(|c| c.map(|k| self.legend[k].color)).collect()
    }
}

/// A chance map: each cell's value, 0–1, `None` where it does not apply (sea, mostly).
#[derive(Debug, Clone, PartialEq)]
pub struct Chances {
    pub cells: Vec<Option<f64>>,
    /// What a value of 1 means, and how the value is made.
    pub meaning: String,
}

/// A quantity map: each cell's value in `units`.
#[derive(Debug, Clone, PartialEq)]
pub struct Quantity {
    pub cells: Vec<f64>,
    pub units: &'static str,
    /// Contour spacing for previews, in `units`.
    pub contour: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Map {
    Classes(Classes),
    Chances(Chances),
    Quantity(Quantity),
}

/// What views are worked out from; the per-cell ground is read once, on first need.
pub struct Source<'a> {
    pub products: &'a Products,
    pub grid: &'a AnyGrid,
    pub body: &'a BodyContext,
    ground: std::cell::OnceCell<Option<Ground<'a>>>,
}

impl<'a> Source<'a> {
    #[must_use]
    pub fn new(products: &'a Products, grid: &'a AnyGrid, body: &'a BodyContext) -> Self {
        Self { products, grid, body, ground: std::cell::OnceCell::new() }
    }

    fn sea_world(&self) -> bool {
        self.body.interior.has_surface_liquid_water()
    }

    fn ground(&self) -> Option<&Ground<'a>> {
        self.ground.get_or_init(|| Ground::read(self.products, self.grid, self.body.exterior.radius_km, self.sea_world())).as_ref()
    }

    /// Cells below sea level on a world with seas.
    #[must_use]
    pub fn sea(&self) -> Vec<bool> {
        let wet = self.sea_world();
        self.products.elevation().or(self.products.rough_elevation()).map_or_else(|| vec![false; self.grid.len()], |e| e.values.iter().map(|&v| wet && v < 0.0).collect())
    }
}

/// Work out a view; `None` while the products it reads are missing.
#[must_use]
pub fn compute(view: View, s: &Source) -> Option<Map> {
    let sea = s.sea_world();
    let chances = |cells: Vec<Option<f64>>, meaning: String| Map::Chances(Chances { cells, meaning });
    Some(match view {
        View::Koppen => Map::Classes(tally(s.grid, &koppen::map(s.products, sea, &KoppenRules::default())?, |k| (k.code(), k.name(), k.color()), &Koppen::ALL)),
        View::Soil => {
            let g = s.ground()?;
            let soils: Vec<Option<Soil>> = (0..g.len()).map(|c| resources::soil(g, c)).collect();
            Map::Classes(tally(s.grid, &soils, |k| (k.code(), k.name(), k.color()), &Soil::ALL))
        }
        View::Lightning => Map::Classes(tally(s.grid, &hazards::lightning(s.ground()?), |k| k.info(), &Lightning::ALL)),
        View::Fertility => {
            let g = s.ground()?;
            chances((0..g.len()).map(|c| resources::soil(g, c).map(Soil::fertility)).collect(), "soil fertility for farming, 0 (none) to 1 (floodplain, black earth)".into())
        }
        View::Timber => {
            let g = s.ground()?;
            chances((0..g.len()).map(|c| resources::timber(g, c)).collect(), "standing timber, 0 (none) to 1 (rainforest), below the tree line".into())
        }
        View::Deposit(d) => {
            let g = s.ground()?;
            let lambda = resources::prospectivity(g, d);
            let name = d.slug().replace('_', "–");
            chances((0..g.len()).map(|c| (!g.sea[c]).then(|| 1.0 - (-lambda[c]).exp())).collect(), format!("chance a 10,000 km² region holds at least one {name} deposit, 1 − e^(−λ)"))
        }
        View::Hazard(h) => {
            let g = s.ground()?;
            let lambda = hazards::rate(g, &hazards::Setting { products: s.products, aurora: aurora_of(s.body) }, h);
            let m = h.measure();
            let meaning = match m {
                Measure::Window(1.0) => format!("chance of {} in a 10,000 km² region in a year, 1 − e^(−λ)", h.event()),
                Measure::Window(years) => format!("chance of {} in a 10,000 km² region in {years} years, 1 − e^(−λT)", h.event()),
                Measure::Days => format!("share of the year's days with {}", h.event()),
            };
            chances((0..g.len()).map(|c| (h.at_sea() || !g.sea[c]).then(|| m.value(lambda[c]))).collect(), meaning)
        }
        View::Pressure { july } => {
            let wind = if july { s.products.wind_jul() } else { s.products.wind_jan() }?;
            let (t, range) = (&s.products.temperature()?.values, &s.products.temperature_range()?.values);
            // July: the north at its warmest month, the south at its coldest.
            let season_t: Vec<f64> = (0..s.grid.len())
                .map(|c| {
                    let north = s.grid.center(c as Cell).to_lat_lon().lat_deg >= 0.0;
                    f64::from(t[c]) + 0.5 * f64::from(range[c]) * if north == july { 1.0 } else { -1.0 }
                })
                .collect();
            use geoforge_surface::pressure::{sea_level_pressure, Air, Planet};
            let land: Vec<bool> = s.sea().iter().map(|&wet| !wet).collect();
            let elevation = &s.products.elevation().or(s.products.rough_elevation())?.values;
            let b = s.body;
            let planet = Planet { radius_km: b.exterior.radius_km, rotation_h: b.exterior.rotation_h, mean_hpa: b.interior.surface_pressure_bar() * 1000.0, gravity_m_s2: b.gravity_m_s2() };
            let cells = sea_level_pressure(s.grid, &planet, &Air { wind: &wind.values, temperature: &season_t, elevation, land: &land });
            Map::Quantity(Quantity { cells, units: "hPa", contour: 4.0 })
        }
        View::SeaIce | View::Reefs | View::Productivity => {
            let sea = geoforge_surface::oceans::Sea::read(s.products, s.grid, s.body.exterior.radius_km, s.sea_world())?;
            if !s.sea_world() {
                return None;
            }
            match view {
                View::SeaIce => Map::Classes(tally(s.grid, &sea.ice().into_iter().map(Some).collect::<Vec<_>>(), SeaIce::info, &SeaIce::ALL)),
                View::Reefs => Map::Classes(tally(s.grid, &sea.reefs().into_iter().map(Some).collect::<Vec<_>>(), Reef::info, &Reef::ALL)),
                _ => chances(sea.productivity(), "how much life the sea holds, 0 (gyre desert) to 1 (upwelling coast)".into()),
            }
        }
    })
}

/// The magnetic pole and the auroral ovals' latitude, if the body has a field.
#[must_use]
pub fn aurora_of(body: &BodyContext) -> Option<(geoforge_core::geom::UnitVec3, f64)> {
    let f = MagneticField::of(body);
    f.aurora_lat_deg.map(|lat| (f.pole(), lat))
}

/// Index cells into a legend of the classes present (in `order`), with area shares.
fn tally<T: Copy + PartialEq>(grid: &AnyGrid, cells: &[Option<T>], info: impl Fn(T) -> (&'static str, &'static str, [u8; 3]), order: &[T]) -> Classes {
    let mut area = vec![0.0; order.len()];
    let index: Vec<Option<usize>> = cells.iter().map(|c| c.map(|v| order.iter().position(|&o| o == v).expect("every class is in the order"))).collect();
    for (c, k) in index.iter().enumerate() {
        if let Some(k) = *k {
            area[k] += grid.area_sr(c as Cell);
        }
    }
    let land: f64 = area.iter().sum();
    let present: Vec<usize> = (0..order.len()).filter(|&k| area[k] > 0.0).collect();
    let legend = present
        .iter()
        .map(|&k| {
            let (code, name, color) = info(order[k]);
            Class { code: code.into(), name: name.into(), color, area_share: area[k] / land }
        })
        .collect();
    let cells = index.into_iter().map(|c| c.map(|k| present.iter().position(|&p| p == k).expect("present"))).collect();
    Classes { cells, legend }
}

#[cfg(test)]
mod tests {
    use super::*;
    use geoforge_core::grid::GridSpec;

    #[test]
    fn names_round_trip() {
        let all = View::all();
        for v in &all {
            assert_eq!(View::from_name(v.name()), Some(*v));
        }
        let names: std::collections::BTreeSet<&str> = all.iter().map(|v| v.name()).collect();
        assert_eq!(names.len(), all.len(), "names are unique");
        let grouped: usize = ["climate_maps", "oceans", "resources", "hazards"].iter().map(|g| View::group(g).unwrap().len()).sum();
        assert_eq!(grouped, all.len(), "every view is in one group");
    }

    #[test]
    fn tally_shares_land_by_area() {
        let grid = GridSpec::Equirect { width: 4, height: 2 }.build().unwrap();
        // North row: two of one class, one of another, one sea; south row all the first class.
        let cells = [Some(1), Some(1), Some(2), None, Some(1), Some(1), Some(1), Some(1)];
        let c = tally(&grid, &cells, |v: i32| (["0", "1", "2"][v as usize], "", [0; 3]), &[0, 1, 2]);
        assert_eq!(c.legend.iter().map(|x| x.code.as_str()).collect::<Vec<_>>(), ["1", "2"]);
        assert!((c.legend[0].area_share - 6.0 / 7.0).abs() < 1e-9 && c.cells[3].is_none() && c.cells[2] == Some(1));
    }
}
