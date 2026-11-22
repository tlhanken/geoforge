//! Derived views: maps that are pure functions of stored products. They are worked out when
//! rendered or exported and never stored, so they cannot go stale and cost nothing to keep.

use geoforge_core::grid::{AnyGrid, Cell, Grid};
use geoforge_core::product::{ProductId, Products};
use geoforge_surface::koppen::{self, Koppen, KoppenRules};
use serde::Serialize;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum View {
    /// Köppen–Geiger climate classes, from the climate layer.
    Koppen,
}

impl View {
    pub const ALL: [Self; 1] = [Self::Koppen];

    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Self::Koppen => "koppen",
        }
    }

    #[must_use]
    pub fn from_name(name: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|v| v.name() == name)
    }

    /// The products the view is worked out from.
    #[must_use]
    pub fn reads(self) -> &'static [ProductId] {
        match self {
            Self::Koppen => &[ProductId::Temperature, ProductId::TemperatureRange, ProductId::Precipitation, ProductId::PrecipitationWarmShare],
        }
    }
}

/// One class of a categorical view, with its share of the land.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Class {
    pub code: String,
    pub name: String,
    pub color: [u8; 3],
    pub land_share: f64,
}

/// A categorical view: each cell's class (an index into `legend`), `None` at sea.
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

/// Work out a view. `sea`: the world has surface water, so cells below 0 m are sea. `None`
/// while the products it reads are missing.
#[must_use]
pub fn classes(view: View, products: &Products, grid: &AnyGrid, sea: bool) -> Option<Classes> {
    match view {
        View::Koppen => {
            let map = koppen::map(products, sea, &KoppenRules::default())?;
            Some(tally(grid, &map, |k| Class { code: k.code().into(), name: k.name().into(), color: k.color(), land_share: 0.0 }, &Koppen::ALL))
        }
    }
}

/// Index cells into a legend of the classes present (in `order`), with area shares.
fn tally<T: Copy + PartialEq>(grid: &AnyGrid, cells: &[Option<T>], class: impl Fn(T) -> Class, order: &[T]) -> Classes {
    let mut area = vec![0.0; order.len()];
    let index: Vec<Option<usize>> = cells.iter().map(|c| c.map(|v| order.iter().position(|&o| o == v).expect("every class is in the order"))).collect();
    for (c, k) in index.iter().enumerate() {
        if let Some(k) = *k {
            area[k] += grid.area_sr(c as Cell);
        }
    }
    let land: f64 = area.iter().sum();
    let present: Vec<usize> = (0..order.len()).filter(|&k| area[k] > 0.0).collect();
    let legend = present.iter().map(|&k| Class { land_share: area[k] / land, ..class(order[k]) }).collect();
    let cells = index.into_iter().map(|c| c.map(|k| present.iter().position(|&p| p == k).expect("present"))).collect();
    Classes { cells, legend }
}

#[cfg(test)]
mod tests {
    use super::*;
    use geoforge_core::grid::GridSpec;

    #[test]
    fn names_round_trip() {
        for v in View::ALL {
            assert_eq!(View::from_name(v.name()), Some(v));
        }
    }

    #[test]
    fn tally_shares_land_by_area() {
        let grid = GridSpec::Equirect { width: 4, height: 2 }.build().unwrap();
        // North row: two of one class, one of another, one sea; south row all the first class.
        let cells = [Some(1), Some(1), Some(2), None, Some(1), Some(1), Some(1), Some(1)];
        let c = tally(&grid, &cells, |v: i32| Class { code: v.to_string(), name: String::new(), color: [0; 3], land_share: 0.0 }, &[0, 1, 2]);
        assert_eq!(c.legend.iter().map(|x| x.code.as_str()).collect::<Vec<_>>(), ["1", "2"]);
        assert!((c.legend[0].land_share - 6.0 / 7.0).abs() < 1e-9 && c.cells[3].is_none() && c.cells[2] == Some(1));
    }
}
