//! Sand seas and their dunes, from the winds, the rain and the land's liquids.
//!
//! Sand gathers into seas of dunes on dry, flat ground away from standing liquid: the Sahara's
//! and Arabia's ergs where the year's rain is under ~150 mm (Lancaster 1995), Titan's belt of
//! dunes across its dry tropics (Lorenz et al. 2006). A cell is a sand sea where the rain is
//! under 150 mm/yr, the ground slopes less than 1 in 50, and no lake or river lies within
//! 150 km. The sand itself is made by a liquid wearing rock down and carrying it (Titan's
//! organic sand too, Barnes et al. 2015), so a world without rivers has none here: Mars's
//! dunes, from older water and the wind's own abrasion, are not modelled, nor are a world's
//! sand supply, the wind's strength against the threshold that lifts sand, or vegetation
//! fixing old dunes.
//!
//! The dunes' form follows how the winds that move the sand agree (Fryberger 1979; Wasson &
//! Hyde 1983). Each season's wind moves sand as the cube of its speed; the share of that drift
//! left in the resultant (RDP/DP) separates one steady wind, which builds transverse ridges and
//! barchans (over 0.8), two winds at an angle, which build linear dunes along their mean
//! (0.3–0.8: the Namib's, the Simpson's, Titan's), and winds from all round, which pile star
//! dunes (under 0.3: the Grand Erg Oriental's). Two seasons (January and July) stand for the
//! year's winds.

use geoforge_core::grid::{Cell, Grid};
use geoforge_core::product::Products;
use std::collections::VecDeque;

/// Rain under which sand seas form, mm/yr.
pub const DRY_MM: f64 = 150.0;
/// Steepest ground sand seas cover, m/m.
pub const MAX_SLOPE: f64 = 0.02;
/// Nearest a lake or river may be, km.
pub const WET_KM: f64 = 150.0;
/// RDP/DP above which one wind dominates (Fryberger 1979).
const UNIMODAL: f64 = 0.8;
/// RDP/DP below which the winds come from all round.
const COMPLEX: f64 = 0.3;

/// What covers the land, as far as sand goes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Dunes {
    /// No sand sea.
    None,
    /// Transverse ridges and barchans: one steady wind.
    Transverse,
    /// Linear dunes along two winds' mean.
    Linear,
    /// Star dunes: winds from all round.
    Star,
}

impl Dunes {
    pub const ALL: [Self; 4] = [Self::None, Self::Transverse, Self::Linear, Self::Star];

    #[must_use]
    pub fn info(self) -> (&'static str, &'static str, [u8; 3]) {
        match self {
            Self::None => ("none", "No sand sea", [205, 200, 190]),
            Self::Transverse => ("transverse", "Transverse dunes and barchans", [240, 200, 110]),
            Self::Linear => ("linear", "Linear dunes", [215, 150, 70]),
            Self::Star => ("star", "Star dunes", [170, 90, 50]),
        }
    }

    /// The dunes that winds `a` and `b` (the year's two seasons, any units) build: by the share
    /// of their drift, ∝ speed³, left in the resultant.
    #[must_use]
    pub fn of_winds(a: [f64; 2], b: [f64; 2]) -> Self {
        let speed = |v: [f64; 2]| v[0].hypot(v[1]);
        let (sa, sb) = (speed(a), speed(b));
        let dp = sa.powi(3) + sb.powi(3);
        if dp <= 0.0 {
            return Self::None;
        }
        let rdp = (sa * sa * a[0] + sb * sb * b[0]).hypot(sa * sa * a[1] + sb * sb * b[1]);
        match rdp / dp {
            r if r > UNIMODAL => Self::Transverse,
            r if r >= COMPLEX => Self::Linear,
            _ => Self::Star,
        }
    }
}

/// Each cell's dunes, `None` at sea (on a world with seas); `None` overall while the winds,
/// rain, heights, rivers or lakes are missing.
#[must_use]
pub fn map(p: &Products, grid: &dyn Grid, radius_km: f64, sea: bool) -> Option<Vec<Option<Dunes>>> {
    let (jan, jul) = (&p.wind_jan()?.values, &p.wind_jul()?.values);
    let rain = &p.precipitation()?.values;
    let elev = &p.elevation()?.values;
    let (rivers, lakes) = (p.rivers()?, &p.lakes()?.owner.values);
    let n = grid.len();
    let mut liquid = vec![false; n];
    for (_, node) in rivers.nodes.iter() {
        liquid[node.cell as usize] = true;
    }
    for (c, &k) in lakes.iter().enumerate() {
        liquid[c] |= k != 0;
    }
    let wet = near(grid, &liquid, WET_KM / radius_km);
    let land = |c: usize| !sea || elev[c] >= 0.0;
    Some(
        (0..n)
            .map(|c| {
                if !land(c) {
                    return None;
                }
                let sandy = f64::from(rain[c]) < DRY_MM && !wet[c] && slope(grid, elev, radius_km, c as Cell) < MAX_SLOPE;
                let wind = |v: [f32; 2]| [f64::from(v[0]), f64::from(v[1])];
                Some(if sandy { Dunes::of_winds(wind(jan[c]), wind(jul[c])) } else { Dunes::None })
            })
            .collect(),
    )
}

/// Steepest drop from `c` to a neighbour, m/m.
fn slope(grid: &dyn Grid, elev: &[f32], radius_km: f64, c: Cell) -> f64 {
    let here = grid.center(c);
    grid.neighbors(c)
        .iter()
        .map(|nb| {
            let d_m = here.angle_to(grid.center(nb)) * radius_km * 1000.0;
            (f64::from(elev[nb as usize]) - f64::from(elev[c as usize])).abs() / d_m.max(1.0)
        })
        .fold(0.0, f64::max)
}

/// Cells within `reach` radians of a `source` cell: a search outward from every source that
/// keeps each cell's nearest one.
fn near(grid: &dyn Grid, source: &[bool], reach: f64) -> Vec<bool> {
    let mut from: Vec<Option<Cell>> = source.iter().enumerate().map(|(c, &s)| s.then_some(c as Cell)).collect();
    let mut queue: VecDeque<Cell> = (0..grid.len() as Cell).filter(|&c| source[c as usize]).collect();
    while let Some(c) = queue.pop_front() {
        let origin = from[c as usize].expect("queued cells have a source");
        for nb in grid.neighbors(c).iter() {
            if from[nb as usize].is_none() && grid.center(origin).angle_to(grid.center(nb)) <= reach {
                from[nb as usize] = Some(origin);
                queue.push_back(nb);
            }
        }
    }
    from.iter().map(Option::is_some).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_winds_shape_the_dunes() {
        // One wind all year: transverse; the same wind at half strength the other season too.
        assert_eq!(Dunes::of_winds([3.0, 0.0], [3.0, 0.0]), Dunes::Transverse);
        assert_eq!(Dunes::of_winds([3.0, 0.0], [1.5, 0.3]), Dunes::Transverse);
        // Two equal winds at right angles: linear, along their mean (RDP/DP = 0.71).
        assert_eq!(Dunes::of_winds([3.0, 0.0], [0.0, 3.0]), Dunes::Linear);
        // Opposite seasons: star dunes.
        assert_eq!(Dunes::of_winds([3.0, 0.0], [-2.8, 0.2]), Dunes::Star);
        assert_eq!(Dunes::of_winds([0.0, 0.0], [0.0, 0.0]), Dunes::None);
    }

    #[test]
    fn nearness_spreads_from_sources() {
        let grid = geoforge_core::grid::GridSpec::Equirect { width: 72, height: 36 }.build().unwrap();
        let mut source = vec![false; grid.len()];
        source[18 * 72] = true;
        let reach = 12f64.to_radians();
        let close = near(&grid, &source, reach);
        for c in 0..grid.len() as Cell {
            let d = grid.center(18 * 72).angle_to(grid.center(c));
            if d < 0.8 * reach {
                assert!(close[c as usize], "{c} at {}°", d.to_degrees());
            }
            if d > reach {
                assert!(!close[c as usize]);
            }
        }
    }
}
