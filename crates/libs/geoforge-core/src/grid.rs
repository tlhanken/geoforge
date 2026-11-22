//! Where surface data lives.
//!
//! Surface algorithms never touch pixel x/y. They ask a [`Grid`] for a cell's position,
//! area and neighbours on the unit sphere. [`EquirectGrid`] is the first implementation;
//! a Goldberg hex grid can be added as a second one without touching the algorithms.

use crate::geom::UnitVec3;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::collections::BinaryHeap;

/// Index of a cell in a grid.
pub type Cell = u32;

/// Serializable description of a grid; enough to rebuild it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum GridSpec {
    Equirect { width: u32, height: u32 },
    Goldberg { n: u32, m: u32 },
}

impl GridSpec {
    /// Build the grid this spec describes.
    ///
    /// # Errors
    /// Returns an error for grid kinds not implemented yet.
    pub fn build(&self) -> Result<AnyGrid, crate::Error> {
        match *self {
            GridSpec::Equirect { width, height } => Ok(AnyGrid::Equirect(EquirectGrid::new(width, height)?)),
            GridSpec::Goldberg { .. } => Err(crate::Error::Unsupported("Goldberg grids are not implemented yet".into())),
        }
    }

    #[must_use]
    pub fn cell_count(&self) -> usize {
        match *self {
            GridSpec::Equirect { width, height } => width as usize * height as usize,
            GridSpec::Goldberg { n, m } => 10 * (n * n + n * m + m * m) as usize + 2,
        }
    }
}

/// Up to eight neighbours of a cell, without allocation.
#[derive(Debug, Clone, Copy)]
pub struct Neighbors {
    buf: [Cell; 8],
    len: u8,
}

impl Neighbors {
    fn new() -> Self {
        Self { buf: [0; 8], len: 0 }
    }
    fn push(&mut self, c: Cell) {
        if self.as_slice().contains(&c) {
            return;
        }
        self.buf[self.len as usize] = c;
        self.len += 1;
    }
    #[must_use]
    pub fn as_slice(&self) -> &[Cell] {
        &self.buf[..self.len as usize]
    }
    pub fn iter(&self) -> std::iter::Copied<std::slice::Iter<'_, Cell>> {
        self.as_slice().iter().copied()
    }
    #[must_use]
    pub fn len(&self) -> usize {
        self.len as usize
    }
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }
}

/// A discretization of the unit sphere.
pub trait Grid: Send + Sync {
    fn len(&self) -> usize;
    fn center(&self, c: Cell) -> UnitVec3;
    /// Cell area in steradians. Multiply by R² for km².
    fn area_sr(&self, c: Cell) -> f64;
    fn neighbors(&self, c: Cell) -> Neighbors;
    fn locate(&self, p: UnitVec3) -> Cell;
    fn spec(&self) -> GridSpec;

    fn is_empty(&self) -> bool {
        self.len() == 0
    }
    /// Typical cell spacing in radians (used to size noise and filters).
    fn spacing_rad(&self) -> f64 {
        (4.0 * std::f64::consts::PI / self.len() as f64).sqrt()
    }
    fn cells(&self) -> std::ops::Range<Cell> {
        0..self.len() as Cell
    }
    /// East–west and north–south extent of a cell, radians. Defaults to a square cell.
    fn cell_extent_rad(&self, c: Cell) -> (f64, f64) {
        let s = self.area_sr(c).sqrt();
        (s, s)
    }
}

/// Equirectangular (plate carrée) grid: `width × height` cells, row-major from the north-west.
///
/// Neighbours are the 8 surrounding cells, wrapping east–west at the seam; cells in the
/// top and bottom rows also neighbour the three cells directly across the pole.
#[derive(Debug, Clone)]
pub struct EquirectGrid {
    width: u32,
    height: u32,
    lon_res: f64,
    lat_res: f64,
    row_area: Vec<f64>,
}

impl EquirectGrid {
    /// # Errors
    /// Width must be even and at least 4; height at least 2.
    pub fn new(width: u32, height: u32) -> Result<Self, crate::Error> {
        if width < 4 || height < 2 || width % 2 != 0 {
            return Err(crate::Error::Invalid(format!("equirect grid {width}x{height}: width must be even and >= 4, height >= 2")));
        }
        let lon_res = 360.0 / f64::from(width);
        let lat_res = 180.0 / f64::from(height);
        let dlon = lon_res.to_radians();
        let row_area = (0..height)
            .map(|y| {
                let top = (90.0 - f64::from(y) * lat_res).to_radians();
                let bot = (90.0 - f64::from(y + 1) * lat_res).to_radians();
                dlon * (top.sin() - bot.sin())
            })
            .collect();
        Ok(Self { width, height, lon_res, lat_res, row_area })
    }

    #[must_use]
    pub fn width(&self) -> u32 {
        self.width
    }
    #[must_use]
    pub fn height(&self) -> u32 {
        self.height
    }
    #[must_use]
    pub fn xy(&self, c: Cell) -> (u32, u32) {
        (c % self.width, c / self.width)
    }
    #[must_use]
    pub fn cell(&self, x: u32, y: u32) -> Cell {
        y * self.width + x
    }
    /// Pixel-centre latitude/longitude, exactly as geoforge v1 computed it.
    #[must_use]
    pub fn lat_lon_xy(&self, x: u32, y: u32) -> (f64, f64) {
        let lon = -180.0 + (f64::from(x) + 0.5) * self.lon_res;
        let lat = 90.0 - (f64::from(y) + 0.5) * self.lat_res;
        (lat, lon)
    }
    #[must_use]
    pub fn locate_lat_lon(&self, lat: f64, lon: f64) -> Cell {
        let x = ((lon + 180.0) / self.lon_res).floor().max(0.0) as u32;
        let y = ((90.0 - lat) / self.lat_res).floor().max(0.0) as u32;
        self.cell(x.min(self.width - 1), y.min(self.height - 1))
    }
    #[must_use]
    pub fn lat_res_deg(&self) -> f64 {
        self.lat_res
    }
    #[must_use]
    pub fn lon_res_deg(&self) -> f64 {
        self.lon_res
    }
}

impl Grid for EquirectGrid {
    fn len(&self) -> usize {
        self.width as usize * self.height as usize
    }

    fn center(&self, c: Cell) -> UnitVec3 {
        let (x, y) = self.xy(c);
        let (lat, lon) = self.lat_lon_xy(x, y);
        UnitVec3::from_lat_lon(lat, lon)
    }

    fn area_sr(&self, c: Cell) -> f64 {
        self.row_area[(c / self.width) as usize]
    }

    fn neighbors(&self, c: Cell) -> Neighbors {
        let (x, y) = self.xy(c);
        let w = self.width as i64;
        let mut n = Neighbors::new();
        for dy in -1i64..=1 {
            for dx in -1i64..=1 {
                if dx == 0 && dy == 0 {
                    continue;
                }
                let ny = i64::from(y) + dy;
                if ny < 0 || ny >= i64::from(self.height) {
                    continue;
                }
                let nx = (i64::from(x) + dx).rem_euclid(w);
                n.push(self.cell(nx as u32, ny as u32));
            }
        }
        if y == 0 || y == self.height - 1 {
            for d in -1i64..=1 {
                let nx = (i64::from(x) + w / 2 + d).rem_euclid(w);
                n.push(self.cell(nx as u32, y));
            }
        }
        n
    }

    fn locate(&self, p: UnitVec3) -> Cell {
        let ll = p.to_lat_lon();
        self.locate_lat_lon(ll.lat_deg, ll.lon_deg)
    }

    fn spec(&self) -> GridSpec {
        GridSpec::Equirect { width: self.width, height: self.height }
    }

    fn cell_extent_rad(&self, c: Cell) -> (f64, f64) {
        let dlat = self.lat_res.to_radians();
        // Width chosen so that width × height equals the exact cell area.
        (self.area_sr(c) / dlat, dlat)
    }
}

/// Any implemented grid (static dispatch).
#[derive(Debug, Clone)]
pub enum AnyGrid {
    Equirect(EquirectGrid),
}

impl AnyGrid {
    #[must_use]
    pub fn as_equirect(&self) -> Option<&EquirectGrid> {
        match self {
            AnyGrid::Equirect(g) => Some(g),
        }
    }
}

impl Grid for AnyGrid {
    fn len(&self) -> usize {
        match self {
            AnyGrid::Equirect(g) => g.len(),
        }
    }
    fn center(&self, c: Cell) -> UnitVec3 {
        match self {
            AnyGrid::Equirect(g) => g.center(c),
        }
    }
    fn area_sr(&self, c: Cell) -> f64 {
        match self {
            AnyGrid::Equirect(g) => g.area_sr(c),
        }
    }
    fn neighbors(&self, c: Cell) -> Neighbors {
        match self {
            AnyGrid::Equirect(g) => g.neighbors(c),
        }
    }
    fn locate(&self, p: UnitVec3) -> Cell {
        match self {
            AnyGrid::Equirect(g) => g.locate(p),
        }
    }
    fn spec(&self) -> GridSpec {
        match self {
            AnyGrid::Equirect(g) => g.spec(),
        }
    }
    fn cell_extent_rad(&self, c: Cell) -> (f64, f64) {
        match self {
            AnyGrid::Equirect(g) => g.cell_extent_rad(c),
        }
    }
}

// ---------------------------------------------------------------------------------------
// Grid algorithms shared by every layer.
// ---------------------------------------------------------------------------------------

#[derive(PartialEq)]
struct Item(f64, Cell);
impl Eq for Item {}
impl Ord for Item {
    fn cmp(&self, o: &Self) -> Ordering {
        o.0.partial_cmp(&self.0).unwrap_or(Ordering::Equal).then_with(|| o.1.cmp(&self.1))
    }
}
impl PartialOrd for Item {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}

/// Multi-source geodesic distance (radians, along neighbour steps) from `sources`, limited to
/// `max_rad` and to cells where `allow(cell)` is true. Unreached cells are `f64::INFINITY`.
/// Sources are always included (distance 0) even if `allow` rejects them.
pub fn distance_field(grid: &dyn Grid, sources: impl IntoIterator<Item = Cell>, max_rad: f64, mut allow: impl FnMut(Cell) -> bool) -> Vec<f64> {
    let mut dist = vec![f64::INFINITY; grid.len()];
    let mut heap = BinaryHeap::new();
    for s in sources {
        if dist[s as usize] > 0.0 {
            dist[s as usize] = 0.0;
            heap.push(Item(0.0, s));
        }
    }
    while let Some(Item(d, c)) = heap.pop() {
        if d > dist[c as usize] {
            continue;
        }
        let pc = grid.center(c);
        for n in grid.neighbors(c).iter() {
            if !allow(n) {
                continue;
            }
            let nd = d + pc.angle_to(grid.center(n));
            if nd <= max_rad && nd < dist[n as usize] {
                dist[n as usize] = nd;
                heap.push(Item(nd, n));
            }
        }
    }
    dist
}

/// Like [`distance_field`] but returns only the reached cells, in ascending cell order.
pub fn cells_within(grid: &dyn Grid, sources: &[Cell], max_rad: f64, allow: impl FnMut(Cell) -> bool) -> Vec<Cell> {
    // Sparse Dijkstra: avoids allocating a full-grid vector for small regions.
    let mut allow = allow;
    let mut best: std::collections::HashMap<Cell, f64> = std::collections::HashMap::new();
    let mut heap = BinaryHeap::new();
    for &s in sources {
        if best.insert(s, 0.0).is_none() {
            heap.push(Item(0.0, s));
        }
    }
    while let Some(Item(d, c)) = heap.pop() {
        if best.get(&c).is_some_and(|&b| d > b) {
            continue;
        }
        let pc = grid.center(c);
        for n in grid.neighbors(c).iter() {
            if !allow(n) {
                continue;
            }
            let nd = d + pc.angle_to(grid.center(n));
            if nd <= max_rad && best.get(&n).is_none_or(|&b| nd < b) {
                best.insert(n, nd);
                heap.push(Item(nd, n));
            }
        }
    }
    let mut out: Vec<Cell> = best.into_keys().collect();
    out.sort_unstable();
    out
}

/// Connected components of the cells where `member(cell)` is true, each sorted, the list
/// ordered by smallest cell.
pub fn components(grid: &dyn Grid, cells: &[Cell], member: impl Fn(Cell) -> bool) -> Vec<Vec<Cell>> {
    let mut seen = std::collections::HashSet::new();
    let mut out = Vec::new();
    let mut sorted = cells.to_vec();
    sorted.sort_unstable();
    for &start in &sorted {
        if !member(start) || !seen.insert(start) {
            continue;
        }
        let mut comp = vec![start];
        let mut i = 0;
        while i < comp.len() {
            let c = comp[i];
            i += 1;
            for n in grid.neighbors(c).iter() {
                if member(n) && seen.insert(n) {
                    comp.push(n);
                }
            }
        }
        comp.sort_unstable();
        out.push(comp);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn areas_sum_to_sphere() {
        let g = EquirectGrid::new(360, 180).unwrap();
        let total: f64 = g.cells().map(|c| g.area_sr(c)).sum();
        assert!((total - 4.0 * PI).abs() < 1e-9, "total {total}");
        // Equator cells are the largest, pole cells the smallest.
        assert!(g.area_sr(g.cell(0, 90)) > g.area_sr(g.cell(0, 0)) * 50.0);
    }

    #[test]
    fn seam_and_pole_neighbours() {
        let g = EquirectGrid::new(36, 18).unwrap();
        let west_edge = g.cell(0, 9);
        let n = g.neighbors(west_edge);
        assert_eq!(n.len(), 8);
        assert!(n.as_slice().contains(&g.cell(35, 9)), "wraps across the seam");
        let pole = g.cell(0, 0);
        let n = g.neighbors(pole);
        assert_eq!(n.len(), 8);
        assert!(n.as_slice().contains(&g.cell(18, 0)), "neighbours across the pole");
        // Symmetry: every neighbour relation goes both ways.
        for c in g.cells() {
            for m in g.neighbors(c).iter() {
                assert!(g.neighbors(m).as_slice().contains(&c), "{c} -> {m}");
            }
        }
    }

    #[test]
    fn locate_inverts_center() {
        let g = EquirectGrid::new(360, 180).unwrap();
        for c in (0..g.len() as Cell).step_by(97) {
            assert_eq!(g.locate(g.center(c)), c);
        }
    }

    #[test]
    fn distance_field_is_geodesic() {
        let g = EquirectGrid::new(360, 180).unwrap();
        let src = g.locate(UnitVec3::from_lat_lon(0.0, 0.0));
        let d = distance_field(&g, [src], PI, |_| true);
        let far = g.locate(UnitVec3::from_lat_lon(0.0, 90.0));
        // Grid paths are a little longer than the great circle, never shorter.
        let exact = g.center(src).angle_to(g.center(far));
        assert!(d[far as usize] >= exact - 1e-12);
        assert!(d[far as usize] < exact * 1.09, "{} vs {}", d[far as usize], exact);
    }
}
