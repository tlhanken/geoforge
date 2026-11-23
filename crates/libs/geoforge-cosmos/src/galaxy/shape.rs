//! Where a component's stars are: unit-mass density shapes, how much of one lies in a box,
//! and draws from it.
//!
//! Positions are parsecs from the galaxy's centre, z toward its north pole.

use crate::draw::Draw;
use serde::{Deserialize, Serialize};
use std::f64::consts::{PI, TAU};

/// The shape of a component, of unit mass.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Shape {
    /// Exponential in radius and in height above the plane; the height may grow with age.
    Disk { scale_length_kpc: f64, height: Height },
    /// Exponential in an ellipsoidal radius √((x/a)² + (y/b)² + (z/c)²), its long axis
    /// `angle_deg` from the x axis (a bar or a triaxial bulge).
    Ellipsoid { scales_kpc: [f64; 3], angle_deg: f64 },
    /// Hernquist's (1990) sphere, flattened to axis ratio `q`: its projection follows de
    /// Vaucouleurs' law with an effective radius of 1.815 `a`.
    Hernquist { a_kpc: f64, q: f64 },
    /// Plummer's sphere (a cluster, a nucleus, a dwarf spheroidal).
    Plummer { a_pc: f64 },
    /// A flattened broken power law in m = √(R² + (z/q)²): (m² + core²)^(inner/2) to the
    /// break, continuing as m^outer, cut at `edge_kpc`.
    PowerLaw { inner: f64, outer: f64, break_kpc: f64, core_kpc: f64, q: f64, edge_kpc: f64 },
}

/// A disk's scale height, pc.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Height {
    Fixed {
        pc: f64,
    },
    /// Stars are born in a thin layer that thickens as they age: young + (old − young)(1 −
    /// e^(−age/τ)).
    Aging {
        young_pc: f64,
        old_pc: f64,
        tau_gyr: f64,
    },
}

impl Height {
    #[must_use]
    pub fn at(&self, age_gyr: f64) -> f64 {
        match *self {
            Self::Fixed { pc } => pc,
            Self::Aging { young_pc, old_pc, tau_gyr } => young_pc + (old_pc - young_pc) * (1.0 - (-age_gyr / tau_gyr).exp()),
        }
    }
}

/// A box of space: its corners and the mass inside, and for a disk the height its stars are
/// spread over vertically (drawn from the exponential, not evenly).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cell {
    pub lo: [f64; 3],
    pub hi: [f64; 3],
    pub mass: f64,
    pub height_pc: Option<f64>,
}

/// A square of a disk's plane: its corners, its surface mass, and that mass weighted by the
/// arms' excess (so the mass under arms of any contrast c is `mass + c · excess`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlaneLeaf {
    pub lo: [f64; 2],
    pub hi: [f64; 2],
    pub mass: f64,
    pub excess: f64,
}

/// A shape made ready to use: normalised, at one stellar age.
#[derive(Debug, Clone)]
pub struct Prepared {
    shape: Shape,
    /// Scale height for a disk at this age, pc.
    height: f64,
    /// Normalisation of the density.
    norm: f64,
    /// For a power law: its value at the break from the inner law, and the cumulative mass
    /// by radius for draws.
    break_value: f64,
    cumulative: Vec<(f64, f64)>,
    /// Disks only: the radius beyond which the galaxy has no stars, pc.
    edge_pc: f64,
}

/// Gauss–Legendre nodes for two points on [−1, 1].
const GAUSS: [f64; 2] = [-0.577_350_269_189_625_8, 0.577_350_269_189_625_8];
/// Deepest subdivision of a box around a concentrated shape's centre.
const MAX_DEPTH: u32 = 20;

impl Prepared {
    /// `shape` for stars of age `age_gyr`, disks cut at `edge_kpc`.
    #[must_use]
    pub fn new(shape: &Shape, age_gyr: f64, edge_kpc: f64) -> Self {
        let mut p = Self { shape: shape.clone(), height: 0.0, norm: 1.0, break_value: 0.0, cumulative: Vec::new(), edge_pc: edge_kpc * 1000.0 };
        match *shape {
            Shape::Disk { scale_length_kpc, height } => {
                p.height = height.at(age_gyr);
                let rd = scale_length_kpc * 1000.0;
                p.norm = 1.0 / (TAU * rd * rd * 2.0 * p.height);
            }
            Shape::Ellipsoid { scales_kpc, .. } => p.norm = 1.0 / (8.0 * PI * scales_kpc.iter().map(|s| s * 1000.0).product::<f64>()),
            Shape::Hernquist { a_kpc, q } => p.norm = a_kpc * 1000.0 / (TAU * q),
            Shape::Plummer { a_pc } => p.norm = 3.0 / (4.0 * PI * a_pc.powi(3)),
            Shape::PowerLaw { inner, break_kpc, core_kpc, q, edge_kpc, .. } => {
                let (rb, c) = (break_kpc * 1000.0, core_kpc * 1000.0);
                p.break_value = (rb * rb + c * c).powf(0.5 * inner);
                // The mass within m, on a log grid: 4π q ∫ m² f(m) dm.
                let (lo, hi) = (c * 1e-2, edge_kpc * 1000.0);
                let n = 400;
                let mut total = 0.0;
                let mut prev = (lo, lo.powi(2) * p.radial(lo));
                p.cumulative.push((0.0, 0.0));
                for i in 1..=n {
                    let m = lo * (hi / lo).powf(f64::from(i) / f64::from(n));
                    let f = m * m * p.radial(m);
                    total += 2.0 * PI * q * (f + prev.1) * (m - prev.0);
                    p.cumulative.push((m, total));
                    prev = (m, f);
                }
                p.norm = 1.0 / total;
                p.cumulative.iter_mut().for_each(|c| c.1 /= total);
            }
        }
        p
    }

    /// The power law's radial profile, unnormalised.
    fn radial(&self, m: f64) -> f64 {
        let Shape::PowerLaw { inner, outer, break_kpc, core_kpc, edge_kpc, .. } = self.shape else { return 0.0 };
        let (rb, c) = (break_kpc * 1000.0, core_kpc * 1000.0);
        if m > edge_kpc * 1000.0 {
            0.0
        } else if m <= rb {
            (m * m + c * c).powf(0.5 * inner)
        } else {
            self.break_value * (m / rb).powf(outer)
        }
    }

    /// The disk's scale height at this age, pc (0 for other shapes).
    #[must_use]
    pub fn height_pc(&self) -> f64 {
        self.height
    }

    /// Unit-mass density at `p`, per pc³.
    #[must_use]
    pub fn density(&self, p: [f64; 3]) -> f64 {
        match self.shape {
            Shape::Disk { .. } => self.in_plane(p[0], p[1]) * self.vertical(p[2]),
            _ => self.spheroid(p),
        }
    }

    /// A disk's surface density (its density is this times its vertical profile).
    fn in_plane(&self, x: f64, y: f64) -> f64 {
        match self.shape {
            Shape::Disk { scale_length_kpc, .. } => {
                let r = x.hypot(y);
                if r > self.edge_pc {
                    0.0
                } else {
                    self.norm * 2.0 * self.height * (-r / (scale_length_kpc * 1000.0)).exp()
                }
            }
            _ => 1.0,
        }
    }

    fn vertical(&self, z: f64) -> f64 {
        match self.shape {
            Shape::Disk { .. } => (-z.abs() / self.height).exp() / (2.0 * self.height),
            _ => 0.0,
        }
    }

    /// Unit-mass density of a shape that is not a disk.
    fn spheroid(&self, p: [f64; 3]) -> f64 {
        match self.shape {
            Shape::Disk { .. } => self.in_plane(p[0], p[1]) * self.vertical(p[2]),
            Shape::Ellipsoid { scales_kpc, angle_deg } => {
                let (s, c) = angle_deg.to_radians().sin_cos();
                let (u, v) = (c * p[0] + s * p[1], -s * p[0] + c * p[1]);
                let m = ((u / (scales_kpc[0] * 1000.0)).powi(2) + (v / (scales_kpc[1] * 1000.0)).powi(2) + (p[2] / (scales_kpc[2] * 1000.0)).powi(2)).sqrt();
                self.norm * (-m).exp()
            }
            Shape::Hernquist { a_kpc, q } => {
                let a = a_kpc * 1000.0;
                let m = (p[0] * p[0] + p[1] * p[1] + (p[2] / q).powi(2)).sqrt().max(1e-3 * a);
                self.norm / (m * (m + a).powi(3))
            }
            Shape::Plummer { a_pc } => self.norm * (1.0 + (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]) / (a_pc * a_pc)).powf(-2.5),
            Shape::PowerLaw { q, .. } => self.norm * self.radial((p[0] * p[0] + p[1] * p[1] + (p[2] / q).powi(2)).sqrt()),
        }
    }

    /// The scale under which a box near the centre is split further, pc.
    fn fine_scale(&self) -> f64 {
        match self.shape {
            Shape::Disk { scale_length_kpc, .. } => 250.0f64.min(500.0 * scale_length_kpc),
            Shape::Ellipsoid { scales_kpc, .. } => 500.0 * scales_kpc.iter().copied().fold(f64::MAX, f64::min),
            Shape::Hernquist { a_kpc, .. } => 50.0 * a_kpc,
            Shape::Plummer { a_pc } => 0.5 * a_pc,
            Shape::PowerLaw { core_kpc, .. } => 500.0 * core_kpc,
        }
    }

    /// Is it a disk (whose cells are made from its plane)?
    #[must_use]
    pub fn is_disk(&self) -> bool {
        matches!(self.shape, Shape::Disk { .. })
    }

    /// The mass in the box from `lo` to `hi` as cells, `factor(x, y)` scaling a disk's density
    /// in the plane (spiral arms; resolved to `detail_pc`).
    pub fn cells(&self, lo: [f64; 3], hi: [f64; 3], factor: &dyn Fn(f64, f64) -> f64, detail_pc: f64) -> Vec<Cell> {
        if self.is_disk() {
            let plane = self.plane([lo[0], lo[1]], [hi[0], hi[1]], &|x, y| factor(x, y) - 1.0, detail_pc);
            return self.disk_cells(&plane, 1.0, lo[2], hi[2]);
        }
        let mut out = Vec::new();
        self.split(lo, hi, 0, &mut out);
        out
    }

    /// A disk's square of plane as leaves: split finely near the centre, and to `detail_pc`
    /// within the disk where the arms shape it; `excess(x, y)` is the arms' excess there. The
    /// plane is the same for every age of a disk, so one serves all its slices.
    pub fn plane(&self, lo: [f64; 2], hi: [f64; 2], excess: &dyn Fn(f64, f64) -> f64, detail_pc: f64) -> Vec<PlaneLeaf> {
        let mut out = Vec::new();
        self.split_plane(lo, hi, excess, detail_pc, 0, &mut out);
        out
    }

    /// A disk's cells between heights `lo_z` and `hi_z` from its plane's leaves, under arms of
    /// `contrast`.
    #[must_use]
    pub fn disk_cells(&self, plane: &[PlaneLeaf], contrast: f64, lo_z: f64, hi_z: f64) -> Vec<Cell> {
        let z = vertical_share(lo_z, hi_z, self.height);
        if z <= 0.0 {
            return Vec::new();
        }
        plane
            .iter()
            .filter_map(|l| {
                let mass = (l.mass + contrast * l.excess).max(0.0) * z;
                (mass > 0.0).then_some(Cell { lo: [l.lo[0], l.lo[1], lo_z], hi: [l.hi[0], l.hi[1], hi_z], mass, height_pc: Some(self.height) })
            })
            .collect()
    }

    fn split_plane(&self, lo: [f64; 2], hi: [f64; 2], excess: &dyn Fn(f64, f64) -> f64, detail_pc: f64, depth: u32, out: &mut Vec<PlaneLeaf>) {
        let side = (hi[0] - lo[0]).max(hi[1] - lo[1]);
        let dist = (lo[0].max(0.0) + (-hi[0]).max(0.0)).hypot(lo[1].max(0.0) + (-hi[1]).max(0.0));
        if dist > self.edge_pc {
            return;
        }
        let centre = side > self.fine_scale() && dist < 2.0 * side;
        let arms = side > detail_pc;
        if depth < MAX_DEPTH && (centre || arms) {
            let mid = [0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1])];
            for (a, b) in [((lo[0], mid[0]), (lo[1], mid[1])), ((mid[0], hi[0]), (lo[1], mid[1])), ((lo[0], mid[0]), (mid[1], hi[1])), ((mid[0], hi[0]), (mid[1], hi[1]))] {
                self.split_plane([a.0, b.0], [a.1, b.1], excess, detail_pc, depth + 1, out);
            }
            return;
        }
        let (dx, dy) = (hi[0] - lo[0], hi[1] - lo[1]);
        let (mut mass, mut weighted) = (0.0, 0.0);
        for gx in GAUSS {
            for gy in GAUSS {
                let (x, y) = (lo[0] + 0.5 * dx * (1.0 + gx), lo[1] + 0.5 * dy * (1.0 + gy));
                let f = self.in_plane(x, y);
                mass += f;
                weighted += f * excess(x, y);
            }
        }
        if mass > 0.0 {
            out.push(PlaneLeaf { lo, hi, mass: mass / 4.0 * dx * dy, excess: weighted / 4.0 * dx * dy });
        }
    }

    fn split(&self, lo: [f64; 3], hi: [f64; 3], depth: u32, out: &mut Vec<Cell>) {
        let side = (0..3).map(|k| hi[k] - lo[k]).fold(0.0, f64::max);
        // Distance from the centre to the box.
        let dist = (0..3).map(|k| (lo[k].max(0.0) + (-hi[k]).max(0.0)).powi(2)).sum::<f64>().sqrt();
        if depth < MAX_DEPTH && side > self.fine_scale() && dist < 2.0 * side {
            let mid = [0, 1, 2].map(|k| 0.5 * (lo[k] + hi[k]));
            for c in 0..8 {
                let pick = |k: usize, bit: usize| if c >> bit & 1 == 0 { (lo[k], mid[k]) } else { (mid[k], hi[k]) };
                let (a, b, cc) = (pick(0, 0), pick(1, 1), pick(2, 2));
                self.split([a.0, b.0, cc.0], [a.1, b.1, cc.1], depth + 1, out);
            }
            return;
        }
        let d = [0, 1, 2].map(|k| hi[k] - lo[k]);
        let mut sum = 0.0;
        for gx in GAUSS {
            for gy in GAUSS {
                for gz in GAUSS {
                    let p = [0, 1, 2].map(|k| lo[k] + 0.5 * d[k] * (1.0 + [gx, gy, gz][k]));
                    sum += self.spheroid(p);
                }
            }
        }
        let mass = sum / 8.0 * d[0] * d[1] * d[2];
        if mass > 0.0 {
            out.push(Cell { lo, hi, mass, height_pc: None });
        }
    }

    /// A position drawn from the whole shape, with `accept(p)` (0 to 1) thinning it (spiral
    /// arms).
    pub(crate) fn draw(&self, d: &mut Draw, accept: &dyn Fn([f64; 3]) -> f64) -> [f64; 3] {
        loop {
            let p = self.draw_once(d);
            if p.iter().all(|v| v.is_finite()) && d.u() < accept(p) {
                return p;
            }
        }
    }

    fn draw_once(&self, d: &mut Draw) -> [f64; 3] {
        let dir = |d: &mut Draw| {
            let z = 2.0 * d.u() - 1.0;
            let phi = TAU * d.u();
            let s = (1.0 - z * z).sqrt();
            [s * phi.cos(), s * phi.sin(), z]
        };
        match self.shape {
            Shape::Disk { scale_length_kpc, .. } => {
                // Surface density ∝ e^(−R/h): radii follow R e^(−R/h), a gamma of shape 2.
                let r = -scale_length_kpc * 1000.0 * (d.u().max(1e-12).ln() + d.u().max(1e-12).ln());
                if r > self.edge_pc {
                    return [f64::NAN; 3];
                }
                let phi = TAU * d.u();
                let z = -self.height * d.u().max(1e-12).ln() * if d.u() < 0.5 { 1.0 } else { -1.0 };
                [r * phi.cos(), r * phi.sin(), z]
            }
            Shape::Ellipsoid { scales_kpc, angle_deg } => {
                // e^(−m) m² dm: a gamma of shape 3.
                let m = -(d.u().max(1e-12).ln() + d.u().max(1e-12).ln() + d.u().max(1e-12).ln());
                let u = dir(d);
                let (s, c) = angle_deg.to_radians().sin_cos();
                let (a, b) = (u[0] * m * scales_kpc[0] * 1000.0, u[1] * m * scales_kpc[1] * 1000.0);
                [c * a - s * b, s * a + c * b, u[2] * m * scales_kpc[2] * 1000.0]
            }
            Shape::Hernquist { a_kpc, q } => {
                let s = d.u().sqrt();
                let r = a_kpc * 1000.0 * s / (1.0 - s).max(1e-12);
                let u = dir(d);
                [u[0] * r, u[1] * r, u[2] * r * q]
            }
            Shape::Plummer { a_pc } => {
                let r = a_pc / (d.u().max(1e-12).powf(-2.0 / 3.0) - 1.0).max(1e-12).sqrt();
                let u = dir(d);
                [u[0] * r, u[1] * r, u[2] * r]
            }
            Shape::PowerLaw { q, .. } => {
                let v = d.u();
                let i = self.cumulative.partition_point(|c| c.1 < v).clamp(1, self.cumulative.len() - 1);
                let (a, b) = (self.cumulative[i - 1], self.cumulative[i]);
                let m = a.0 + (b.0 - a.0) * ((v - a.1) / (b.1 - a.1).max(1e-300));
                let u = dir(d);
                [u[0] * m, u[1] * m, u[2] * m * q]
            }
        }
    }
}

/// The share of an exponential layer of scale height `h` between heights `lo` and `hi`.
fn vertical_share(lo: f64, hi: f64, h: f64) -> f64 {
    // Cumulative of e^(−|z|/h)/(2h).
    let cdf = |z: f64| if z < 0.0 { 0.5 * (z / h).exp() } else { 1.0 - 0.5 * (-z / h).exp() };
    (cdf(hi) - cdf(lo)).max(0.0)
}

/// A height drawn from an exponential layer of scale height `h` between `lo` and `hi`.
pub(crate) fn draw_height(d: &mut Draw, lo: f64, hi: f64, h: f64) -> f64 {
    let cdf = |z: f64| if z < 0.0 { 0.5 * (z / h).exp() } else { 1.0 - 0.5 * (-z / h).exp() };
    let v = cdf(lo) + d.u() * (cdf(hi) - cdf(lo));
    let z = if v < 0.5 { h * (2.0 * v).max(1e-300).ln() } else { -h * (2.0 * (1.0 - v)).max(1e-300).ln() };
    z.clamp(lo, hi)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn shapes() -> Vec<Shape> {
        vec![
            Shape::Disk { scale_length_kpc: 2.6, height: Height::Aging { young_pc: 50.0, old_pc: 350.0, tau_gyr: 3.0 } },
            Shape::Ellipsoid { scales_kpc: [0.7, 0.44, 0.18], angle_deg: -27.0 },
            Shape::Hernquist { a_kpc: 2.75, q: 0.7 },
            Shape::Plummer { a_pc: 3.2 },
            Shape::PowerLaw { inner: -2.5, outer: -3.7, break_kpc: 25.0, core_kpc: 1.0, q: 0.65, edge_kpc: 100.0 },
        ]
    }

    /// The cells of a box holding the whole shape add up to its whole mass, and a box's
    /// mass splits between its halves.
    #[test]
    fn boxes_hold_the_mass() {
        for shape in shapes() {
            let p = Prepared::new(&shape, 5.0, 1000.0);
            let big = 400_000.0;
            let whole: f64 = p.cells([-big; 3], [big; 3], &|_, _| 1.0, 1000.0).iter().map(|c| c.mass).sum();
            assert!((whole - 1.0).abs() < 0.03, "{shape:?}: {whole}");
            let (lo, hi) = ([-3000.0, -1000.0, -500.0], [5000.0, 1500.0, 700.0]);
            let mass = |lo, hi| p.cells(lo, hi, &|_, _| 1.0, 100.0).iter().map(|c| c.mass).sum::<f64>();
            let halves = mass(lo, [1000.0, hi[1], hi[2]]) + mass([1000.0, lo[1], lo[2]], hi);
            assert!((mass(lo, hi) - halves).abs() < 0.02 * mass(lo, hi) + 1e-12, "{shape:?}");
        }
    }

    /// Draws follow the density: the share drawn into a box matches its mass.
    #[test]
    fn draws_follow_the_density() {
        for shape in shapes() {
            let p = Prepared::new(&shape, 5.0, 1000.0);
            let (lo, hi) = ([0.0, 0.0, -200.0], [3000.0, 3000.0, 200.0]);
            let mass: f64 = p.cells(lo, hi, &|_, _| 1.0, 100.0).iter().map(|c| c.mass).sum();
            let mut d = Draw::new(geoforge_core::seed::Seed(3));
            let n = 40_000;
            let inside = (0..n).filter(|_| {
                let q = p.draw(&mut d, &|_| 1.0);
                (0..3).all(|k| (lo[k]..hi[k]).contains(&q[k]))
            });
            let share = inside.count() as f64 / f64::from(n);
            assert!((share - mass).abs() < 4.0 * (mass / f64::from(n)).sqrt() + 0.01 * mass, "{shape:?}: drawn {share}, mass {mass}");
        }
    }

    #[test]
    fn heights_grow_with_age() {
        let h = Height::Aging { young_pc: 50.0, old_pc: 350.0, tau_gyr: 3.0 };
        assert!(h.at(0.0) == 50.0 && h.at(3.0) > 200.0 && h.at(11.0) < 350.0);
        let mut d = Draw::new(geoforge_core::seed::Seed(1));
        let zs: Vec<f64> = (0..5000).map(|_| draw_height(&mut d, 0.0, 100.0, 50.0)).collect();
        assert!(zs.iter().all(|z| (0.0..=100.0).contains(z)));
        let mean = zs.iter().sum::<f64>() / zs.len() as f64;
        // The mean of an exponential of 50 pc cut at 100 pc: 50(1 − 3e⁻²)/(1 − e⁻²) ≈ 34.3.
        assert!((mean - 34.3).abs() < 1.5, "{mean}");
    }
}
