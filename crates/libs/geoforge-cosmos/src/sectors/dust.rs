//! Dust: a smooth exponential layer, and clouds where stars form.
//!
//! - **The layer** dims by the galaxy's `mag_per_kpc` in the plane at its reference radius,
//!   falling with its scale length and height; along a line of sight it is integrated
//!   exactly in height, piece by piece in radius.
//! - **Clouds** form with the young stars: per 512 pc cube, the cube's young stellar mass
//!   times the galaxy's `cloud_gas` in clouds of dN/dM ∝ M⁻¹·⁶ from 10³ to 3×10⁶ M☉ (the
//!   molecular cloud mass function's slope of 1.5–1.8, Rosolowsky 2005), each a Gaussian of
//!   σ = 3 pc × (M/10⁴ M☉)^½, so about 90 M☉ pc⁻² within 2σ, the near-constant column of
//!   molecular clouds (Larson 1981). Gas dims by 0.047 mag per M☉ pc⁻²: A_V = N_H / 1.9×10²¹
//!   cm⁻² (Bohlin et al. 1978), with 1 M☉ pc⁻² ≈ 8.9×10¹⁹ hydrogen atoms cm⁻².

use super::bins;
use super::field::{cube_seed, place};
use super::groups::LEVEL;
use super::Sectors;
use crate::draw::Draw;
use std::f64::consts::{PI, SQRT_2};

const CLOUD_MASS: (f64, f64) = (1e3, 3e6);
const CLOUD_SLOPE: f64 = 1.6;
const MAG_PER_MSUN_PC2: f64 = 0.047;
/// Clouds whose cubes lie farther than this from the plane are not looked for, pc.
const CLOUD_LAYER_PC: f64 = 400.0;
/// Cells of the angular index of clouds, degrees.
const INDEX_DEG: f64 = 2.0;

/// A molecular cloud.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cloud {
    pub centre: [f64; 3],
    pub mass_msun: f64,
    pub sigma_pc: f64,
}

impl Cloud {
    /// Its dimming along a line from `from` in unit direction `u` out to `length`, mag.
    fn along(&self, from: [f64; 3], u: [f64; 3], length: f64) -> f64 {
        let c = [0, 1, 2].map(|k| self.centre[k] - from[k]);
        let t = c[0] * u[0] + c[1] * u[1] + c[2] * u[2];
        let b2 = (c[0] * c[0] + c[1] * c[1] + c[2] * c[2] - t * t).max(0.0);
        let s = self.sigma_pc;
        if b2 > 16.0 * s * s {
            return 0.0;
        }
        let k = SQRT_2 * s;
        MAG_PER_MSUN_PC2 * self.mass_msun * (-b2 / (2.0 * s * s)).exp() / (4.0 * PI * s * s) * (erf((length - t) / k) + erf(t / k))
    }
}

/// The mean of dN/dM ∝ M^(−slope) between the cloud masses.
fn mean_cloud_mass() -> f64 {
    let (lo, hi) = CLOUD_MASS;
    let a = CLOUD_SLOPE;
    ((hi.powf(2.0 - a) - lo.powf(2.0 - a)) / (2.0 - a)) / ((hi.powf(1.0 - a) - lo.powf(1.0 - a)) / (1.0 - a))
}

impl Sectors {
    /// The clouds of group cube `cube`.
    #[must_use]
    pub fn clouds(&self, cube: [i32; 3]) -> Vec<Cloud> {
        let gas = self.galaxy().dust.cloud_gas;
        if gas <= 0.0 {
            return Vec::new();
        }
        let (lo, hi) = bins::corners(cube, LEVEL);
        let cells = self.young_cells(lo, hi);
        let masses: Vec<f64> = self.young.iter().zip(&cells).map(|(yp, c)| yp.mass_msun * c.iter().map(|c| c.mass).sum::<f64>()).collect();
        let total: f64 = masses.iter().sum();
        let mut d = Draw::new(cube_seed(self.seed(), b"clouds", 0, cube));
        let (a, (lo_m, hi_m)) = (1.0 - CLOUD_SLOPE, CLOUD_MASS);
        (0..d.poisson(total * gas / mean_cloud_mass()))
            .map(|_| {
                let mut u = d.u() * total;
                let p = masses.iter().position(|&m| {
                    u -= m;
                    u <= 0.0
                });
                let p = p.unwrap_or(masses.len() - 1);
                let share = masses[p] / self.young[p].mass_msun;
                let centre = place(&mut d, &cells[p], share);
                let mass_msun = (lo_m.powf(a) + d.u() * (hi_m.powf(a) - lo_m.powf(a))).powf(1.0 / a);
                Cloud { centre, mass_msun, sigma_pc: 3.0 * (mass_msun / 1e4).sqrt() }
            })
            .collect()
    }
}

/// The dust as seen from one place: the smooth layer, and the clouds within a radius indexed
/// by direction.
pub struct DustView {
    from: [f64; 3],
    /// Smooth layer: mag per pc at the reference radius in the plane, scale length and height,
    /// reference radius (pc).
    per_pc: f64,
    length: f64,
    height: f64,
    reference: f64,
    clouds: Vec<Cloud>,
    index: Vec<Vec<u32>>,
}

const LON_CELLS: usize = (360.0 / INDEX_DEG) as usize;
const LAT_CELLS: usize = (180.0 / INDEX_DEG) as usize;

/// The index cell of a direction (in the galaxy's frame).
fn cell_of(u: [f64; 3]) -> (usize, usize) {
    let lat = u[2].clamp(-1.0, 1.0).asin().to_degrees();
    let lon = u[1].atan2(u[0]).to_degrees().rem_euclid(360.0);
    (((lon / INDEX_DEG) as usize).min(LON_CELLS - 1), (((lat + 90.0) / INDEX_DEG) as usize).min(LAT_CELLS - 1))
}

impl DustView {
    /// The dust seen from `from` (pc), clouds within `radius_pc`.
    #[must_use]
    pub fn new(sectors: &Sectors, from: [f64; 3], radius_pc: f64) -> Self {
        let g = sectors.galaxy();
        let mut clouds = Vec::new();
        for cube in bins::cubes_near(from, radius_pc, LEVEL) {
            let (lo, hi) = bins::corners(cube, LEVEL);
            if lo[2] > CLOUD_LAYER_PC || hi[2] < -CLOUD_LAYER_PC {
                continue;
            }
            clouds.extend(sectors.clouds(cube));
        }
        let mut index = vec![Vec::new(); LON_CELLS * LAT_CELLS];
        for (i, c) in clouds.iter().enumerate() {
            let v = [0, 1, 2].map(|k| c.centre[k] - from[k]);
            let d = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            let reach = 4.0 * c.sigma_pc;
            if d <= reach {
                index.iter_mut().for_each(|cell| cell.push(i as u32));
                continue;
            }
            let u = v.map(|x| x / d);
            let alpha = (reach / d).asin().to_degrees() + INDEX_DEG;
            let lat = u[2].asin().to_degrees();
            let lon = u[1].atan2(u[0]).to_degrees();
            let (lat_lo, lat_hi) = ((lat - alpha).max(-90.0), (lat + alpha).min(90.0));
            let mut row = lat_lo;
            while row <= lat_hi + 1e-9 {
                let y = (((row + 90.0) / INDEX_DEG) as usize).min(LAT_CELLS - 1);
                let half = if row.abs() + INDEX_DEG >= 90.0 { 180.0 } else { (alpha / row.to_radians().cos().max(1e-3)).min(180.0) };
                let mut col = lon - half;
                while col <= lon + half + 1e-9 {
                    let x = ((col.rem_euclid(360.0) / INDEX_DEG) as usize).min(LON_CELLS - 1);
                    let cell = &mut index[y * LON_CELLS + x];
                    if cell.last() != Some(&(i as u32)) {
                        cell.push(i as u32);
                    }
                    col += INDEX_DEG;
                }
                row += INDEX_DEG;
            }
        }
        Self { from, per_pc: g.dust.mag_per_kpc / 1000.0, length: g.dust.scale_length_kpc * 1000.0, height: g.dust.height_pc, reference: g.reference_kpc * 1000.0, clouds, index }
    }

    /// How many clouds it holds.
    #[must_use]
    pub fn cloud_count(&self) -> usize {
        self.clouds.len()
    }

    /// Dimming from here to `p`, mag.
    #[must_use]
    pub fn to(&self, p: [f64; 3]) -> f64 {
        let v = [0, 1, 2].map(|k| p[k] - self.from[k]);
        let d = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        if d <= 0.0 {
            return 0.0;
        }
        let u = v.map(|x| x / d);
        self.smooth(u, 0.0, d) + self.cloudy(u, d)
    }

    /// Dimming along unit direction `u` at each of `distances` (rising), mag.
    #[must_use]
    pub fn along(&self, u: [f64; 3], distances: &[f64]) -> Vec<f64> {
        let (x, y) = cell_of(u);
        let near: Vec<&Cloud> = self.index[y * LON_CELLS + x].iter().map(|&i| &self.clouds[i as usize]).collect();
        let mut smooth = 0.0;
        let mut prev = 0.0;
        distances
            .iter()
            .map(|&s| {
                smooth += self.smooth(u, prev, s);
                prev = s;
                smooth + near.iter().map(|c| c.along(self.from, u, s)).sum::<f64>()
            })
            .collect()
    }

    fn cloudy(&self, u: [f64; 3], length: f64) -> f64 {
        let (x, y) = cell_of(u);
        self.index[y * LON_CELLS + x].iter().map(|&i| self.clouds[i as usize].along(self.from, u, length)).sum()
    }

    /// The smooth layer between distances `s0` and `s1` along `u`: in pieces short against
    /// the scale length, each exact in height.
    fn smooth(&self, u: [f64; 3], s0: f64, s1: f64) -> f64 {
        if self.per_pc <= 0.0 || s1 <= s0 {
            return 0.0;
        }
        let pieces = ((s1 - s0) / (0.25 * self.length)).ceil().clamp(1.0, 32.0) as usize;
        let step = (s1 - s0) / pieces as f64;
        let h = self.height;
        let g = |z: f64| z.signum() * h * (1.0 - (-z.abs() / h).exp());
        (0..pieces)
            .map(|i| {
                let (a, b) = (s0 + step * i as f64, s0 + step * (i + 1) as f64);
                let mid = 0.5 * (a + b);
                let r = (self.from[0] + mid * u[0]).hypot(self.from[1] + mid * u[1]);
                let radial = (-(r - self.reference) / self.length).exp();
                let (za, zb) = (self.from[2] + a * u[2], self.from[2] + b * u[2]);
                let vertical = if u[2].abs() < 1e-9 { step * (-za.abs() / h).exp() } else { (g(zb) - g(za)).abs() / u[2].abs() };
                self.per_pc * radial * vertical
            })
            .sum()
    }
}

/// The error function (Abramowitz & Stegun 7.1.26, to 1.5×10⁻⁷).
fn erf(x: f64) -> f64 {
    let t = 1.0 / (1.0 + 0.327_591_1 * x.abs());
    let y = 1.0 - (((((1.061_405_429 * t - 1.453_152_027) * t) + 1.421_413_741) * t - 0.284_496_736) * t + 0.254_829_592) * t * (-x * x).exp();
    y.copysign(x)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::galaxy::Galaxy;
    use geoforge_core::seed::Seed;

    #[test]
    fn a_cloud_dims_what_is_behind_it() {
        let c = Cloud { centre: [100.0, 0.0, 0.0], mass_msun: 1e5, sigma_pc: 3.0 * 10f64.sqrt() };
        let through = c.along([0.0; 3], [1.0, 0.0, 0.0], 1000.0);
        let central = MAG_PER_MSUN_PC2 * 1e5 / (2.0 * PI * c.sigma_pc.powi(2));
        assert!((through - central).abs() < 1e-6 * central && through > 3.0, "{through}");
        assert!(c.along([0.0; 3], [1.0, 0.0, 0.0], 50.0) < 1e-3, "in front");
        assert!((c.along([0.0; 3], [1.0, 0.0, 0.0], 100.0) - 0.5 * central).abs() < 1e-6 * central, "half way in");
        assert!(c.along([0.0; 3], [0.0, 1.0, 0.0], 1000.0) == 0.0, "beside");
        assert!((erf(0.5) - 0.520_499_9).abs() < 1e-6 && (erf(-2.0) + 0.995_322_3).abs() < 1e-6);
    }

    /// In the Sun's neighbourhood the dust dims by about a magnitude per kiloparsec in the
    /// plane, patchily, and little toward the poles.
    #[test]
    fn dust_near_the_sun() {
        let s = Sectors::new(&Galaxy::milky_way(), Seed(4));
        let sun = [8200.0, 0.0, 25.0];
        let view = DustView::new(&s, sun, 3000.0);
        assert!(view.cloud_count() > 100, "{} clouds", view.cloud_count());
        let mut plane: Vec<f64> = (0..360)
            .map(|i| {
                let a = f64::from(i).to_radians();
                view.to([sun[0] + 1000.0 * a.cos(), sun[1] + 1000.0 * a.sin(), sun[2]])
            })
            .collect();
        plane.sort_by(f64::total_cmp);
        let (median, mean) = (plane[180], plane.iter().sum::<f64>() / 360.0);
        assert!((0.4..2.0).contains(&median) && mean > median, "median {median}, mean {mean} mag per kpc");
        let pole = view.to([sun[0], sun[1], sun[2] + 2000.0]);
        assert!(pole < 0.15, "{pole} toward the pole");
        // `along` agrees with `to`.
        let u = [0.6, 0.8, 0.0];
        let a = view.along(u, &[200.0, 700.0, 1500.0]);
        let direct = view.to([sun[0] + 1500.0 * 0.6, sun[1] + 1500.0 * 0.8, sun[2]]);
        assert!((a[2] - direct).abs() < 0.01 * direct + 1e-6, "{} vs {direct}", a[2]);
    }
}
