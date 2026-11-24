//! Dark-matter halos: how many of each mass, and where, made cube by cube.
//!
//! The mass function is fitted to the standard cosmology's abundances (h = 0.7, Ω_m = 0.3;
//! Tinker et al. 2008): N(>M) = 8.9×10⁻⁴ Mpc⁻³ (M/1.43×10¹² M☉)^−0.9
//! exp(−(M/2.54×10¹⁴ M☉)^0.676), so ~9×10⁻⁴ Mpc⁻³ above a Milky Way's halo, ~7×10⁻⁶ above
//! 1.4×10¹⁴ M☉ (rich groups and clusters) and ~7×10⁻⁸ above 1.4×10¹⁵. Each half-decade of
//! mass from 10⁹ M☉ has its own grid of cubes, sized to hold a few halos: a cube's halos are
//! a pure function of the universe's seed and the cube's address.
//!
//! Halos cluster: they fall by a lognormal density field (Coles & Jones 1991), the exponential
//! of a Gaussian field of σ = 0.9 built from noise on scales of 4 to 32 Mpc, weighted by the
//! halo's bias: e^(bσg − (bσ)²/2), whose mean is 1, so the mass function holds on average.
//! The bias rises with mass, from 0.7 for dwarfs' halos to ~1 for a Milky Way's and ~3 for a
//! cluster's (Tinker et al. 2010). Within a cube, candidates are drawn at the field's ceiling
//! and kept in proportion to the field where they fall (thinning a Poisson process).

use crate::draw::Draw;
use geoforge_core::seed::Seed;

/// Lightest halo with a galaxy, M☉.
pub const LIGHTEST: f64 = 1e9;
/// Mass bins: half a decade each from [`LIGHTEST`].
pub const BINS: usize = 13;
/// The field's spread, in the Gaussian's units.
const SIGMA: f64 = 0.9;
/// Halos are drawn at the field's value this many σ up and thinned below it.
const CEILING_SIGMA: f64 = 3.0;
/// Smallest cube, Mpc.
const BASE_MPC: f64 = 0.25;

/// Halos above `m` M☉ per Mpc³.
#[must_use]
pub fn above(m: f64) -> f64 {
    8.9e-4 * (m / 1.43e12).powf(-0.9) * (-(m / 2.54e14).powf(0.676)).exp()
}

/// Mass range of bin `b`, M☉.
#[must_use]
pub fn bin_range(b: usize) -> (f64, f64) {
    let lo = LIGHTEST * 10f64.powf(0.5 * b as f64);
    (lo, lo * 10f64.sqrt())
}

/// Halos of bin `b` per Mpc³.
fn density(b: usize) -> f64 {
    let (lo, hi) = bin_range(b);
    above(lo) - above(hi)
}

/// Large-scale bias of a halo of `m` M☉.
#[must_use]
pub fn bias(m: f64) -> f64 {
    (0.7 + 0.35 * (m / 1e12).powf(0.35)).min(4.0)
}

/// Side of bin `b`'s cubes, Mpc: a power of two times [`BASE_MPC`], holding at least two
/// halos on average.
#[must_use]
pub fn cube_mpc(b: usize) -> f64 {
    let want = (2.0 / density(b)).cbrt();
    BASE_MPC * 2f64.powi((want / BASE_MPC).log2().ceil().max(0.0) as i32)
}

/// A halo.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Halo {
    pub bin: u8,
    pub cube: [i32; 3],
    pub n: u32,
    pub mass_msun: f64,
    pub position_mpc: [f64; 3],
}

/// The Gaussian field the halos follow: unit variance, structure from 4 to 32 Mpc.
#[derive(Debug, Clone)]
pub struct DensityField {
    seed: u64,
    norm: f64,
}

/// Octaves of the field: wavelength (Mpc) and weight.
const OCTAVES: [(f64, f64); 4] = [(32.0, 1.0), (16.0, 0.75), (8.0, 0.55), (4.0, 0.4)];

impl DensityField {
    #[must_use]
    pub fn new(seed: Seed) -> Self {
        let mut f = Self { seed: seed.child(b"density").0, norm: 1.0 };
        // Normalise to unit variance by sampling it.
        let mut d = Draw::new(seed.child(b"density norm"));
        let n = 4000;
        let (mut s, mut s2) = (0.0, 0.0);
        for _ in 0..n {
            let v = f.raw([d.range(-500.0, 500.0), d.range(-500.0, 500.0), d.range(-500.0, 500.0)]);
            s += v;
            s2 += v * v;
        }
        let mean = s / f64::from(n);
        f.norm = (s2 / f64::from(n) - mean * mean).sqrt().max(1e-9);
        f
    }

    fn raw(&self, p: [f64; 3]) -> f64 {
        OCTAVES.iter().enumerate().map(|(k, &(wave, w))| w * value_noise(self.seed.wrapping_add(k as u64 * 0x9E37_79B9_7F4A_7C15), p.map(|x| x / wave))).sum()
    }

    /// The Gaussian field at `p` (Mpc).
    #[must_use]
    pub fn gaussian(&self, p: [f64; 3]) -> f64 {
        self.raw(p) / self.norm
    }

    /// Relative density of halos of bias `b` at `p`: mean 1.
    #[must_use]
    pub fn weight(&self, p: [f64; 3], b: f64) -> f64 {
        let s = b * SIGMA;
        (s * self.gaussian(p) - 0.5 * s * s).exp()
    }
}

/// Value noise on the integer lattice, smoothly interpolated: −1 to 1.
fn value_noise(seed: u64, p: [f64; 3]) -> f64 {
    let i = p.map(|x| x.floor());
    let f = [0, 1, 2].map(|k| p[k] - i[k]);
    let fade = f.map(|t| t * t * t * (t * (t * 6.0 - 15.0) + 10.0));
    let at = |dx: i64, dy: i64, dz: i64| {
        let h = mix(seed ^ mix((i[0] as i64 + dx) as u64) ^ mix(((i[1] as i64 + dy) as u64).rotate_left(21)) ^ mix(((i[2] as i64 + dz) as u64).rotate_left(42)));
        (h >> 11) as f64 / (1u64 << 52) as f64 - 1.0
    };
    let lerp = |a: f64, b: f64, t: f64| a + (b - a) * t;
    let x00 = lerp(at(0, 0, 0), at(1, 0, 0), fade[0]);
    let x10 = lerp(at(0, 1, 0), at(1, 1, 0), fade[0]);
    let x01 = lerp(at(0, 0, 1), at(1, 0, 1), fade[0]);
    let x11 = lerp(at(0, 1, 1), at(1, 1, 1), fade[0]);
    lerp(lerp(x00, x10, fade[1]), lerp(x01, x11, fade[1]), fade[2])
}

/// SplitMix64's finaliser.
fn mix(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// The halos of bin `b` in cube `cube`.
#[must_use]
pub fn in_cube(seed: Seed, field: &DensityField, b: usize, cube: [i32; 3]) -> Vec<Halo> {
    let side = cube_mpc(b);
    let (lo, hi) = bin_range(b);
    let s = cube.iter().fold(seed.child_u64(b"halo bin", b as u64), |s, &v| s.child_u64(b"axis", u64::from(v as u32)));
    let mut d = Draw::new(s.child(b"halos"));
    // The heaviest of the bin sets the ceiling: its bias is the bin's highest.
    let top = bias(hi) * SIGMA;
    let ceiling = (top * CEILING_SIGMA - 0.5 * top * top).exp().max(1.0);
    let count = d.poisson(density(b) * side.powi(3) * ceiling);
    (0..count as u32)
        .filter_map(|n| {
            let p = [0, 1, 2].map(|a| (f64::from(cube[a]) + d.u()) * side);
            // dN/dlnM ∝ M^−0.9 within the bin.
            let u = d.u();
            let mass = lo * (1.0 - u * (1.0 - (lo / hi).powf(0.9))).powf(-1.0 / 0.9);
            let keep = d.u() * ceiling < field.weight(p, bias(mass));
            keep.then_some(Halo { bin: b as u8, cube, n, mass_msun: mass, position_mpc: p })
        })
        .collect()
}

/// The halos of bins `bins` within `radius` Mpc of `centre`.
#[must_use]
pub fn within(seed: Seed, field: &DensityField, centre: [f64; 3], radius: f64, bins: impl Iterator<Item = usize>) -> Vec<Halo> {
    let mut out = Vec::new();
    for b in bins {
        let side = cube_mpc(b);
        let lo = centre.map(|c| ((c - radius) / side).floor() as i32);
        let hi = centre.map(|c| ((c + radius) / side).floor() as i32);
        for i in lo[0]..=hi[0] {
            for j in lo[1]..=hi[1] {
                for k in lo[2]..=hi[2] {
                    out.extend(in_cube(seed, field, b, [i, j, k]).into_iter().filter(|h| dist(h.position_mpc, centre) <= radius));
                }
            }
        }
    }
    out
}

pub(crate) fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_mass_function_matches_the_abundances() {
        assert!((above(1.43e12) / 8.6e-4 - 1.0).abs() < 0.05);
        assert!((above(1.43e14) / 6.9e-6 - 1.0).abs() < 0.1);
        assert!((above(1.43e15) / 6.9e-8 - 1.0).abs() < 0.15);
        for b in 0..BINS {
            let n = density(b) * cube_mpc(b).powi(3);
            assert!((2.0..17.0).contains(&n), "bin {b}: {n} a cube");
        }
    }

    #[test]
    fn halos_are_as_many_as_the_mass_function_and_clustered() {
        let seed = Seed(9);
        let field = DensityField::new(seed);
        // Halos of 10¹¹–10¹¹·⁵ M☉ in an 80 Mpc box: the mean holds, and their counts in
        // 10 Mpc cells scatter beyond Poisson's by σ ≈ 0.7, like galaxies' σ₈ (~0.9 for L*
        // galaxies, ~0.7 for these fainter ones; 8 h⁻¹ Mpc spheres are 11 Mpc).
        let b = 4;
        let halos = within(seed, &field, [0.0; 3], 80.0, std::iter::once(b));
        let volume = 4.0 / 3.0 * std::f64::consts::PI * 80f64.powi(3);
        let ratio = halos.len() as f64 / (density(b) * volume);
        assert!((0.85..1.15).contains(&ratio), "{} halos, {ratio} of the expected", halos.len());
        let mut cells = std::collections::HashMap::new();
        for h in &halos {
            *cells.entry(h.position_mpc.map(|x| (x / 10.0).floor() as i32)).or_insert(0.0) += 1.0;
        }
        let inside: Vec<f64> = (-5..5).flat_map(|i| (-5..5).flat_map(move |j| (-5..5).map(move |k| [i, j, k]))).map(|c| *cells.get(&c).unwrap_or(&0.0)).collect();
        let mean = inside.iter().sum::<f64>() / inside.len() as f64;
        let var = inside.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / inside.len() as f64;
        let sigma = ((var - mean).max(0.0)).sqrt() / mean;
        assert!((0.5..1.0).contains(&sigma), "σ in 10 Mpc cells {sigma} (mean {mean})");
        // The same cube gives the same halos.
        assert_eq!(in_cube(seed, &field, b, [1, 2, 3]), in_cube(seed, &field, b, [1, 2, 3]));
    }
}
