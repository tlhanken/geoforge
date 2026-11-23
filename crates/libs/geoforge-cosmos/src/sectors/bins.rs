//! Brightness bins and the cubes each lives in.
//!
//! Stars are binned by absolute visual magnitude, one magnitude a bin from −12 to +18 (the
//! faintest bin takes everything fainter). A bin's reach is how far its brightest stars stay
//! brighter than a limit with no dust; its cubes are about half that reach on a side, in
//! powers of two of 8 pc, so a query for a bin around a point visits about 5³ cubes.

/// The brightest bin's bright edge, absolute V.
pub const BRIGHTEST: f64 = -12.0;
/// Number of bins.
pub const BINS: usize = 30;
/// The smallest cube's side, pc.
pub const BASE_PC: f64 = 8.0;
/// The largest cube's side: 8 pc × 2¹¹ ≈ 16 kpc.
const MAX_LEVEL: u32 = 11;
/// The limit cube sizes are set for: the naked eye.
const SIZING_LIMIT: f64 = 6.5;

/// The bin of absolute magnitude `v`.
#[must_use]
pub fn bin_of(v: f64) -> usize {
    ((v - BRIGHTEST).floor().max(0.0) as usize).min(BINS - 1)
}

/// The bright edge of bin `b`, absolute V.
#[must_use]
pub fn bright_edge(b: usize) -> f64 {
    BRIGHTEST + b as f64
}

/// How far a star of absolute magnitude `v` stays brighter than `limit` without dust, pc.
#[must_use]
pub fn reach_pc(v: f64, limit: f64) -> f64 {
    10f64.powf((limit - v) / 5.0 + 1.0)
}

/// The level of bin `b`'s cubes.
#[must_use]
pub fn level(b: usize) -> u32 {
    let half = reach_pc(bright_edge(b), SIZING_LIMIT) / 2.0;
    ((half / BASE_PC).log2().round().max(0.0) as u32).min(MAX_LEVEL)
}

/// The side of cubes at `level`, pc.
#[must_use]
pub fn side_pc(level: u32) -> f64 {
    BASE_PC * f64::from(1u32 << level)
}

/// The cube at `level` holding `p`.
#[must_use]
pub fn cube_of(p: [f64; 3], level: u32) -> [i32; 3] {
    let s = side_pc(level);
    p.map(|v| (v / s).floor() as i32)
}

/// The corners of cube `c` at `level`.
#[must_use]
pub fn corners(c: [i32; 3], level: u32) -> ([f64; 3], [f64; 3]) {
    let s = side_pc(level);
    (c.map(|v| f64::from(v) * s), c.map(|v| f64::from(v + 1) * s))
}

/// Distance from `p` to the nearest point of cube `c` at `level`, pc (0 inside).
#[must_use]
pub fn distance_to(p: [f64; 3], c: [i32; 3], level: u32) -> f64 {
    let (a, b) = corners(c, level);
    (0..3).map(|k| (a[k] - p[k]).max(0.0).max(p[k] - b[k]).powi(2)).sum::<f64>().sqrt()
}

/// The cubes at `level` that reach within `radius` pc of `p`.
pub fn cubes_near(p: [f64; 3], radius: f64, level: u32) -> impl Iterator<Item = [i32; 3]> {
    let s = side_pc(level);
    let lo = p.map(|v| ((v - radius) / s).floor() as i32);
    let hi = p.map(|v| ((v + radius) / s).floor() as i32);
    (lo[0]..=hi[0]).flat_map(move |i| (lo[1]..=hi[1]).flat_map(move |j| (lo[2]..=hi[2]).map(move |k| [i, j, k]))).filter(move |&c| distance_to(p, c, level) <= radius)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bins_and_cubes() {
        assert_eq!(bin_of(4.83), 16);
        assert_eq!(bin_of(-40.0), 0);
        assert_eq!(bin_of(25.0), BINS - 1);
        // The Sun is seen to 22 pc, a supergiant of −7 to 5 kpc.
        assert!((reach_pc(4.83, 6.5) - 21.6).abs() < 0.2);
        assert!((reach_pc(-7.0, 6.5) / 5012.0 - 1.0).abs() < 0.01);
        // Cubes grow with reach and are nested.
        assert!(level(bin_of(-7.0)) > level(bin_of(1.0)) && level(bin_of(1.0)) > level(bin_of(10.0)));
        assert_eq!(level(BINS - 1), 0);
        let p = [8200.5, -3.2, 25.0];
        let (a, b) = corners(cube_of(p, 4), 4);
        assert!((0..3).all(|k| a[k] <= p[k] && p[k] < b[k]));
        assert_eq!(cube_of(p, 5), cube_of(p, 4).map(|v| v.div_euclid(2)));
        // Every point within the radius lies in a listed cube.
        let near: Vec<[i32; 3]> = cubes_near(p, 30.0, 1).collect();
        for q in [[8230.0, -3.2, 25.0], [8200.5, 26.0, 25.0], [8200.5, -3.2, -4.0]] {
            assert!(near.contains(&cube_of(q, 1)), "{q:?}");
        }
        assert!(near.len() < 150);
    }
}
