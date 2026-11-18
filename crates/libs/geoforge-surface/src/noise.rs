//! Deterministic 3D value noise.
//!
//! [`ValueNoise`] is a line-for-line port of geoforge v1's boundary-refinement noise so that
//! the ported domain warp reproduces v1 plate maps exactly.

use geoforge_core::geom::UnitVec3;

/// v1's hashed trilinear value noise with smoothstep, summed over octaves.
#[derive(Debug, Clone, Copy)]
pub struct ValueNoise {
    pub seed: u64,
    pub octaves: usize,
    pub persistence: f64,
}

impl ValueNoise {
    #[must_use]
    pub fn new(seed: u64, octaves: usize, persistence: f64) -> Self {
        Self { seed, octaves, persistence }
    }

    /// Multi-octave noise in roughly `[-1, 1]`.
    #[must_use]
    pub fn fbm(&self, x: f64, y: f64, z: f64, channel: u64) -> f64 {
        let mut total = 0.0;
        let mut amplitude = 1.0;
        let mut frequency = 1.0;
        let mut max_value = 0.0;
        for _ in 0..self.octaves {
            total += self.noise(x * frequency, y * frequency, z * frequency, channel) * amplitude;
            max_value += amplitude;
            amplitude *= self.persistence;
            frequency *= 2.0;
        }
        total / max_value
    }

    /// Noise sampled on the sphere at `p` with spatial `scale` (features per radian, roughly).
    #[must_use]
    pub fn at(&self, p: UnitVec3, scale: f64, channel: u64) -> f64 {
        self.fbm(p.x * scale, p.y * scale, p.z * scale, channel)
    }

    fn noise(&self, x: f64, y: f64, z: f64, channel: u64) -> f64 {
        let x0 = x.floor() as i64;
        let y0 = y.floor() as i64;
        let z0 = z.floor() as i64;
        let x1 = x0 + 1;
        let y1 = y0 + 1;
        let z1 = z0 + 1;
        let fx = x - x0 as f64;
        let fy = y - y0 as f64;
        let fz = z - z0 as f64;
        let sx = smoothstep(fx);
        let sy = smoothstep(fy);
        let sz = smoothstep(fz);
        let n000 = self.hash(x0, y0, z0, channel);
        let n100 = self.hash(x1, y0, z0, channel);
        let n010 = self.hash(x0, y1, z0, channel);
        let n110 = self.hash(x1, y1, z0, channel);
        let n001 = self.hash(x0, y0, z1, channel);
        let n101 = self.hash(x1, y0, z1, channel);
        let n011 = self.hash(x0, y1, z1, channel);
        let n111 = self.hash(x1, y1, z1, channel);
        let nx00 = lerp(n000, n100, sx);
        let nx10 = lerp(n010, n110, sx);
        let nx01 = lerp(n001, n101, sx);
        let nx11 = lerp(n011, n111, sx);
        let ny0 = lerp(nx00, nx10, sy);
        let ny1 = lerp(nx01, nx11, sy);
        lerp(ny0, ny1, sz)
    }

    fn hash(&self, x: i64, y: i64, z: i64, channel: u64) -> f64 {
        let mut hash = self.seed.wrapping_add(channel.wrapping_mul(1_000_000_007)).wrapping_mul(2_654_435_761);
        hash = hash.wrapping_add(x as u64).wrapping_mul(2_654_435_761);
        hash = hash.wrapping_add(y as u64).wrapping_mul(2_654_435_761);
        hash = hash.wrapping_add(z as u64).wrapping_mul(2_654_435_761);
        (hash as f64 / u64::MAX as f64) * 2.0 - 1.0
    }
}

fn smoothstep(t: f64) -> f64 {
    t * t * (3.0 - 2.0 * t)
}

fn lerp(a: f64, b: f64, t: f64) -> f64 {
    a + (b - a) * t
}
