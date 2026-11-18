//! Plate motion as rigid rotation about Euler poles.
//!
//! v1 gave each plate one velocity (direction + speed) at its seed. Here each plate gets the
//! Euler pole that reproduces exactly that velocity at the seed point, so v1's motions carry
//! over, but velocity now varies correctly across the plate as a rotation on a sphere.

use geoforge_core::geom::UnitVec3;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

/// v1's random motion draw: per plate, a direction (degrees counter-clockwise from east) and
/// a speed in cm/yr, uniformly distributed.
#[must_use]
pub fn v1_random_motion(n: usize, seed: u64, min_cm_yr: f64, max_cm_yr: f64) -> Vec<(f64, f64)> {
    let mut rng = StdRng::seed_from_u64(seed);
    (0..n)
        .map(|_| {
            let dir = rng.gen_range(0.0..360.0);
            let speed = rng.gen_range(min_cm_yr..=max_cm_yr);
            (dir, speed)
        })
        .collect()
}

/// The Euler pole and rate (deg/Myr) whose rotation moves point `p` at `speed_cm_yr` in
/// direction `dir_deg` (counter-clockwise from east).
#[must_use]
pub fn euler_pole_for(p: UnitVec3, dir_deg: f64, speed_cm_yr: f64, radius_km: f64) -> (UnitVec3, f64) {
    let d = dir_deg.to_radians();
    let v = p.east() * d.cos() + p.north() * d.sin();
    let axis = p.v().cross(v).normalized().unwrap_or(UnitVec3::NORTH);
    // 1 cm/yr = 10 km/Myr.
    let rate_rad_per_myr = speed_cm_yr * 10.0 / radius_km;
    (axis, rate_rad_per_myr.to_degrees())
}

#[cfg(test)]
mod tests {
    use super::*;
    use geoforge_core::surface::Plate;

    #[test]
    fn pole_reproduces_seed_velocity() {
        let p = UnitVec3::from_lat_lon(20.0, 40.0);
        let (axis, rate) = euler_pole_for(p, 30.0, 5.0, 6371.0);
        let plate = Plate { name: "t".into(), euler_pole: axis, rate_deg_per_myr: rate, seed_point: p };
        let v = plate.velocity_km_per_myr(p, 6371.0);
        let speed_cm_yr = v.norm() / 10.0;
        assert!((speed_cm_yr - 5.0).abs() < 1e-9);
        let east = v.dot(p.east());
        let north = v.dot(p.north());
        assert!((north.atan2(east).to_degrees() - 30.0).abs() < 1e-9);
        // Earth-like magnitude: 5 cm/yr is about 0.45°/Myr.
        assert!((rate - 0.45).abs() < 0.01, "{rate}");
    }
}
