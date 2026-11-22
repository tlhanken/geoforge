//! Positions on the unit sphere.
//!
//! Simulation uses [`UnitVec3`] (z = north pole, x = longitude 0); [`LatLon`] exists only for
//! files and display. The conversions match geoforge v1 bit for bit so ported algorithms
//! reproduce v1 output.

use serde::{Deserialize, Serialize};
use std::ops::{Add, Mul, Neg, Sub};

/// A general 3-vector (velocities, angular velocities, offsets).
#[derive(Debug, Clone, Copy, PartialEq, Default, Serialize, Deserialize)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3 {
    #[must_use]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }
    #[must_use]
    pub fn dot(self, o: Self) -> f64 {
        self.x * o.x + self.y * o.y + self.z * o.z
    }
    #[must_use]
    pub fn cross(self, o: Self) -> Self {
        Self::new(self.y * o.z - self.z * o.y, self.z * o.x - self.x * o.z, self.x * o.y - self.y * o.x)
    }
    #[must_use]
    pub fn norm(self) -> f64 {
        self.dot(self).sqrt()
    }
    /// Unit vector in this direction, or `None` for (near) zero.
    #[must_use]
    pub fn normalized(self) -> Option<UnitVec3> {
        let n = self.norm();
        (n > 1e-15).then(|| UnitVec3 { x: self.x / n, y: self.y / n, z: self.z / n })
    }
}

impl Add for Vec3 {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self::new(self.x + o.x, self.y + o.y, self.z + o.z)
    }
}
impl Sub for Vec3 {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Self::new(self.x - o.x, self.y - o.y, self.z - o.z)
    }
}
impl Mul<f64> for Vec3 {
    type Output = Self;
    fn mul(self, s: f64) -> Self {
        Self::new(self.x * s, self.y * s, self.z * s)
    }
}
impl Neg for Vec3 {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(-self.x, -self.y, -self.z)
    }
}

/// A point on the unit sphere.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct UnitVec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

/// Latitude and longitude in degrees (files and display only).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct LatLon {
    pub lat_deg: f64,
    pub lon_deg: f64,
}

impl LatLon {
    #[must_use]
    pub const fn new(lat_deg: f64, lon_deg: f64) -> Self {
        Self { lat_deg, lon_deg }
    }
}

impl UnitVec3 {
    pub const NORTH: Self = Self { x: 0.0, y: 0.0, z: 1.0 };

    #[must_use]
    pub fn from_lat_lon(lat_deg: f64, lon_deg: f64) -> Self {
        let lat = lat_deg.to_radians();
        let lon = lon_deg.to_radians();
        Self { x: lat.cos() * lon.cos(), y: lat.cos() * lon.sin(), z: lat.sin() }
    }

    #[must_use]
    pub fn to_lat_lon(self) -> LatLon {
        LatLon { lat_deg: self.z.asin().to_degrees(), lon_deg: self.y.atan2(self.x).to_degrees() }
    }

    #[must_use]
    pub fn v(self) -> Vec3 {
        Vec3::new(self.x, self.y, self.z)
    }

    #[must_use]
    pub fn dot(self, o: Self) -> f64 {
        self.x * o.x + self.y * o.y + self.z * o.z
    }

    /// Great-circle angle to another point, radians.
    #[must_use]
    pub fn angle_to(self, o: Self) -> f64 {
        self.dot(o).clamp(-1.0, 1.0).acos()
    }

    /// Remove the radial component of `v` at this point.
    #[must_use]
    pub fn tangent(self, v: Vec3) -> Vec3 {
        v - self.v() * v.dot(self.v())
    }

    /// Local unit east vector (undefined exactly at the poles; returns +y there).
    #[must_use]
    pub fn east(self) -> Vec3 {
        Vec3::new(0.0, 0.0, 1.0).cross(self.v()).normalized().map_or(Vec3::new(0.0, 1.0, 0.0), UnitVec3::v)
    }

    /// Local unit north vector.
    #[must_use]
    pub fn north(self) -> Vec3 {
        self.v().cross(self.east())
    }

    /// Rotate this point about `axis` by `angle` radians (right-hand rule).
    #[must_use]
    pub fn rotate_about(self, axis: UnitVec3, angle: f64) -> UnitVec3 {
        let k = axis.v();
        let v = self.v();
        let (s, c) = angle.sin_cos();
        let r = v * c + k.cross(v) * s + k * (k.dot(v) * (1.0 - c));
        r.normalized().unwrap_or(self)
    }

    /// Move along the great circle with initial `bearing_deg` (clockwise from north) by
    /// `angle` radians.
    #[must_use]
    pub fn travel(self, bearing_deg: f64, angle: f64) -> UnitVec3 {
        let b = bearing_deg.to_radians();
        let dir = self.north() * b.cos() + self.east() * b.sin();
        let r = self.v() * angle.cos() + dir * angle.sin();
        r.normalized().unwrap_or(self)
    }

    /// Initial bearing to `to`, degrees clockwise from north in `[0, 360)` (v1 formula).
    #[must_use]
    pub fn bearing_to(self, to: UnitVec3) -> f64 {
        let a = self.to_lat_lon();
        let b = to.to_lat_lon();
        let lat1 = a.lat_deg.to_radians();
        let lat2 = b.lat_deg.to_radians();
        let dlon = (b.lon_deg - a.lon_deg).to_radians();
        let y = dlon.sin() * lat2.cos();
        let x = lat1.cos() * lat2.sin() - lat1.sin() * lat2.cos() * dlon.cos();
        (y.atan2(x).to_degrees() + 360.0) % 360.0
    }

    /// Spherical centroid of a set of points (normalized mean), if well defined.
    #[must_use]
    pub fn centroid(points: impl IntoIterator<Item = UnitVec3>) -> Option<UnitVec3> {
        let mut s = Vec3::default();
        for p in points {
            s = s + p.v();
        }
        s.normalized()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn lat_lon_roundtrip() {
        let p = UnitVec3::from_lat_lon(45.0, -120.0);
        let ll = p.to_lat_lon();
        assert!((ll.lat_deg - 45.0).abs() < 1e-10);
        assert!((ll.lon_deg + 120.0).abs() < 1e-10);
    }

    #[test]
    fn angles_and_rotation() {
        let a = UnitVec3::from_lat_lon(0.0, 0.0);
        let b = UnitVec3::from_lat_lon(0.0, 90.0);
        assert!((a.angle_to(b) - PI / 2.0).abs() < 1e-12);
        let r = a.rotate_about(UnitVec3::NORTH, PI / 2.0);
        assert!(r.angle_to(b) < 1e-12);
    }

    #[test]
    fn local_basis_and_travel() {
        let p = UnitVec3::from_lat_lon(30.0, 10.0);
        assert!(p.east().dot(p.v()).abs() < 1e-12);
        assert!(p.north().dot(p.v()).abs() < 1e-12);
        assert!(p.north().z > 0.0);
        let q = p.travel(0.0, 10f64.to_radians());
        assert!((q.to_lat_lon().lat_deg - 40.0).abs() < 1e-9);
        let e = p.travel(90.0, 1e-3);
        assert!((p.bearing_to(e) - 90.0).abs() < 0.01);
    }
}
