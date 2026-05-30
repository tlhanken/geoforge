//! Surface coordinates: latitude/longitude and raster pixels.

use serde::{Deserialize, Serialize};

use crate::coordinates::{DistanceContext, DistanceKm};
use crate::error::{Result, TypesError};

/// Geographic position on a planetary surface (degrees).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct LatLon {
    /// Latitude in degrees, north positive.
    pub lat_deg: f64,
    /// Longitude in degrees, east positive.
    pub lon_deg: f64,
}

impl LatLon {
    /// Create a position, wrapping longitude to `[-180, 180]` and clamping latitude.
    #[must_use]
    pub fn new(lat_deg: f64, lon_deg: f64) -> Self {
        Self {
            lat_deg: lat_deg.clamp(-90.0, 90.0),
            lon_deg: normalize_longitude(lon_deg),
        }
    }

    /// Validate strict bounds without normalization.
    pub fn try_new(lat_deg: f64, lon_deg: f64) -> Result<Self> {
        if !(-90.0..=90.0).contains(&lat_deg) {
            return Err(TypesError::InvalidLatitude { lat: lat_deg });
        }
        if !(-180.0..=180.0).contains(&lon_deg) {
            return Err(TypesError::InvalidLongitude { lon: lon_deg });
        }
        Ok(Self {
            lat_deg,
            lon_deg,
        })
    }

    /// Convert to unit-sphere radians `(lat, lon)`.
    #[must_use]
    pub fn to_radians(self) -> (f64, f64) {
        (
            self.lat_deg.to_radians(),
            self.lon_deg.to_radians(),
        )
    }
}

impl DistanceKm for LatLon {
    fn distance_km(&self, other: &Self, ctx: &DistanceContext) -> f64 {
        let (lat1, lon1) = self.to_radians();
        let (lat2, lon2) = other.to_radians();
        haversine_rad(lat1, lon1, lat2, lon2) * ctx.radius_km
    }
}

/// Integer pixel index in a row-major raster.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct PixelCoord {
    /// Column (0 = western edge for global equirectangular).
    pub x: u32,
    /// Row (0 = northern edge).
    pub y: u32,
}

impl PixelCoord {
    /// Check bounds against a map extent.
    pub fn validate(self, extent: MapExtent) -> Result<Self> {
        if self.x >= extent.width || self.y >= extent.height {
            return Err(TypesError::PixelOutOfBounds {
                x: self.x,
                y: self.y,
                width: extent.width,
                height: extent.height,
            });
        }
        Ok(self)
    }
}

/// Geographic bounds and size of a planetary raster.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MapExtent {
    /// Width in pixels.
    pub width: u32,
    /// Height in pixels.
    pub height: u32,
    /// Western longitude bound (degrees).
    pub west_deg: f64,
    /// Eastern longitude bound (degrees).
    pub east_deg: f64,
    /// Northern latitude bound (degrees).
    pub north_deg: f64,
    /// Southern latitude bound (degrees).
    pub south_deg: f64,
}

impl MapExtent {
    /// Global equirectangular map (lon -180..180, lat 90..-90).
    pub fn try_global_equirectangular(width: u32, height: u32) -> Result<Self> {
        if width == 0 || height == 0 {
            return Err(TypesError::InvalidMapExtent { width, height });
        }
        Ok(Self {
            width,
            height,
            west_deg: -180.0,
            east_deg: 180.0,
            north_deg: 90.0,
            south_deg: -90.0,
        })
    }

    /// Longitude degrees per pixel.
    #[must_use]
    pub fn lon_resolution_deg(&self) -> f64 {
        (self.east_deg - self.west_deg) / f64::from(self.width)
    }

    /// Latitude degrees per pixel.
    #[must_use]
    pub fn lat_resolution_deg(&self) -> f64 {
        (self.north_deg - self.south_deg) / f64::from(self.height)
    }

    /// Center of pixel `(x, y)` as geographic coordinates.
    #[must_use]
    pub fn pixel_to_latlon(self, pixel: PixelCoord) -> LatLon {
        let lon = self.west_deg + (f64::from(pixel.x) + 0.5) * self.lon_resolution_deg();
        let lat = self.north_deg - (f64::from(pixel.y) + 0.5) * self.lat_resolution_deg();
        LatLon::new(lat, lon)
    }

    /// Nearest pixel for a geographic position (clamped to map edges).
    #[must_use]
    pub fn latlon_to_pixel(self, pos: LatLon) -> PixelCoord {
        let x = ((pos.lon_deg - self.west_deg) / self.lon_resolution_deg())
            .floor()
            .clamp(0.0, f64::from(self.width.saturating_sub(1))) as u32;
        let y = ((self.north_deg - pos.lat_deg) / self.lat_resolution_deg())
            .floor()
            .clamp(0.0, f64::from(self.height.saturating_sub(1))) as u32;
        PixelCoord { x, y }
    }

    /// Approximate kilometers per pixel at a given latitude (equirectangular).
    #[must_use]
    pub fn km_per_pixel_at_lat(&self, lat_deg: f64, radius_km: f64) -> f64 {
        let lat_rad = lat_deg.to_radians();
        let dx = self.lon_resolution_deg().to_radians() * radius_km * lat_rad.cos().abs();
        let dy = self.lat_resolution_deg().to_radians() * radius_km;
        (dx * dy).sqrt().max(0.01)
    }
}

/// Haversine great-circle distance in radians.
#[must_use]
pub fn haversine_rad(lat1: f64, lon1: f64, lat2: f64, lon2: f64) -> f64 {
    let dlat = lat2 - lat1;
    let dlon = lon2 - lon1;
    let a = (dlat / 2.0).sin().powi(2)
        + lat1.cos() * lat2.cos() * (dlon / 2.0).sin().powi(2);
    2.0 * a.sqrt().asin()
}

fn normalize_longitude(lon_deg: f64) -> f64 {
    let mut lon = lon_deg % 360.0;
    if lon > 180.0 {
        lon -= 360.0;
    } else if lon <= -180.0 {
        lon += 360.0;
    }
    lon
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn longitude_wraps() {
        let p = LatLon::new(0.0, 270.0);
        assert!((p.lon_deg - (-90.0)).abs() < 1e-10);
    }

    #[test]
    fn pixel_round_trip() {
        let extent = MapExtent::try_global_equirectangular(360, 180).unwrap();
        let pixel = PixelCoord { x: 180, y: 90 };
        let ll = extent.pixel_to_latlon(pixel);
        let back = extent.latlon_to_pixel(ll);
        assert_eq!(pixel, back);
    }

    #[test]
    fn haversine_equator_quarter_earth() {
        let ctx = DistanceContext::earth();
        let a = LatLon::new(0.0, 0.0);
        let b = LatLon::new(0.0, 90.0);
        let d = a.distance_km(&b, &ctx);
        // quarter circumference ~ 10_000 km
        assert!((d - 10_001.966).abs() < 50.0);
    }
}
