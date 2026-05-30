//! Error types for coordinate and parameter validation.

use thiserror::Error;

/// Errors raised when constructing or validating domain types.
#[derive(Debug, Error, PartialEq)]
pub enum TypesError {
    /// Map dimensions must be positive.
    #[error("map extent must have width and height > 0, got {width}x{height}")]
    InvalidMapExtent {
        /// Width in pixels.
        width: u32,
        /// Height in pixels.
        height: u32,
    },

    /// Latitude must lie in [-90, 90] degrees.
    #[error("latitude {lat}° is out of range [-90, 90]")]
    InvalidLatitude {
        /// Offending latitude in degrees.
        lat: f64,
    },

    /// Longitude must lie in [-180, 180] degrees when strict validation is used.
    #[error("longitude {lon}° is out of range [-180, 180]")]
    InvalidLongitude {
        /// Offending longitude in degrees.
        lon: f64,
    },

    /// Pixel coordinate outside map bounds.
    #[error("pixel ({x}, {y}) is outside map {width}x{height}")]
    PixelOutOfBounds {
        /// X pixel index.
        x: u32,
        /// Y pixel index.
        y: u32,
        /// Map width.
        width: u32,
        /// Map height.
        height: u32,
    },
}

pub(crate) type Result<T> = std::result::Result<T, TypesError>;
