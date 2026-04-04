use thiserror::Error;

#[derive(Error, Debug)]
pub enum GeoforgeError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Image error: {0}")]
    #[cfg(feature = "export-png")]
    Image(#[from] image::ImageError),

    #[error("Dimensions mismatch: expected {expected}, got {actual}")]
    DimensionsMismatch { expected: String, actual: String },

    #[error("Invalid configuration: {0}")]
    Config(String),

    #[error("Simulation error: {0}")]
    Simulation(String),
}

pub type Result<T> = std::result::Result<T, GeoforgeError>;
