//! Utility functions and helpers
//!
//! This module contains generic utility functions that aren't tied to specific
//! world generation logic.

pub mod color;

// Re-export commonly used utilities
pub use color::{hsv_to_rgb, rgb_to_hsv};
