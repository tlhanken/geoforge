//! Color conversion utilities
//!
//! Provides conversions between RGB and HSV color spaces.

/// Convert HSV to RGB
///
/// # Arguments
/// * `h` - Hue in degrees (0-360)
/// * `s` - Saturation (0.0-1.0, representing 0%-100%)
/// * `v` - Value/Brightness (0.0-1.0, representing 0%-100%)
///
/// # Returns
/// RGB values as [r, g, b] where each component is 0-255
///
/// # Example
/// ```
/// use geoforge::utils::color::hsv_to_rgb;
///
/// // Pure red
/// let rgb = hsv_to_rgb(0.0, 1.0, 1.0);
/// assert_eq!(rgb, [255, 0, 0]);
///
/// // Pure green
/// let rgb = hsv_to_rgb(120.0, 1.0, 1.0);
/// assert_eq!(rgb, [0, 255, 0]);
///
/// // Pure blue
/// let rgb = hsv_to_rgb(240.0, 1.0, 1.0);
/// assert_eq!(rgb, [0, 0, 255]);
/// ```
pub fn hsv_to_rgb(h: f64, s: f64, v: f64) -> [u8; 3] {
    let c = v * s;
    let h_prime = h / 60.0;
    let x = c * (1.0 - ((h_prime % 2.0) - 1.0).abs());
    let m = v - c;

    let (r, g, b) = if h_prime < 1.0 {
        (c, x, 0.0)
    } else if h_prime < 2.0 {
        (x, c, 0.0)
    } else if h_prime < 3.0 {
        (0.0, c, x)
    } else if h_prime < 4.0 {
        (0.0, x, c)
    } else if h_prime < 5.0 {
        (x, 0.0, c)
    } else {
        (c, 0.0, x)
    };

    [
        ((r + m) * 255.0).round() as u8,
        ((g + m) * 255.0).round() as u8,
        ((b + m) * 255.0).round() as u8,
    ]
}

/// Convert RGB to HSV
///
/// # Arguments
/// * `r` - Red component (0-255)
/// * `g` - Green component (0-255)
/// * `b` - Blue component (0-255)
///
/// # Returns
/// Tuple of (hue, saturation, value) where:
/// * Hue: 0-360 degrees
/// * Saturation: 0.0-1.0 (0%-100%)
/// * Value: 0.0-1.0 (0%-100%)
///
/// # Example
/// ```
/// use geoforge::utils::color::rgb_to_hsv;
///
/// // Pure red
/// let (h, s, v) = rgb_to_hsv(255, 0, 0);
/// assert!((h - 0.0).abs() < 0.1);
/// assert!((s - 1.0).abs() < 0.01);
/// assert!((v - 1.0).abs() < 0.01);
///
/// // Gray (low saturation)
/// let (_, s, _) = rgb_to_hsv(128, 128, 128);
/// assert!((s - 0.0).abs() < 0.01);
/// ```
pub fn rgb_to_hsv(r: u8, g: u8, b: u8) -> (f64, f64, f64) {
    let r = r as f64 / 255.0;
    let g = g as f64 / 255.0;
    let b = b as f64 / 255.0;

    let max = r.max(g).max(b);
    let min = r.min(g).min(b);
    let delta = max - min;

    // Value
    let v = max;

    // Saturation
    let s = if max == 0.0 {
        0.0
    } else {
        delta / max
    };

    // Hue
    let h = if delta == 0.0 {
        0.0
    } else if max == r {
        60.0 * (((g - b) / delta) % 6.0)
    } else if max == g {
        60.0 * (((b - r) / delta) + 2.0)
    } else {
        60.0 * (((r - g) / delta) + 4.0)
    };

    // Normalize hue to 0-360
    let h = if h < 0.0 { h + 360.0 } else { h };

    (h, s, v)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hsv_to_rgb_primary_colors() {
        // Red
        assert_eq!(hsv_to_rgb(0.0, 1.0, 1.0), [255, 0, 0]);

        // Green
        assert_eq!(hsv_to_rgb(120.0, 1.0, 1.0), [0, 255, 0]);

        // Blue
        assert_eq!(hsv_to_rgb(240.0, 1.0, 1.0), [0, 0, 255]);
    }

    #[test]
    fn test_hsv_to_rgb_grayscale() {
        // Black (value = 0)
        assert_eq!(hsv_to_rgb(0.0, 0.0, 0.0), [0, 0, 0]);

        // White (saturation = 0, value = 1)
        assert_eq!(hsv_to_rgb(0.0, 0.0, 1.0), [255, 255, 255]);

        // Gray (saturation = 0, value = 0.5)
        assert_eq!(hsv_to_rgb(0.0, 0.0, 0.5), [128, 128, 128]);
    }

    #[test]
    fn test_rgb_to_hsv_primary_colors() {
        // Red
        let (h, s, v) = rgb_to_hsv(255, 0, 0);
        assert!((h - 0.0).abs() < 0.1);
        assert!((s - 1.0).abs() < 0.01);
        assert!((v - 1.0).abs() < 0.01);

        // Green
        let (h, s, v) = rgb_to_hsv(0, 255, 0);
        assert!((h - 120.0).abs() < 0.1);
        assert!((s - 1.0).abs() < 0.01);
        assert!((v - 1.0).abs() < 0.01);

        // Blue
        let (h, s, v) = rgb_to_hsv(0, 0, 255);
        assert!((h - 240.0).abs() < 0.1);
        assert!((s - 1.0).abs() < 0.01);
        assert!((v - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_rgb_to_hsv_grayscale() {
        // Black
        let (_, s, v) = rgb_to_hsv(0, 0, 0);
        assert!((s - 0.0).abs() < 0.01);
        assert!((v - 0.0).abs() < 0.01);

        // White
        let (_, s, v) = rgb_to_hsv(255, 255, 255);
        assert!((s - 0.0).abs() < 0.01);
        assert!((v - 1.0).abs() < 0.01);

        // Gray
        let (_, s, v) = rgb_to_hsv(128, 128, 128);
        assert!((s - 0.0).abs() < 0.01);
        assert!((v - 0.502).abs() < 0.01); // 128/255 ≈ 0.502
    }

    #[test]
    fn test_roundtrip_conversion() {
        // Test that HSV -> RGB -> HSV preserves values
        let original_hsv = (45.0, 0.75, 0.85);
        let rgb = hsv_to_rgb(original_hsv.0, original_hsv.1, original_hsv.2);
        let (h, s, v) = rgb_to_hsv(rgb[0], rgb[1], rgb[2]);

        assert!((h - original_hsv.0).abs() < 1.0); // Hue within 1 degree
        assert!((s - original_hsv.1).abs() < 0.02); // Saturation within 2%
        assert!((v - original_hsv.2).abs() < 0.02); // Value within 2%
    }
}
