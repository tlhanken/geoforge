//! A world seen whole: an equirectangular map wrapped on a lit sphere, as the viewer's globe
//! shows it.

use image::{Rgb, RgbImage};

/// `map` (equirectangular, north up, longitude −180° at its left edge) on a sphere seen from
/// above latitude `lat_deg`, longitude `lon_deg`, `size` pixels square on `background`, lit
/// from the upper left with a little ambient light; the limb is antialiased.
#[must_use]
pub fn render_globe(map: &RgbImage, size: u32, lat_deg: f64, lon_deg: f64, background: [u8; 3]) -> RgbImage {
    let (mw, mh) = (f64::from(map.width()), f64::from(map.height()));
    let r = f64::from(size) * 0.47;
    let c = f64::from(size) / 2.0;
    let (sl, cl) = lat_deg.to_radians().sin_cos();
    let light = normalize([-0.55, 0.45, 0.70]);
    RgbImage::from_fn(size, size, |px, py| {
        // Four samples a pixel, for a smooth limb.
        let mut acc = [0.0; 3];
        for (dx, dy) in [(0.25, 0.25), (0.75, 0.25), (0.25, 0.75), (0.75, 0.75)] {
            let x = (f64::from(px) + dx - c) / r;
            let y = (c - f64::from(py) - dy) / r;
            let rr = x * x + y * y;
            let col = if rr > 1.0 {
                background.map(f64::from)
            } else {
                let z = (1.0 - rr).sqrt();
                // Rotate the view's (x, y, z) back to the world: tilt by the latitude, then turn.
                let (wy, wz) = (y * cl + z * sl, -y * sl + z * cl);
                let lat = wy.asin();
                let lon = lon_deg.to_radians() + x.atan2(wz);
                let u = ((lon.to_degrees() + 180.0).rem_euclid(360.0) / 360.0 * mw).min(mw - 1.0);
                let v = ((90.0 - lat.to_degrees()) / 180.0 * mh).clamp(0.0, mh - 1.0);
                let p = map.get_pixel(u as u32, v as u32).0;
                let shade = 0.2 + 0.9 * (x * light[0] + y * light[1] + z * light[2]).max(0.0);
                p.map(|v| (f64::from(v) * shade).min(255.0))
            };
            for k in 0..3 {
                acc[k] += col[k] / 4.0;
            }
        }
        Rgb(acc.map(|v| v.round() as u8))
    })
}

fn normalize(v: [f64; 3]) -> [f64; 3] {
    let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    v.map(|x| x / n)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The map's point under the centre is the one asked for, the corners are background, and
    /// the lit side is brighter than the shaded one.
    #[test]
    fn a_globe_shows_the_map_where_asked() {
        // Left half red (west), right half blue (east).
        let map = RgbImage::from_fn(360, 180, |x, _| if x < 180 { Rgb([200, 0, 0]) } else { Rgb([0, 0, 200]) });
        let g = render_globe(&map, 101, 0.0, 90.0, [1, 2, 3]);
        let mid = g.get_pixel(50, 50).0;
        assert!(mid[2] > 100 && mid[0] == 0, "east under 90°E: {mid:?}");
        assert_eq!(g.get_pixel(0, 0).0, [1, 2, 3]);
        let west = render_globe(&map, 101, 0.0, -90.0, [0; 3]);
        assert!(west.get_pixel(50, 50).0[0] > 100);
        let (lit, dark) = (g.get_pixel(30, 40).0[2], g.get_pixel(80, 70).0[2]);
        assert!(lit > dark + 40, "lit {lit}, shaded {dark}");
    }
}
