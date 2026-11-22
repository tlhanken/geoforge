//! Polar cap views: the regions poleward of 60° in an azimuthal equidistant projection, for
//! painting the poles in flat image editors without the equirectangular view's stretching.
//!
//! A cap is a square image with the pole at its centre and 60° on the inscribed circle;
//! distance from the centre is distance from the pole. The north cap has 0° longitude
//! straight down and east anticlockwise, as seen from above the pole; the south cap has 0°
//! straight up and east clockwise. A world file and projection definition come with each,
//! so GIS software places it.
//!
//! On import a cap is compared with the cap the export made: every grid cell whose nearest
//! cap pixel was painted takes that pixel, over whatever the main view says there.

use crate::edit::BodyRef;
use crate::raster::{each, Pixels};
use crate::{Result, StoreError};
use std::path::Path;

/// Caps reach from the pole to this latitude, degrees.
pub const CAP_EDGE_DEG: f64 = 60.0;

/// A polar cap view `side` pixels square.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Cap {
    pub north: bool,
    pub side: u32,
}

impl Cap {
    /// The cap for a grid `height` rows tall, at `scale`: twice as many pixels from pole to
    /// edge as the grid has rows over that span, so the edge circle has about as many pixels
    /// as a grid row there.
    #[must_use]
    pub fn for_grid(north: bool, height: u32, scale: u32) -> Self {
        let rows = f64::from(height) * (90.0 - CAP_EDGE_DEG) / 180.0;
        let side = (4.0 * rows * f64::from(scale.max(1))).round() as u32;
        Self { north, side: side.max(2) }
    }

    /// "north" or "south", as in the file names.
    #[must_use]
    pub fn name(self) -> &'static str {
        if self.north {
            "north"
        } else {
            "south"
        }
    }

    fn radius(self) -> f64 {
        f64::from(self.side) / 2.0
    }

    /// The pixel (column, row; fractional, from the top-left corner) of a point, if on the cap.
    #[must_use]
    pub fn pixel_of(self, lat: f64, lon: f64) -> Option<(f64, f64)> {
        let from_pole = 90.0 - if self.north { lat } else { -lat };
        if from_pole > 90.0 - CAP_EDGE_DEG {
            return None;
        }
        let rho = from_pole / (90.0 - CAP_EDGE_DEG) * self.radius();
        let (s, c) = lon.to_radians().sin_cos();
        let down = if self.north { c } else { -c };
        Some((self.radius() + rho * s, self.radius() + rho * down))
    }

    /// The latitude and longitude at a pixel position, if inside the cap.
    #[must_use]
    pub fn lat_lon_of(self, x: f64, y: f64) -> Option<(f64, f64)> {
        let (dx, dy) = (x - self.radius(), y - self.radius());
        let rho = dx.hypot(dy);
        if rho > self.radius() {
            return None;
        }
        let from_pole = rho / self.radius() * (90.0 - CAP_EDGE_DEG);
        let down = if self.north { dy } else { -dy };
        let lon = dx.atan2(down).to_degrees();
        Some((if self.north { 90.0 - from_pole } else { from_pole - 90.0 }, lon))
    }

    /// The cap of an equirectangular view `width` × `height`: each pixel the nearest cell's;
    /// zero (black) outside the circle.
    #[must_use]
    pub fn render(self, view: &Pixels, width: u32, height: u32) -> Pixels {
        let n = self.side as usize;
        let at = |i: usize| self.lat_lon_of(i as f64 % n as f64 + 0.5, (i / n) as f64 + 0.5).map(|(lat, lon)| cell_of(lat, lon, width, height));
        each!(view, |v| (0..n * n).map(|i| at(i).map_or_else(Default::default, |c| v[c])).collect())
    }

    /// Fold an edited cap into an equirectangular view `width` × `height`: each cell poleward
    /// of the edge whose nearest pixel in `edited` differs from `exported` takes that pixel.
    /// Returns how many cells changed.
    ///
    /// # Errors
    /// If the three are not the same kind of pixels or the caps' sizes differ.
    pub fn merge(self, edited: &Pixels, exported: &Pixels, into: &mut Pixels, width: u32, height: u32) -> Result<usize> {
        let n = self.side as usize;
        if edited.len() != n * n || exported.len() != n * n {
            return Err(StoreError::Invalid(format!("the {} cap is not {n}x{n}", self.name())));
        }
        let cells: Vec<(usize, usize)> = (0..(width * height) as usize)
            .filter_map(|c| {
                let (lat, lon) = centre(c, width, height);
                let (x, y) = self.pixel_of(lat, lon)?;
                let (x, y) = ((x as usize).min(n - 1), (y as usize).min(n - 1));
                Some((c, y * n + x))
            })
            .collect();
        match (edited, exported, into) {
            (Pixels::Gray16(e), Pixels::Gray16(x), Pixels::Gray16(v)) => Ok(fold_changes(&cells, e, x, v)),
            (Pixels::Rgb16(e), Pixels::Rgb16(x), Pixels::Rgb16(v)) => Ok(fold_changes(&cells, e, x, v)),
            (Pixels::Rgb8(e), Pixels::Rgb8(x), Pixels::Rgb8(v)) => Ok(fold_changes(&cells, e, x, v)),
            _ => Err(StoreError::Invalid(format!("the {} cap's pixels are not the view's kind", self.name()))),
        }
    }

    /// Write `<stem>.<north|south>.pgw` and `.prj`: an azimuthal equidistant projection
    /// about the pole, in metres on the body's sphere.
    ///
    /// # Errors
    /// On IO errors.
    pub fn write_georef(self, dir: &Path, stem: &str, body: &BodyRef) -> Result<()> {
        let metres = (90.0 - CAP_EDGE_DEG).to_radians() * body.radius_km * 1000.0 / self.radius();
        let corner = (self.radius() - 0.5) * metres;
        let pgw = format!("{metres:.6}\n0.0\n0.0\n{:.6}\n{:.6}\n{corner:.6}\n", -metres, -corner);
        let file = format!("{stem}.{}", self.name());
        std::fs::write(dir.join(format!("{file}.pgw")), pgw)?;
        let n = body.name.replace('"', "");
        let prj = format!(
            "PROJCS[\"{n} {} polar azimuthal equidistant\",GEOGCS[\"{n} sphere\",DATUM[\"D_{n}\",SPHEROID[\"{n}\",{:.1},0]],PRIMEM[\"Reference_Meridian\",0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Azimuthal_Equidistant\"],PARAMETER[\"latitude_of_center\",{}],PARAMETER[\"longitude_of_center\",0],PARAMETER[\"false_easting\",0],PARAMETER[\"false_northing\",0],UNIT[\"Meter\",1]]",
            self.name(),
            body.radius_km * 1000.0,
            if self.north { 90 } else { -90 }
        );
        std::fs::write(dir.join(format!("{file}.prj")), prj)?;
        Ok(())
    }
}

fn fold_changes<T: Copy + PartialEq>(cells: &[(usize, usize)], edited: &[T], exported: &[T], into: &mut [T]) -> usize {
    let mut n = 0;
    for &(c, px) in cells {
        if edited[px] != exported[px] {
            into[c] = edited[px];
            n += 1;
        }
    }
    n
}

/// The centre of cell `c` of an equirectangular grid, (lat, lon) degrees.
fn centre(c: usize, width: u32, height: u32) -> (f64, f64) {
    let (w, h) = (width as usize, height as usize);
    (90.0 - (((c / w) as f64) + 0.5) * 180.0 / h as f64, -180.0 + (((c % w) as f64) + 0.5) * 360.0 / w as f64)
}

/// The cell holding (lat, lon).
fn cell_of(lat: f64, lon: f64, width: u32, height: u32) -> usize {
    let row = (((90.0 - lat) / 180.0 * f64::from(height)) as usize).min(height as usize - 1);
    let col = (((lon + 180.0) / 360.0 * f64::from(width)).floor() as i64).rem_euclid(i64::from(width)) as usize;
    row * width as usize + col
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn projection_round_trips() {
        for north in [true, false] {
            let cap = Cap { north, side: 400 };
            for (lat, lon) in [(75.0, 0.0), (61.0, 135.0), (89.0, -90.0), (70.0, -179.0)] {
                let lat = if north { lat } else { -lat };
                let (x, y) = cap.pixel_of(lat, lon).unwrap();
                let (la, lo) = cap.lat_lon_of(x, y).unwrap();
                assert!((la - lat).abs() < 1e-9 && ((lo - lon + 540.0) % 360.0 - 180.0).abs() < 1e-9, "{lat},{lon} -> {la},{lo}");
            }
            assert!(cap.pixel_of(if north { 50.0 } else { -50.0 }, 0.0).is_none());
        }
        // North: 0° straight down, 90°E to the right; south: 0° straight up.
        let (n, s) = (Cap { north: true, side: 100 }, Cap { north: false, side: 100 });
        assert!(n.pixel_of(75.0, 0.0).unwrap().1 > 50.0 && n.pixel_of(75.0, 90.0).unwrap().0 > 50.0);
        assert!(s.pixel_of(-75.0, 0.0).unwrap().1 < 50.0 && s.pixel_of(-75.0, 90.0).unwrap().0 > 50.0);
    }

    /// An untouched cap changes nothing; a painted disc around the pole changes only cells
    /// near the pole, to the paint.
    #[test]
    fn painting_a_cap() {
        let (w, h) = (64u32, 32u32);
        let view = Pixels::Gray16((0..w * h).map(|i| (i % 977) as u16).collect());
        let cap = Cap::for_grid(true, h, 1);
        let exported = cap.render(&view, w, h);
        let mut same = view.clone();
        assert_eq!(cap.merge(&exported, &exported, &mut same, w, h).unwrap(), 0);
        assert_eq!(same, view);
        let Pixels::Gray16(mut px) = exported.clone() else { unreachable!() };
        let n = cap.side as usize;
        for (i, p) in px.iter_mut().enumerate() {
            let (dx, dy) = ((i % n) as f64 + 0.5 - n as f64 / 2.0, (i / n) as f64 + 0.5 - n as f64 / 2.0);
            if dx.hypot(dy) < n as f64 / 6.0 {
                *p = 60_000;
            }
        }
        let mut edited = view.clone();
        let changed = cap.merge(&Pixels::Gray16(px), &exported, &mut edited, w, h).unwrap();
        let Pixels::Gray16(e) = edited else { unreachable!() };
        // Poleward of 80°: the top row and a bit (the disc is a third of the cap's radius).
        assert!(changed >= w as usize, "{changed}");
        for (c, &v) in e.iter().enumerate() {
            let (lat, _) = centre(c, w, h);
            if v == 60_000 {
                assert!(lat > 79.0, "cell at {lat}°");
            } else if let Pixels::Gray16(orig) = &view {
                assert_eq!(v, orig[c]);
            }
        }
    }
}
