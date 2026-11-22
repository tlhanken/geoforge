//! The pixels of a raster edit view, read and written at any size.
//!
//! Views can be exported larger than the grid (each cell repeated, so unpainted areas come
//! back exactly) and imported at any 2:1 size: shrinking averages continuous values and
//! takes the commonest colour of palette views; enlarging takes the nearest pixel.

use crate::{Result, StoreError};
use image::{ImageBuffer, Luma, Rgb};
use std::path::Path;

/// A view's pixels, row by row from the north-west corner.
#[derive(Debug, Clone, PartialEq)]
pub enum Pixels {
    /// Continuous fields.
    Gray16(Vec<u16>),
    /// Vector fields (red east, green north, blue unused).
    Rgb16(Vec<[u16; 3]>),
    /// Palette views: categories, partitions and catalogs.
    Rgb8(Vec<[u8; 3]>),
}

/// Which kind of pixels a view holds.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Depth {
    Gray16,
    Rgb16,
    Rgb8,
}

/// Apply a generic expression to whichever pixel vector `$p` holds, keeping its kind.
macro_rules! each {
    ($p:expr, |$v:ident| $e:expr) => {
        match $p {
            Pixels::Gray16($v) => Pixels::Gray16($e),
            Pixels::Rgb16($v) => Pixels::Rgb16($e),
            Pixels::Rgb8($v) => Pixels::Rgb8($e),
        }
    };
}
pub(crate) use each;

impl Pixels {
    #[must_use]
    pub fn len(&self) -> usize {
        match self {
            Self::Gray16(v) => v.len(),
            Self::Rgb16(v) => v.len(),
            Self::Rgb8(v) => v.len(),
        }
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    #[must_use]
    pub fn depth(&self) -> Depth {
        match self {
            Self::Gray16(_) => Depth::Gray16,
            Self::Rgb16(_) => Depth::Rgb16,
            Self::Rgb8(_) => Depth::Rgb8,
        }
    }

    /// Read `path` as `depth` pixels; returns them with the image's width and height.
    ///
    /// # Errors
    /// If the file cannot be read as an image.
    pub fn load(path: &Path, depth: Depth) -> Result<(Self, u32, u32)> {
        let img = image::open(path)?;
        let (w, h) = (img.width(), img.height());
        let px = match depth {
            Depth::Gray16 => Self::Gray16(img.into_luma16().into_raw()),
            Depth::Rgb16 => Self::Rgb16(img.into_rgb16().pixels().map(|p| p.0).collect()),
            Depth::Rgb8 => Self::Rgb8(img.into_rgb8().pixels().map(|p| p.0).collect()),
        };
        Ok((px, w, h))
    }

    /// Write as a PNG `width` × `height`.
    ///
    /// # Errors
    /// On IO errors.
    pub fn save(&self, width: u32, height: u32, path: &Path) -> Result<()> {
        let size = || StoreError::Invalid(format!("{} pixels do not make a {width}x{height} image", self.len()));
        match self {
            Self::Gray16(v) => ImageBuffer::<Luma<u16>, Vec<u16>>::from_raw(width, height, v.clone()).ok_or_else(size)?.save(path)?,
            Self::Rgb16(v) => ImageBuffer::<Rgb<u16>, Vec<u16>>::from_raw(width, height, v.iter().flatten().copied().collect()).ok_or_else(size)?.save(path)?,
            Self::Rgb8(v) => ImageBuffer::<Rgb<u8>, Vec<u8>>::from_raw(width, height, v.iter().flatten().copied().collect()).ok_or_else(size)?.save(path)?,
        }
        Ok(())
    }

    /// Each pixel repeated `k` × `k` times.
    #[must_use]
    pub fn enlarge(&self, width: u32, height: u32, k: u32) -> Self {
        let (w, h, k) = (width as usize, height as usize, k as usize);
        each!(self, |v| (0..h * k).flat_map(|y| (0..w * k).map(move |x| (y / k) * w + x / k)).map(|i| v[i]).collect())
    }

    /// Resampled from `from` (width, height) to `to`: averages (continuous) or the commonest
    /// colour (palette) over each target pixel's footprint when shrinking, the nearest pixel
    /// when enlarging.
    #[must_use]
    pub fn resample(&self, from: (u32, u32), to: (u32, u32)) -> Self {
        match self {
            Self::Gray16(v) => Self::Gray16(resample(v, from, to, |s| average(s.iter().map(|&x| u32::from(x))) as u16)),
            Self::Rgb16(v) => Self::Rgb16(resample(v, from, to, |s| [0, 1, 2].map(|c| average(s.iter().map(|p| u32::from(p[c]))) as u16))),
            Self::Rgb8(v) => Self::Rgb8(resample(v, from, to, commonest)),
        }
    }
}

/// For each target pixel, `pick` over the source pixels whose centres fall in its footprint
/// (or the nearest one, along an axis that grows).
fn resample<T: Copy>(src: &[T], (w0, h0): (u32, u32), (w, h): (u32, u32), pick: impl Fn(&[T]) -> T) -> Vec<T> {
    let span = |i: u32, n: u32, n0: u32| -> std::ops::Range<usize> {
        let s = f64::from(n0) / f64::from(n);
        if s <= 1.0 {
            let j = ((f64::from(i) + 0.5) * s) as usize;
            return j..j + 1;
        }
        let (a, b) = ((f64::from(i) * s - 0.5).ceil().max(0.0) as usize, ((f64::from(i + 1) * s - 0.5).ceil() as usize).min(n0 as usize));
        a..b.max(a + 1)
    };
    let mut out = Vec::with_capacity((w * h) as usize);
    let mut buf = Vec::new();
    for y in 0..h {
        let rows = span(y, h, h0);
        for x in 0..w {
            buf.clear();
            for r in rows.clone() {
                buf.extend(span(x, w, w0).map(|c| src[r * w0 as usize + c]));
            }
            out.push(pick(&buf));
        }
    }
    out
}

fn average(values: impl Iterator<Item = u32>) -> u32 {
    let (sum, n) = values.fold((0u64, 0u64), |(s, n), v| (s + u64::from(v), n + 1));
    ((sum + n / 2) / n.max(1)) as u32
}

/// The commonest colour; ties go to the one seen first.
fn commonest(colours: &[[u8; 3]]) -> [u8; 3] {
    let mut counts: Vec<([u8; 3], usize)> = Vec::new();
    for &c in colours {
        match counts.iter_mut().find(|(k, _)| *k == c) {
            Some(e) => e.1 += 1,
            None => counts.push((c, 1)),
        }
    }
    counts.iter().fold(counts[0], |best, &e| if e.1 > best.1 { e } else { best }).0
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Enlarging by repetition and shrinking back returns the original exactly.
    #[test]
    fn enlarge_then_shrink_is_exact() {
        let g = Pixels::Gray16((0..32).map(|i| i * 997 % 65_535).collect());
        let c = Pixels::Rgb8((0..32u8).map(|i| [i, 255 - i, i / 2]).collect());
        for p in [g, c] {
            let big = p.enlarge(8, 4, 3);
            assert_eq!(big.len(), 32 * 9);
            assert_eq!(big.resample((24, 12), (8, 4)), p);
        }
    }

    /// Shrinking averages continuous values and keeps the commonest colour.
    #[test]
    fn shrinking() {
        let g = Pixels::Gray16(vec![0, 10, 20, 30, 100, 110, 120, 130]);
        assert_eq!(g.resample((4, 2), (2, 1)), Pixels::Gray16(vec![55, 75]));
        let (a, b) = ([1, 2, 3], [9, 9, 9]);
        let c = Pixels::Rgb8(vec![a, b, b, a, a, b, b, b]);
        assert_eq!(c.resample((4, 2), (2, 1)), Pixels::Rgb8(vec![a, b]));
        // Enlarging takes the nearest.
        assert_eq!(Pixels::Gray16(vec![1, 2]).resample((2, 1), (4, 2)), Pixels::Gray16(vec![1, 1, 2, 2, 1, 1, 2, 2]));
    }
}
