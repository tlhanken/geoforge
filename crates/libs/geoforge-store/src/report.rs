//! A report: one HTML page and its images, written to a folder (`index.html`, `img/*.jpg`)
//! that opens in any browser or publishes as is.

use crate::{Result, StoreError};
use image::codecs::jpeg::JpegEncoder;
use image::RgbImage;
use std::path::{Path, PathBuf};

/// JPEG quality of the report's images.
const QUALITY: u8 = 85;

/// Write `html` as `dir/index.html` and each image as `dir/img/<name>.jpg`, where the page
/// refers to it. Returns the files written.
///
/// # Errors
/// On IO or encoding errors.
pub fn write(dir: &Path, html: &str, images: &[(String, RgbImage)]) -> Result<Vec<PathBuf>> {
    let img_dir = dir.join("img");
    std::fs::create_dir_all(&img_dir)?;
    let mut out = Vec::with_capacity(images.len() + 1);
    for (name, img) in images {
        let path = img_dir.join(format!("{name}.jpg"));
        let file = std::io::BufWriter::new(std::fs::File::create(&path)?);
        JpegEncoder::new_with_quality(file, QUALITY).encode_image(img).map_err(StoreError::from)?;
        out.push(path);
    }
    let page = dir.join("index.html");
    std::fs::write(&page, html)?;
    out.push(page);
    Ok(out)
}
