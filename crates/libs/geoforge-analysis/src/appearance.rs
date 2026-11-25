//! What a body looks like from space, in true colour, without clouds.
//!
//! - **A generated surface** ([`surface`]): seas darkening with depth (turquoise over
//!   shelves, navy over the abyss), sea ice where the sea freezes (white where it lasts the
//!   year), land in the colour each biome shows from orbit (forests dark green, deserts sand,
//!   tundra olive-brown), snow and ice where the year is cold enough. A dry world's bare
//!   ground is grey regolith where there is no air to redden it (darker on basalt, lighter on
//!   anorthosite and granite) and rust where there is. A sky thick enough to hide the ground
//!   (over 30 bar, as Venus's) shows its cloud deck.
//! - **A body with no generated surface** ([`body`]): drawn from what it is and its seed — a
//!   gas giant's bands coloured by its temperature (dark and red when hot, tan and cream like
//!   Jupiter's, pale gold like Saturn's when cold), an ice giant's haze cyan or blue, an ocean
//!   world's seas with scattered land, an icy body's frost, an airless rock's grey or rusty
//!   cratered plains, a cloud-wrapped world's deck.
//! - **A star** ([`star_rgb`]): its colour by temperature, as a blackbody's to the eye.
//!
//! Colours are sRGB as seen in daylight, with no shading: a viewer lights and shades them.

use crate::oceans::{Sea, SeaIce};
use geoforge_core::body::{BodyClass, BodyContext, ShellKind, Species};
use geoforge_core::grid::{AnyGrid, Grid};
use geoforge_core::product::Products;
use geoforge_core::seed::Seed;

/// Above this surface pressure the cloud deck hides the ground (Venus: 92 bar), bar.
const HIDDEN_BELOW_BAR: f64 = 30.0;
/// Snow lies all year where the mean annual temperature is below the first and covers the
/// ground completely below the second (the snow line of Earth's mountains), °C.
const SNOW_FROM_C: f64 = -6.0;
const SNOW_FULL_C: f64 = -14.0;

const SEA_SHELF: [u8; 3] = [36, 98, 128];
const SEA_DEEP: [u8; 3] = [9, 24, 62];
const LAKE: [u8; 3] = [28, 58, 86];
const ICE: [u8; 3] = [236, 240, 246];
const REGOLITH: [u8; 3] = [128, 124, 118];
const RUST: [u8; 3] = [170, 104, 66];
const CLOUD_DECK: [u8; 3] = [228, 212, 168];

/// Each biome's colour from orbit; `None` for one this table does not know (drawn as its
/// own map colour, darkened).
fn biome_rgb(slug: &str) -> Option<[u8; 3]> {
    Some(match slug {
        "ocean" => SEA_DEEP,
        "lake" => LAKE,
        "ice_sheet" => ICE,
        "tundra" => [118, 112, 92],
        "boreal_forest" => [36, 56, 38],
        "temperate_rainforest" => [32, 64, 36],
        "temperate_deciduous_forest" => [54, 78, 40],
        "temperate_grassland" => [128, 124, 80],
        "mediterranean_woodland" => [108, 100, 66],
        "cold_desert" => [166, 148, 116],
        "hot_desert" => [212, 182, 136],
        "savanna" => [148, 132, 78],
        "tropical_seasonal_forest" => [68, 88, 40],
        "tropical_rainforest" => [28, 60, 28],
        "alpine" => [126, 120, 112],
        "wetland" => [58, 76, 50],
        _ => return None,
    })
}

/// Each cell's colour on a generated surface; `None` until there is a heightmap to read.
#[must_use]
pub fn surface(products: &Products, grid: &AnyGrid, body: &BodyContext) -> Option<Vec<[u8; 3]>> {
    let heights = products.elevation().or(products.rough_elevation())?;
    let n = grid.len();
    if body.interior.surface_pressure_bar() > HIDDEN_BELOW_BAR {
        return Some(vec![CLOUD_DECK; n]);
    }
    let wet = body.interior.has_surface_liquid_water();
    let ice = Sea::read(products, grid as &dyn Grid, body.exterior.radius_km, wet).map(|s| s.ice());
    let temperature = products.temperature();
    let biomes = products.biome();
    let rocks = products.lithology();
    let bare = if body.interior.has_atmosphere() { RUST } else { REGOLITH };
    Some(
        (0..n)
            .map(|c| {
                let h = f64::from(heights.values[c]);
                let t = temperature.map(|t| f64::from(t.values[c]));
                if wet && h < 0.0 {
                    let deep = ((-h).max(1.0).log10() / 3.7).clamp(0.0, 1.0);
                    let sea = mix(SEA_SHELF, SEA_DEEP, deep);
                    return match ice.as_ref().map(|i| i[c]) {
                        Some(SeaIce::Perennial) => ICE,
                        Some(SeaIce::Seasonal) => mix(sea, ICE, 0.35),
                        _ => sea,
                    };
                }
                let ground = match biomes.and_then(|b| b.row_at(c as u32)) {
                    Some(b) => biome_rgb(&b.slug).unwrap_or_else(|| b.color.map(|v| (f64::from(v) * 0.6) as u8)),
                    None => {
                        let tone = rocks.and_then(|r| r.row_at(c as u32)).map_or(1.0, |r| match () {
                            () if r.slug.contains("basalt") => 0.75,
                            () if r.slug.contains("anorthosite") || r.slug.contains("granite") => 1.2,
                            () => 1.0,
                        });
                        bare.map(|v| (f64::from(v) * tone).min(255.0) as u8)
                    }
                };
                match t {
                    Some(t) if wet || body.interior.has_atmosphere() => mix(ground, ICE, ((SNOW_FROM_C - t) / (SNOW_FROM_C - SNOW_FULL_C)).clamp(0.0, 1.0)),
                    _ => ground,
                }
            })
            .collect(),
    )
}

/// An equirectangular image (`width` × `width`/2, rows from the north) of a body with no
/// generated surface, drawn from what it is; `seed` places its bands, lands and craters.
#[must_use]
pub fn body(b: &BodyContext, seed: Seed, width: usize) -> Vec<[u8; 3]> {
    let height = (width / 2).max(1);
    let t_k = crate::appearance::equilibrium_k(b);
    let noise = Noise(seed.0);
    let ice_share = b.interior.surface_shell().map_or(0.0, |s| if s.kind == ShellKind::IceShell { 1.0 } else { s.composition.fraction(Species::Ice) });
    let kind = match b.exterior.class {
        BodyClass::GasGiant => Look::Bands(giant_palette(t_k)),
        BodyClass::IceGiant => Look::Haze(if noise.at([0.3, 0.1, 0.7], 1.0) > 0.5 { [176, 222, 228] } else { [92, 128, 206] }),
        _ if b.interior.surface_pressure_bar() > HIDDEN_BELOW_BAR => Look::Deck,
        _ if b.interior.has_surface_liquid_water() => Look::Ocean,
        _ if ice_share > 0.5 => Look::Frost,
        _ => Look::Rock(if b.interior.has_atmosphere() { RUST } else { REGOLITH }),
    };
    let mut out = Vec::with_capacity(width * height);
    for y in 0..height {
        let lat = std::f64::consts::FRAC_PI_2 - (y as f64 + 0.5) / height as f64 * std::f64::consts::PI;
        for x in 0..width {
            let lon = (x as f64 + 0.5) / width as f64 * std::f64::consts::TAU;
            let p = [lat.cos() * lon.cos(), lat.cos() * lon.sin(), lat.sin()];
            out.push(kind.at(&noise, p, lat));
        }
    }
    out
}

/// A body's equilibrium temperature with its starlight (albedo 0.3, spread over the
/// sphere), K.
fn equilibrium_k(b: &BodyContext) -> f64 {
    278.6 * (b.irradiance.flux_w_m2 / 1361.0 * 0.7).max(0.0).powf(0.25)
}

/// A gas giant's band colours by its temperature: three tones, light to dark.
fn giant_palette(t_k: f64) -> [[u8; 3]; 3] {
    match t_k {
        t if t > 900.0 => [[120, 70, 52], [84, 44, 34], [52, 28, 24]],
        t if t > 300.0 => [[214, 196, 168], [176, 128, 92], [140, 100, 76]],
        t if t > 120.0 => [[234, 222, 196], [196, 160, 120], [150, 112, 82]],
        _ => [[232, 214, 166], [212, 190, 140], [186, 164, 118]],
    }
}

enum Look {
    Bands([[u8; 3]; 3]),
    Haze([u8; 3]),
    Deck,
    Ocean,
    Frost,
    Rock([u8; 3]),
}

impl Look {
    fn at(&self, noise: &Noise, p: [f64; 3], lat: f64) -> [u8; 3] {
        match *self {
            Self::Bands([light, mid, dark]) => {
                // Zones and belts: bands by latitude, their edges wavy.
                let wobble = 0.08 * noise.at(p, 6.0) - 0.04;
                let b = (lat + wobble) * 9.0 + 3.0 * noise.at([0.0, 0.0, lat.sin()], 4.0);
                let s = 0.5 + 0.5 * b.sin();
                let base = if s > 0.6 { mix(mid, light, (s - 0.6) / 0.4) } else { mix(dark, mid, s / 0.6) };
                shade(base, 0.9 + 0.2 * noise.at(p, 20.0))
            }
            Self::Haze(c) => shade(c, 0.96 + 0.06 * (lat * 6.0).sin().abs() + 0.03 * noise.at(p, 8.0)),
            Self::Deck => shade(CLOUD_DECK, 0.92 + 0.1 * noise.at([p[0], p[1], p[2] * 4.0], 3.0)),
            Self::Ocean => {
                let land = noise.at(p, 2.5) * 0.5 + noise.at(p, 6.0) * 0.3 + noise.at(p, 15.0) * 0.2;
                let cold = lat.abs() > 1.2;
                match () {
                    () if cold => ICE,
                    () if land > 0.62 => mix([92, 104, 62], [150, 132, 96], noise.at(p, 12.0)),
                    () => mix(SEA_SHELF, SEA_DEEP, ((0.62 - land) * 4.0).clamp(0.0, 1.0)),
                }
            }
            Self::Frost => shade(mix([226, 228, 232], [186, 164, 140], noise.at(p, 4.0).powi(2)), 0.9 + 0.15 * noise.at(p, 18.0)),
            Self::Rock(c) => {
                let crater = (noise.at(p, 16.0) - 0.5).abs() < 0.03;
                let tone = 0.75 + 0.35 * noise.at(p, 3.0) + 0.1 * noise.at(p, 24.0);
                shade(c, if crater { tone * 1.25 } else { tone })
            }
        }
    }
}

/// Smooth value noise on the unit sphere, 0–1, at `scale` cells per radius.
struct Noise(u64);

impl Noise {
    fn at(&self, p: [f64; 3], scale: f64) -> f64 {
        let q = p.map(|v| v * scale + 100.0);
        let i = q.map(|v| v.floor());
        let f = [0, 1, 2].map(|k| {
            let t = q[k] - i[k];
            t * t * (3.0 - 2.0 * t)
        });
        let mut acc = 0.0;
        for corner in 0..8u32 {
            let o = [corner & 1, (corner >> 1) & 1, (corner >> 2) & 1];
            let w: f64 = (0..3).map(|k| if o[k] == 1 { f[k] } else { 1.0 - f[k] }).product();
            let h = (0..3).fold(self.0 ^ 0x9E37_79B9_7F4A_7C15, |h, k| mix64(h ^ ((i[k] as i64 + i64::from(o[k])) as u64).wrapping_mul(0x100_0000_01B3 + k as u64)));
            acc += w * (h >> 11) as f64 / (1u64 << 53) as f64;
        }
        acc
    }
}

fn mix64(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// How a star of any kind is drawn: living stars and white dwarfs as [`star_rgb`]; brown
/// dwarfs by their spectra, not a blackbody's (sodium and potassium swallow their yellow and
/// green, and methane their near infrared): L dwarfs a deep red-orange, T dwarfs a dim
/// magenta, Y dwarfs nearly black-purple (Burgasser 2009's rendered colours); a neutron star
/// a hard blue-white; a black hole the violet of the lensed light round it.
#[must_use]
pub fn star_colour(star: &geoforge_core::body::Star) -> [u8; 3] {
    use geoforge_core::body::StellarStage;
    match star.stage {
        StellarStage::BrownDwarf if star.temperature_k >= 1300.0 => [200, 70, 30],
        StellarStage::BrownDwarf if star.temperature_k >= 550.0 => [150, 50, 120],
        StellarStage::BrownDwarf => [80, 35, 75],
        StellarStage::NeutronStar => [190, 210, 255],
        StellarStage::BlackHole => [150, 110, 255],
        _ => star_rgb(star.temperature_k),
    }
}

/// A star's colour at `temperature_k`: a blackbody's, as the eye sees it in daylight (the
/// Sun near white, red dwarfs orange, the hottest blue-white). Tanner Helland's fit to
/// Charity's blackbody table, good to a few percent from 1,000 to 40,000 K.
#[must_use]
pub fn star_rgb(temperature_k: f64) -> [u8; 3] {
    let t = temperature_k.clamp(1000.0, 40000.0) / 100.0;
    let r = if t <= 66.0 { 255.0 } else { 329.698_727_446 * (t - 60.0).powf(-0.133_204_759_2) };
    let g = if t <= 66.0 { 99.470_802_586_1 * t.ln() - 161.119_568_166_1 } else { 288.122_169_528_3 * (t - 60.0).powf(-0.075_514_849_2) };
    let b = if t >= 66.0 {
        255.0
    } else if t <= 19.0 {
        0.0
    } else {
        138.517_731_223_1 * (t - 10.0).ln() - 305.044_792_730_7
    };
    [r, g, b].map(|v| v.clamp(0.0, 255.0) as u8)
}

fn mix(a: [u8; 3], b: [u8; 3], t: f64) -> [u8; 3] {
    let t = t.clamp(0.0, 1.0);
    [0, 1, 2].map(|k| (f64::from(a[k]) * (1.0 - t) + f64::from(b[k]) * t).round() as u8)
}

fn shade(c: [u8; 3], f: f64) -> [u8; 3] {
    c.map(|v| (f64::from(v) * f).clamp(0.0, 255.0) as u8)
}

#[cfg(test)]
mod tests {
    use super::*;
    use geoforge_core::body::Archetype;
    use geoforge_core::grid::GridSpec;
    use geoforge_core::key::Key;
    use geoforge_core::product::ProductData;
    use geoforge_core::shapes::Field;

    fn luminance(c: [u8; 3]) -> f64 {
        0.2126 * f64::from(c[0]) + 0.7152 * f64::from(c[1]) + 0.0722 * f64::from(c[2])
    }

    #[test]
    fn stars_run_from_orange_to_blue() {
        let (m, sun, a) = (star_rgb(3200.0), star_rgb(5772.0), star_rgb(10000.0));
        assert!(m[0] == 255 && m[1] < 200 && m[2] < 140, "red dwarf {m:?} (Charity: 255, 180, 107)");
        assert!(sun.iter().all(|&v| v > 225) && sun[0] >= sun[2], "the Sun {sun:?}");
        assert!(a[2] == 255 && a[0] < 230, "an A star {a:?}");
    }

    /// Sea darkens with depth, the ground snows over when cold, land takes its biome's
    /// colour, and an unknown planet's face is the same for the same seed.
    #[test]
    fn a_surface_in_true_colour() {
        let spec = GridSpec::Equirect { width: 4, height: 2 };
        let grid = spec.build().unwrap();
        let earth = Archetype::Earth.context(Key(1), "Earth");
        let mut p = Products::new();
        p.insert(ProductData::Elevation(Field::from_values(spec, vec![-50.0, -5000.0, 300.0, 3000.0, 100.0, 100.0, 100.0, 100.0])));
        p.insert(ProductData::Temperature(Field::from_values(spec, vec![15.0, 15.0, 15.0, -25.0, 15.0, 15.0, 15.0, 15.0])));
        let c = surface(&p, &grid, &earth).unwrap();
        assert!(c[0][2] > c[0][0] && luminance(c[0]) > luminance(c[1]), "shelf {:?} brighter than the abyss {:?}", c[0], c[1]);
        assert!(luminance(c[3]) > 220.0, "snow {:?}", c[3]);
        assert!(luminance(c[2]) < 200.0, "bare ground {:?}", c[2]);
        let moon = Archetype::Moon.context(Key(2), "Moon");
        let face = body(&moon, Seed(3), 64);
        assert_eq!(face.len(), 64 * 32);
        assert_eq!(face, body(&moon, Seed(3), 64));
        assert!(face.iter().all(|c| (c[0] as i32 - c[2] as i32).abs() < 30), "grey");
    }

    /// A gas giant is banded: colour changes far more down the latitudes than along them.
    #[test]
    fn giants_are_banded() {
        let jupiter = Archetype::GasGiant.context(Key(3), "Jupiter");
        let (w, h) = (128, 64);
        let face = body(&jupiter, Seed(5), w);
        let lum: Vec<f64> = face.iter().map(|&c| luminance(c)).collect();
        let along: f64 = (0..h).map(|y| (1..w).map(|x| (lum[y * w + x] - lum[y * w + x - 1]).abs()).sum::<f64>()).sum::<f64>() / (h * (w - 1)) as f64;
        let down: f64 = (0..w).map(|x| (1..h).map(|y| (lum[y * w + x] - lum[(y - 1) * w + x]).abs()).sum::<f64>()).sum::<f64>() / (w * (h - 1)) as f64;
        assert!(down > 1.5 * along, "down {down}, along {along}");
    }
}
