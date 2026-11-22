//! The stars around a system and the glow of its galaxy, as its night sky shows them.
//!
//! The neighbourhood is drawn from the galaxy's own population, as the system was. Stars near
//! the system's orbit around the galactic centre number 0.09 per cubic parsec at the
//! reference radius (the Sun's neighbourhood: 0.08–0.1 stars per pc³ counting white dwarfs,
//! RECONS 10 pc census), thinning outward with the disk's scale length. Old stars sit in a
//! thicker layer than young ones: the scale height is 50 pc at birth growing to 300 pc
//! (OB stars 45–60 pc, the old thin disk ~300 pc; Bovy et al. 2012). Dust dims light by
//! 1 magnitude per kiloparsec (0.7–1.8 in the plane, Allen's Astrophysical Quantities).
//!
//! Only stars the eye can see (to magnitude 6.5) are drawn. Bright stars are rare and seen
//! from far, faint ones common and seen only nearby, so a census by draw-and-reject would
//! need hundreds of millions of draws. Instead kinds of star are sampled with masses
//! log-uniform and ages within each kind's life (importance sampling), each kind standing
//! for a density of stars; the volume within which a kind is visible gives the expected
//! number seen, and that many stars are placed. The Sun's sky has 9,110 such stars (Yale
//! Bright Star Catalogue) and 22 brighter than magnitude 1.5.
//!
//! Each star may have an unresolved companion of the same age, as often as the galaxy's
//! systems do, adding its light (up to 0.75 magnitude for a twin).
//!
//! The galaxy's diffuse glow is the starlight of its disk along each line of sight, dimmed
//! by dust in a thinner layer (scale height 100 pc) that also thickens toward the centre.

use crate::draw::{name, Draw};
use crate::galaxy::Galaxy;
use crate::population::multiple_fraction;
use crate::stars::{end_of_life_gyr, evolve};
use geoforge_core::body::{Star, StellarStage};
use geoforge_core::seed::Seed;

/// Faintest magnitude the eye sees on a dark night.
pub const NAKED_EYE: f64 = 6.5;
/// Stars per cubic parsec at the galaxy's reference radius.
const LOCAL_DENSITY_PC3: f64 = 0.09;
/// Dust dimming in the plane, magnitudes per kiloparsec.
const EXTINCTION_MAG_KPC: f64 = 1.0;
/// Dust layer scale height, pc.
const DUST_HEIGHT_PC: f64 = 100.0;
/// Kinds of star sampled.
const KINDS: usize = 200_000;
/// Stars named, brightest first.
const NAMED: usize = 60;
/// Light-years per parsec.
pub const LY_PER_PC: f64 = 3.261_56;

/// A star of the neighbourhood.
#[derive(Debug, Clone, PartialEq, serde::Serialize)]
pub struct Neighbour {
    /// The brightest have names.
    pub name: Option<String>,
    // `star` is the brighter of a pair; `absolute_v` and `v` count both.
    /// Unit vector toward the star in the system's frame (that of its orbital elements).
    pub direction: [f64; 3],
    pub distance_pc: f64,
    pub star: Star,
    /// Absolute and apparent visual magnitudes.
    pub absolute_v: f64,
    pub v: f64,
}

/// The stars a system sees, and its place in the galaxy.
#[derive(Debug, Clone, PartialEq)]
pub struct Neighbourhood {
    /// Brightest first.
    pub stars: Vec<Neighbour>,
    /// Rows: the system frame's axes in galactic coordinates (x toward the centre, z toward
    /// the north galactic pole). A galactic vector g is `to_system · g` in the system frame.
    pub to_system: [[f64; 3]; 3],
    pub radius_kpc: f64,
    galaxy: Galaxy,
}

/// The naked-eye stars around a system `radius_kpc` from the centre of `galaxy`.
#[must_use]
pub fn neighbourhood(galaxy: &Galaxy, radius_kpc: f64, seed: Seed) -> Neighbourhood {
    let mut d = Draw::new(seed.child(b"neighbourhood"));
    let to_system = random_rotation(&mut d);
    let density = LOCAL_DENSITY_PC3 * (-(radius_kpc - galaxy.layout.reference_kpc) / galaxy.layout.scale_length_kpc).exp();
    let (m_lo, m_hi) = (galaxy.imf.breaks[0], galaxy.imf.breaks[galaxy.imf.breaks.len() - 1]);
    let span = (m_hi / m_lo).ln();
    let history = &galaxy.history;
    let metals = &galaxy.metals;

    // Each kind: (star, absolute V, scale height, reach in pc, expected number seen).
    let mut kinds = Vec::with_capacity(KINDS);
    let mut total = 0.0;
    for _ in 0..KINDS {
        let m = m_lo * (d.u() * span).exp();
        // Ages within the star's shining life (remnants are too faint to see).
        let life = end_of_life_gyr(m).min(history.age_gyr);
        let age = d.u() * life;
        let weight = galaxy.imf.density(m) * m * span * history.lookback_density(age) * life;
        let feh = (metals.at_reference_dex + metals.gradient_dex_per_kpc * (radius_kpc - galaxy.layout.reference_kpc) + metals.per_gyr_dex * age + d.normal(0.0, metals.scatter_dex)).clamp(-2.5, 0.6);
        let star = evolve(m, age, feh);
        let Some(mut abs_v) = absolute_v(&star) else {
            continue;
        };
        // A companion of the same age, unresolved: mass ratios flat (Duchêne & Kraus 2013).
        if d.u() < multiple_fraction(m) {
            let mate = evolve((m * d.u()).max(m_lo), age, feh);
            if let Some(v) = absolute_v(&mate) {
                abs_v = -2.5 * (10f64.powf(-0.4 * abs_v) + 10f64.powf(-0.4 * v)).log10();
            }
        }
        let reach = reach_pc(abs_v);
        let height = scale_height_pc(age);
        let expected = density * weight / KINDS as f64 * disk_volume_pc3(reach, height);
        if expected > 0.0 {
            total += expected;
            kinds.push((star, abs_v, height, reach, total));
        }
    }

    let count = poisson(&mut d, total);
    let mut stars: Vec<Neighbour> = (0..count)
        .map(|_| {
            let u = d.u() * total;
            let i = kinds.partition_point(|k| k.4 < u).min(kinds.len() - 1);
            let (star, abs_v, height, reach, _) = &kinds[i];
            let (g, dist) = place(&mut d, *reach, *height);
            let direction = apply(&to_system, g);
            Neighbour { name: None, direction, distance_pc: dist, star: star.clone(), absolute_v: *abs_v, v: apparent(*abs_v, dist) }
        })
        .collect();
    stars.sort_by(|a, b| a.v.total_cmp(&b.v));
    for (i, s) in stars.iter_mut().take(NAMED).enumerate() {
        s.name = Some(name(seed.child_u64(b"star name", i as u64)));
    }
    Neighbourhood { stars, to_system, radius_kpc, galaxy: galaxy.clone() }
}

impl Neighbourhood {
    /// Relative surface brightness of the galaxy's glow toward `direction` (system frame):
    /// the disk's starlight along the line of sight, dimmed by dust, in units of the
    /// brightness toward the galactic pole.
    #[must_use]
    pub fn glow(&self, direction: [f64; 3]) -> f64 {
        let t = &self.to_system;
        // The transpose takes the system frame back to galactic coordinates.
        let g = [
            t[0][0] * direction[0] + t[1][0] * direction[1] + t[2][0] * direction[2],
            t[0][1] * direction[0] + t[1][1] * direction[1] + t[2][1] * direction[2],
            t[0][2] * direction[0] + t[1][2] * direction[1] + t[2][2] * direction[2],
        ];
        self.line_of_sight(g) / self.line_of_sight([0.0, 0.0, 1.0])
    }

    fn line_of_sight(&self, g: [f64; 3]) -> f64 {
        const STEP_PC: f64 = 50.0;
        const HEIGHT_PC: f64 = 300.0;
        let (r0, h_r) = (self.radius_kpc * 1000.0, self.galaxy.layout.scale_length_kpc * 1000.0);
        let (mut light, mut dimming) = (0.0, 0.0);
        let mut s = 0.5 * STEP_PC;
        while s < 30_000.0 {
            // The system sits at (r0, 0, 0) from the centre; x points toward the centre.
            let (x, y, z) = (r0 - s * g[0], s * g[1], s * g[2]);
            let radial = (-((x * x + y * y).sqrt() - r0) / h_r).exp();
            let stars = radial * (-z.abs() / HEIGHT_PC).exp();
            dimming += EXTINCTION_MAG_KPC * STEP_PC / 1000.0 * radial * (-z.abs() / DUST_HEIGHT_PC).exp();
            light += stars * 10f64.powf(-0.4 * dimming) * STEP_PC;
            if stars < 1e-4 && dimming > 0.0 {
                break;
            }
            s += STEP_PC;
        }
        light
    }
}

/// Absolute visual magnitude, if the star shines.
#[must_use]
pub fn absolute_v(star: &Star) -> Option<f64> {
    if star.luminosity_solar <= 0.0 || matches!(star.stage, StellarStage::BlackHole | StellarStage::NeutronStar) {
        return None;
    }
    Some(SOLAR_BOLOMETRIC - 2.5 * star.luminosity_solar.log10() - bolometric_correction(star.temperature_k))
}

/// The Sun's absolute bolometric magnitude (IAU 2015 B2).
const SOLAR_BOLOMETRIC: f64 = 4.74;

/// Visual bolometric correction at effective temperature `t` K: Torres (2010)'s fit to
/// Flower (1996), −0.07 for the Sun.
#[must_use]
pub fn bolometric_correction(t: f64) -> f64 {
    let x = t.clamp(2000.0, 60_000.0).log10();
    let c: &[f64] = if x < 3.70 {
        &[-0.190_537_291_496_456e5, 0.155_144_866_764_412e5, -0.421_278_819_301_717e4, 0.381_476_328_422_343e3]
    } else if x < 3.90 {
        &[-0.370_510_203_809_015e5, 0.385_672_629_965_804e5, -0.150_651_486_316_025e5, 0.261_724_637_119_416e4, -0.170_623_810_323_864e3]
    } else {
        &[-0.118_115_450_538_963e6, 0.137_145_973_583_929e6, -0.636_233_812_100_225e5, 0.147_412_923_562_646e5, -0.170_587_278_406_872e4, 0.788_731_721_804_990e2]
    };
    c.iter().rev().fold(0.0, |acc, k| acc * x + k)
}

/// Apparent magnitude at `distance_pc`, dimmed by dust.
#[must_use]
pub fn apparent(absolute_v: f64, distance_pc: f64) -> f64 {
    absolute_v + 5.0 * (distance_pc / 10.0).log10() + EXTINCTION_MAG_KPC * distance_pc / 1000.0
}

/// How far a star of `absolute_v` is seen with the naked eye, pc.
fn reach_pc(absolute_v: f64) -> f64 {
    let (mut lo, mut hi) = (1e-3f64.ln(), 1e5f64.ln());
    for _ in 0..50 {
        let mid = 0.5 * (lo + hi);
        if apparent(absolute_v, mid.exp()) < NAKED_EYE {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    lo.exp()
}

/// Scale height of stars born `age_gyr` ago, pc.
fn scale_height_pc(age_gyr: f64) -> f64 {
    50.0 + 250.0 * (1.0 - (-age_gyr / 2.0).exp())
}

/// ∫ e^(−|z|/h) dV over a sphere of radius `r` centred in the plane, pc³.
fn disk_volume_pc3(r: f64, h: f64) -> f64 {
    let a = r / h;
    if a < 1e-3 {
        return 4.0 / 3.0 * std::f64::consts::PI * r.powi(3);
    }
    let e = (-a).exp();
    2.0 * std::f64::consts::PI * (r * r * h * (1.0 - e) - h.powi(3) * (2.0 - e * (a * a + 2.0 * a + 2.0)))
}

/// A point within `reach` pc, weighted by e^(−|z|/h): its galactic direction and distance.
fn place(d: &mut Draw, reach: f64, h: f64) -> ([f64; 3], f64) {
    loop {
        let r = reach * d.u().cbrt();
        let z = 2.0 * d.u() - 1.0;
        let phi = std::f64::consts::TAU * d.u();
        let s = (1.0 - z * z).sqrt();
        let dir = [s * phi.cos(), s * phi.sin(), z];
        if d.u() < (-(r * z).abs() / h).exp() {
            return (dir, r.max(1.0));
        }
    }
}

/// A Poisson count with mean `mean`.
fn poisson(d: &mut Draw, mean: f64) -> usize {
    if mean > 50.0 {
        return d.normal(mean, mean.sqrt()).round().max(0.0) as usize;
    }
    let (limit, mut p, mut k) = ((-mean).exp(), 1.0, 0);
    loop {
        p *= d.u();
        if p <= limit {
            return k;
        }
        k += 1;
    }
}

/// A rotation drawn evenly over all orientations (Shoemake 1992), as a matrix.
fn random_rotation(d: &mut Draw) -> [[f64; 3]; 3] {
    let (u1, u2, u3) = (d.u(), d.u(), d.u());
    let tau = std::f64::consts::TAU;
    let (a, b) = ((1.0 - u1).sqrt(), u1.sqrt());
    let (x, y, z, w) = (a * (tau * u2).sin(), a * (tau * u2).cos(), b * (tau * u3).sin(), b * (tau * u3).cos());
    [
        [1.0 - 2.0 * (y * y + z * z), 2.0 * (x * y - z * w), 2.0 * (x * z + y * w)],
        [2.0 * (x * y + z * w), 1.0 - 2.0 * (x * x + z * z), 2.0 * (y * z - x * w)],
        [2.0 * (x * z - y * w), 2.0 * (y * z + x * w), 1.0 - 2.0 * (x * x + y * y)],
    ]
}

fn apply(m: &[[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [0, 1, 2].map(|i| m[i][0] * v[0] + m[i][1] * v[1] + m[i][2] * v[2])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bolometric_corrections() {
        for (t, bc) in [(5772.0, -0.07), (3500.0, -2.0), (10_000.0, -0.4), (30_000.0, -2.9)] {
            assert!((bolometric_correction(t) - bc).abs() < 0.35, "{t} K: {}", bolometric_correction(t));
        }
    }

    /// A sky like the Sun's: about 9,000 naked-eye stars, about 20 of the first magnitude,
    /// and the glow brightest along a band.
    #[test]
    fn a_sky_like_ours() {
        let mw = Galaxy::milky_way();
        let n = neighbourhood(&mw, mw.layout.reference_kpc, Seed(7));
        let count = n.stars.len();
        let first = n.stars.iter().filter(|s| s.v < 1.5).count();
        let giants = n.stars.iter().filter(|s| s.star.stage == StellarStage::Giant).count();
        assert!((6_000..13_000).contains(&count), "{count} naked-eye stars");
        // Giants are the commonest naked-eye stars, but not all of them.
        assert!((0.4..0.75).contains(&(giants as f64 / count as f64)), "{giants} giants");
        assert!((8..45).contains(&first), "{first} of the first magnitude");
        // The plane outshines the pole several times over.
        let t = &n.to_system;
        let plane = apply(t, [0.0, 1.0, 0.0]);
        assert!(n.glow(plane) > 3.0, "{}", n.glow(plane));
        let toward_centre = apply(t, [1.0, 0.0, 0.0]);
        assert!(n.glow(toward_centre) > n.glow(apply(t, [-1.0, 0.0, 0.0])));
    }

    #[test]
    fn disk_volume_limits() {
        let small = disk_volume_pc3(1.0, 300.0);
        assert!((small / (4.0 / 3.0 * std::f64::consts::PI) - 1.0).abs() < 0.01);
        // A wide sphere holds a slab: 2πr²h.
        let wide = disk_volume_pc3(1e5, 100.0);
        assert!((wide / (2.0 * std::f64::consts::PI * 1e10 * 100.0) - 1.0).abs() < 0.01);
    }
}
