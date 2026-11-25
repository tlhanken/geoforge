//! The stars around a system and the glow of its galaxy, as its night sky shows them.
//!
//! The stars are the galaxy's own ([`crate::starfield`]): the generated stars and a world's
//! edits around the system's place, each dimmed by the dust along its line of sight
//! ([`crate::sectors::dust`]), those brighter than the limit kept. The system's own star is
//! left out (a system made at a star), or any star within half a parsec of it (a system
//! placed on its own). The Sun's sky has 9,110 stars to magnitude 6.5 (Yale Bright Star
//! Catalogue) and 22 brighter than 1.5.
//!
//! The glow is the light of the stars too faint to see one by one: every population's
//! density along each line of sight times its light from stars below the limit at that
//! distance, dimmed by the dust, on a 3° grid of galactic longitude and latitude.
//!
//! Where stars crowd far more than near the Sun (a bulge, a galaxy's core), the stars listed
//! one by one are cut to a brighter limit, keeping about as many as four Suns' skies would
//! hold; the rest join the glow.

use crate::draw::Draw;
use crate::sectors::bins::{bright_edge, BINS, LIT_BINS};
use crate::sectors::dust::DustView;
use crate::starfield::Starfield;
use geoforge_core::body::{SolarSystem, Star, StellarStage};
use geoforge_core::galactic::StarId;
use rayon::prelude::*;

/// Faintest magnitude the eye sees on a dark night.
pub const NAKED_EYE: f64 = 6.5;
/// Light-years per parsec.
pub const LY_PER_PC: f64 = 3.261_56;
/// Stars named, brightest first.
const NAMED: usize = 60;
/// Farthest a star is looked for, pc.
const MAX_PC: f64 = 10_000.0;
/// Generated stars this close to a system placed on its own are taken to be it, pc.
const CLEAR_PC: f64 = 0.5;
/// Systems per pc³ near the Sun (the 10 pc census), and how many times that a sky lists in
/// full before cutting its limit.
const SUN_DENSITY: f64 = 0.08;
const CROWDING: f64 = 4.0;
/// The glow grid, degrees.
const GLOW_DEG: f64 = 3.0;
/// The glow's lines of sight: from 1 pc out by this factor a step, to 30 kpc.
const GLOW_STEP: f64 = 1.06;
const GLOW_REACH_PC: f64 = 30_000.0;

/// Where a sky is seen from and how its frame sits.
#[derive(Debug, Clone, PartialEq)]
pub struct Observer {
    /// pc from the galaxy's centre.
    pub position_pc: [f64; 3],
    /// Rows: the system frame's axes in galactic coordinates as seen from here (x toward the
    /// centre, z toward the north pole). A galactic vector g is `to_system · g` in the
    /// system's frame.
    pub to_system: [[f64; 3]; 3],
    /// The system's own star, if it was made at one.
    pub exclude: Option<StarId>,
    /// Generated stars this close are taken to be the system itself, pc.
    pub clear_pc: f64,
}

impl Observer {
    /// The view from `system`'s place: its site's position and orientation, or for a system
    /// without them, one at its site's radius (the reference radius if it has no site) in the
    /// plane, both drawn from its seed.
    #[must_use]
    pub fn of_system(system: &SolarSystem, reference_kpc: f64) -> Self {
        let site = system.root.site.as_ref();
        let mut d = Draw::new(system.root.seed.child(b"place"));
        let azimuth = std::f64::consts::TAU * d.u();
        let r = site.map_or(reference_kpc, |s| s.radius_kpc) * 1000.0;
        let position_pc = site.and_then(|s| s.position_pc).unwrap_or([r * azimuth.cos(), r * azimuth.sin(), 0.0]);
        let to_system = site.and_then(|s| s.orientation).unwrap_or_else(|| d.rotation());
        let exclude = site.and_then(|s| s.star);
        Self { position_pc, to_system, exclude, clear_pc: if exclude.is_some() { 0.0 } else { CLEAR_PC } }
    }

    /// The galactic axes seen from here, in the galaxy's frame: toward the centre, along the
    /// rotation (longitude 90°), toward the north pole.
    #[must_use]
    pub fn axes(&self) -> [[f64; 3]; 3] {
        let [x, y, _] = self.position_pc;
        let r = x.hypot(y);
        let toward = if r > 1e-6 { [-x / r, -y / r, 0.0] } else { [1.0, 0.0, 0.0] };
        let up = [0.0, 0.0, 1.0];
        [toward, cross(up, toward), up]
    }

    /// A vector of the galaxy's frame in the system's frame.
    #[must_use]
    pub fn to_system_frame(&self, v: [f64; 3]) -> [f64; 3] {
        let a = self.axes();
        let g = a.map(|axis| dot(axis, v));
        apply(&self.to_system, g)
    }

    /// A vector of the system's frame in the galaxy's frame.
    #[must_use]
    pub fn to_galaxy_frame(&self, v: [f64; 3]) -> [f64; 3] {
        let t = &self.to_system;
        let g = [0, 1, 2].map(|k| t[0][k] * v[0] + t[1][k] * v[1] + t[2][k] * v[2]);
        let a = self.axes();
        [0, 1, 2].map(|k| a[0][k] * g[0] + a[1][k] * g[1] + a[2][k] * g[2])
    }
}

/// A star of the neighbourhood.
#[derive(Debug, Clone, PartialEq, serde::Serialize)]
pub struct Neighbour {
    pub id: StarId,
    /// The brightest have names.
    pub name: Option<String>,
    /// Unit vector toward the star in the system's frame.
    pub direction: [f64; 3],
    pub distance_pc: f64,
    /// Its place in the galaxy, pc.
    pub position_pc: [f64; 3],
    /// The brighter star of the system; `absolute_v` and `v` count all of it.
    pub star: Star,
    pub absolute_v: f64,
    pub v: f64,
    /// Dust's dimming of it, mag.
    pub extinction_mag: f64,
}

/// The stars a system sees, and its place in the galaxy.
#[derive(Debug, Clone)]
pub struct Neighbourhood {
    /// Brightest first.
    pub stars: Vec<Neighbour>,
    pub observer: Observer,
    pub galaxy: String,
    pub radius_kpc: f64,
    /// The limit its stars were listed to: the one asked for, or brighter where stars crowd.
    pub limit: f64,
    /// Glow on the grid of galactic longitude (columns) and latitude (rows), relative to the
    /// north galactic pole's.
    glow: Vec<f64>,
    /// The north galactic pole's glow, L☉ (V) per pc².
    pole_sigma: f64,
    /// Extinction out of the galaxy on the same grid, magnitudes.
    out: Vec<f64>,
}

/// The stars brighter than `limit` seen by `observer` among `field`'s, and the glow.
#[must_use]
pub fn neighbourhood(field: &Starfield, observer: &Observer, limit: f64) -> Neighbourhood {
    let here = observer.position_pc;
    // A sky's count grows as density × 10^(0.6 m): a crowd costs magnitudes off the limit.
    let crowd = field.sectors().density(here, 0..LIT_BINS) / (CROWDING * SUN_DENSITY);
    let limit = if crowd > 1.0 { limit - crowd.log10() / 0.6 } else { limit };
    let candidates = field.visible(here, limit, MAX_PC);
    let reach = candidates.iter().map(|s| s.distance_pc(here)).fold(0.0, f64::max).min(MAX_PC);
    let dust = DustView::new(field.sectors(), here, reach);
    let mut stars: Vec<Neighbour> = candidates
        .into_par_iter()
        .filter(|s| Some(s.id) != observer.exclude && s.distance_pc(here) > observer.clear_pc)
        .filter_map(|s| {
            let d = s.distance_pc(here);
            let extinction_mag = dust.to(s.position_pc);
            let v = s.absolute_v + 5.0 * (d / 10.0).log10() + extinction_mag;
            let to = [0, 1, 2].map(|k| (s.position_pc[k] - here[k]) / d);
            (v <= limit).then(|| Neighbour {
                id: s.id,
                name: None,
                direction: observer.to_system_frame(to),
                distance_pc: d,
                position_pc: s.position_pc,
                star: s.star,
                absolute_v: s.absolute_v,
                v,
                extinction_mag,
            })
        })
        .collect();
    stars.sort_by(|a, b| a.v.total_cmp(&b.v).then(a.id.cmp(&b.id)));
    for s in stars.iter_mut().take(NAMED) {
        s.name = Some(field.name_of(&s.id));
    }
    let (glow, pole_sigma, out) = glow(field, observer, &dust, limit);
    let r = here[0].hypot(here[1]) / 1000.0;
    Neighbourhood { stars, observer: observer.clone(), galaxy: field.galaxy().name.clone(), radius_kpc: r, limit, glow, pole_sigma, out }
}

const GLOW_LON: usize = (360.0 / GLOW_DEG) as usize;
const GLOW_LAT: usize = (180.0 / GLOW_DEG) as usize + 1;

/// The glow toward each grid direction, relative to the north galactic pole's; the pole's in
/// L☉ (V) per pc² (the light of the stars along the line of sight, projected: the unit of
/// surface brightness that does not depend on distance, 26.4 mag/arcsec² for 1); and the
/// extinction out of the galaxy toward each.
fn glow(field: &Starfield, observer: &Observer, dust: &DustView, limit: f64) -> (Vec<f64>, f64, Vec<f64>) {
    let s = field.sectors();
    let g = s.galaxy();
    let here = observer.position_pc;
    let mut distances = Vec::new();
    let mut d = 1.0;
    while d < GLOW_REACH_PC {
        distances.push(d);
        d *= GLOW_STEP;
    }
    // Each population's light from stars too faint to see at each distance.
    let faint = |light: &[f64; BINS], at: f64| -> f64 { (0..BINS).filter(|&b| bright_edge(b) + 0.5 + 5.0 * (at / 10.0).log10() > limit).map(|b| light[b]).sum() };
    let field_light: Vec<Vec<f64>> = s.field.iter().map(|p| distances.iter().map(|&at| faint(&p.table.light, at)).collect()).collect();
    let young_light: Vec<Vec<f64>> = s
        .young
        .iter()
        .map(|y| {
            let mut mean = [0.0; BINS];
            for t in &y.tables {
                (0..BINS).for_each(|b| mean[b] += t.light[b] / y.tables.len() as f64);
            }
            distances.iter().map(|&at| faint(&mean, at)).collect()
        })
        .collect();
    let axes = observer.axes();
    let ray = |lon: f64, lat: f64| -> (f64, f64) {
        let (l, b) = (lon.to_radians(), lat.to_radians());
        let local = [b.cos() * l.cos(), b.cos() * l.sin(), b.sin()];
        let u = [0, 1, 2].map(|k| axes[0][k] * local[0] + axes[1][k] * local[1] + axes[2][k] * local[2]);
        let dim = dust.along(u, &distances);
        let mut total = 0.0;
        let mut prev = 0.0;
        for (i, &at) in distances.iter().enumerate() {
            let step = at - prev;
            prev = at;
            let p = [0, 1, 2].map(|k| here[k] + at * u[k]);
            let excess = g.arms.as_ref().map_or(0.0, |a| a.excess(p[0], p[1]));
            let arm = |c: usize, age: f64| (1.0 + g.arm_contrast(c, age) * excess).max(0.0);
            let mut e: f64 = s.field.iter().zip(&field_light).map(|(pop, light)| pop.mass_msun * pop.shape.density(p) * arm(pop.component, pop.age) * light[i]).sum();
            e += s.young.iter().zip(&young_light).map(|(y, light)| y.mass_msun * y.shape.density(p) * arm(y.component, y.age) * light[i]).sum::<f64>();
            total += e * step * 10f64.powf(-0.4 * dim[i]);
        }
        (total, dim.last().copied().unwrap_or(0.0))
    };
    let (grid, out): (Vec<f64>, Vec<f64>) = (0..GLOW_LAT * GLOW_LON).into_par_iter().map(|i| ray((i % GLOW_LON) as f64 * GLOW_DEG, -90.0 + (i / GLOW_LON) as f64 * GLOW_DEG)).unzip();
    let pole = grid[(GLOW_LAT - 1) * GLOW_LON].max(1e-300);
    (grid.into_iter().map(|v| v / pole).collect(), pole, out)
}

impl Neighbourhood {
    /// The glow toward `direction` (the system's frame), relative to the north galactic
    /// pole's.
    #[must_use]
    pub fn glow(&self, direction: [f64; 3]) -> f64 {
        self.on_grid(&self.glow, direction)
    }

    /// `grid` (galactic longitude by latitude) toward `direction` (the system's frame).
    fn on_grid(&self, grid: &[f64], direction: [f64; 3]) -> f64 {
        let t = &self.observer.to_system;
        let g = [0, 1, 2].map(|k| t[0][k] * direction[0] + t[1][k] * direction[1] + t[2][k] * direction[2]);
        let lat = g[2].clamp(-1.0, 1.0).asin().to_degrees();
        let lon = g[1].atan2(g[0]).to_degrees().rem_euclid(360.0);
        let (fx, fy) = (lon / GLOW_DEG, (lat + 90.0) / GLOW_DEG);
        let (x0, y0) = (fx.floor() as usize % GLOW_LON, (fy.floor() as usize).min(GLOW_LAT - 2));
        let (tx, ty) = (fx - fx.floor(), (fy - y0 as f64).clamp(0.0, 1.0));
        let at = |x: usize, y: usize| grid[y * GLOW_LON + x % GLOW_LON];
        (1.0 - ty) * ((1.0 - tx) * at(x0, y0) + tx * at(x0 + 1, y0)) + ty * ((1.0 - tx) * at(x0, y0 + 1) + tx * at(x0 + 1, y0 + 1))
    }

    /// The north galactic pole's glow, L☉ (V) per pc².
    #[must_use]
    pub fn pole_sigma(&self) -> f64 {
        self.pole_sigma
    }

    /// Extinction (magnitudes) out of the galaxy toward `direction` (the system's frame): what
    /// dims another galaxy seen that way.
    #[must_use]
    pub fn extinction_out(&self, direction: [f64; 3]) -> f64 {
        self.on_grid(&self.out, direction)
    }

    /// The place in the galaxy seen `distance_pc` away in `direction` (the system's frame).
    #[must_use]
    pub fn galactic(&self, direction: [f64; 3], distance_pc: f64) -> [f64; 3] {
        let u = self.observer.to_galaxy_frame(direction);
        [0, 1, 2].map(|k| self.observer.position_pc[k] + distance_pc * u[k])
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
    // Brown dwarfs, below 2000 K, from Pecaut & Mamajek (2013)'s table: −6 at the L
    // dwarfs' start, falling 5 mag a halving of temperature (a 1000 K T dwarf at M_V ≈ 28).
    // (The fit itself is good only above ~3200 K: late M dwarfs come out up to a few
    // magnitudes too faint.)
    if t < 2000.0 {
        return -6.0 - 5.0 * (2000.0 / t.max(100.0)).log2();
    }
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

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

fn apply(m: &[[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [0, 1, 2].map(|i| dot(m[i], v))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::galaxy::Galaxy;
    use geoforge_core::galactic::StarEdits;
    use geoforge_core::seed::Seed;

    #[test]
    fn bolometric_corrections() {
        for (t, bc) in [(5772.0, -0.07), (3500.0, -2.0), (10_000.0, -0.4), (30_000.0, -2.9)] {
            assert!((bolometric_correction(t) - bc).abs() < 0.35, "{t} K: {}", bolometric_correction(t));
        }
    }

    fn sun(seed: u64) -> Observer {
        let mut d = Draw::new(Seed(seed));
        Observer { position_pc: [8200.0, 0.0, 25.0], to_system: d.rotation(), exclude: None, clear_pc: CLEAR_PC }
    }

    /// A sky like the Sun's: about 9,000 naked-eye stars, about 20 of the first magnitude,
    /// giants the commonest, and the glow brightest along a band and toward the centre.
    #[test]
    fn a_sky_like_ours() {
        let field = Starfield::new(&Galaxy::milky_way(), Seed(7), StarEdits::default());
        let o = sun(7);
        let n = neighbourhood(&field, &o, NAKED_EYE);
        let count = n.stars.len();
        let first = n.stars.iter().filter(|s| s.v < 1.5).count();
        let giants = n.stars.iter().filter(|s| s.star.stage == StellarStage::Giant).count();
        assert!((6_000..13_000).contains(&count), "{count} naked-eye stars");
        assert!((0.3..0.75).contains(&(giants as f64 / count as f64)), "{giants} giants of {count}");
        assert!((8..45).contains(&first), "{first} of the first magnitude");
        // The farthest are a few kiloparsecs off; a few are dimmed by dust a magnitude or more.
        let far = n.stars.iter().map(|s| s.distance_pc).fold(0.0, f64::max);
        assert!((1500.0..8000.0).contains(&far), "{far} pc");
        assert!(n.stars.iter().filter(|s| s.extinction_mag > 1.0).count() > 50);
        // The plane outshines the pole several times over, and the brightest of the band is
        // toward the centre, not the anticentre (where dust allows, off the plane's lane).
        let galactic = |l: f64, b: f64| {
            let (l, b) = (l.to_radians(), b.to_radians());
            let a = o.axes();
            let local = [b.cos() * l.cos(), b.cos() * l.sin(), b.sin()];
            o.to_system_frame([0, 1, 2].map(|k| a[0][k] * local[0] + a[1][k] * local[1] + a[2][k] * local[2]))
        };
        assert!(n.glow(galactic(90.0, 0.0)) > 3.0, "{}", n.glow(galactic(90.0, 0.0)));
        let brightest = |l0: f64| (-30..=30).flat_map(|l| (-12..=12).map(move |b| (l0 + f64::from(l), f64::from(b)))).map(|(l, b)| n.glow(galactic(l, b))).fold(0.0, f64::max);
        assert!(brightest(0.0) > 1.5 * brightest(180.0), "centre {} vs anticentre {}", brightest(0.0), brightest(180.0));
        // Frames turn both ways.
        let v = [0.3, -0.5, 0.81];
        let back = o.to_galaxy_frame(o.to_system_frame(v));
        assert!((0..3).all(|k| (back[k] - v[k]).abs() < 1e-12));
        let s = &n.stars[0];
        let p = n.galactic(s.direction, s.distance_pc);
        assert!((0..3).all(|k| (p[k] - s.position_pc[k]).abs() < 1e-6));
    }

    /// Inside the bulge the sky is listed to a brighter limit, its count held near four
    /// Suns'.
    #[test]
    fn a_crowded_sky_is_cut() {
        let field = Starfield::new(&Galaxy::milky_way(), Seed(7), StarEdits::default());
        let o = Observer { position_pc: [600.0, 200.0, 50.0], ..sun(3) };
        let n = neighbourhood(&field, &o, NAKED_EYE);
        assert!(n.limit < NAKED_EYE - 1.0, "limit {}", n.limit);
        assert!((5_000..80_000).contains(&n.stars.len()), "{} stars", n.stars.len());
    }

    /// Two systems a few parsecs apart see the same stars, each other among them.
    #[test]
    fn neighbours_see_each_other() {
        let field = Starfield::new(&Galaxy::milky_way(), Seed(7), StarEdits::default());
        let a = field.near([8200.0, 0.0, 25.0], 12.0, 4.0);
        let (x, y) = (&a[0], &a[1]);
        let view = |s: &crate::sectors::GalacticStar| Observer { position_pc: s.position_pc, to_system: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], exclude: Some(s.id), clear_pc: 0.0 };
        let from_x = neighbourhood(&field, &view(x), 7.5);
        let seen = from_x.stars.iter().find(|s| s.id == y.id).expect("y in x's sky");
        assert!((seen.distance_pc - x.distance_pc(y.position_pc)).abs() < 1e-9);
        assert!(!from_x.stars.iter().any(|s| s.id == x.id), "no star sees itself");
    }
}
