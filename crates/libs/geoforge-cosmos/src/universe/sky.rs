//! Other galaxies as a night sky shows them.
//!
//! Each is drawn as an exponential disk of light on the sky: its half-light radius over 1.678
//! the scale length, inclined by its axis to the line of sight (a disk of thickness t at
//! inclination i looks √(cos²i + t² sin²i) as wide as long), its central surface brightness
//! Σ₀ = L / (2π h² q) in L☉ per pc², which, like the glow's, does not depend on distance. The
//! size given is the diameter to 3.5 scale lengths, about where catalogues' 25 mag/arcsec²
//! isophote falls. The home galaxy's dust is applied by the sky.
//!
//! The Solar System sees the real ones brighter than about 7.5 (V, the integrated light;
//! sizes and position angles from the RC3, distances from the NASA/IPAC Extragalactic
//! Database): the Magellanic Clouds, Andromeda and Triangulum, and M81, Centaurus A, the
//! Sculptor Galaxy and M83.

use super::{Universe, UniverseGalaxy};
use crate::draw::Draw;
use crate::galaxy::Galaxy;
use crate::sky::Observer;
use crate::sol::GALACTIC;
use geoforge_core::body::GalacticSite;
use geoforge_core::galactic::GalaxyId;
use geoforge_core::seed::Seed;

/// Scale lengths to a galaxy's catalogued edge.
const EDGE_SCALES: f64 = 3.5;
/// Galaxies are looked for this far, Mpc.
const REACH_MPC: f64 = 40.0;

/// A galaxy in a sky.
#[derive(Debug, Clone, PartialEq)]
pub struct SkyGalaxy {
    pub id: Option<GalaxyId>,
    pub name: String,
    /// Hubble type.
    pub kind: String,
    /// Toward it, a unit vector in the home galaxy's frame.
    pub direction: [f64; 3],
    pub distance_pc: f64,
    pub absolute_v: f64,
    /// Apparent V, integrated.
    pub v: f64,
    /// Its light's scale length along the major axis, degrees.
    pub scale_deg: f64,
    pub axis_ratio: f64,
    /// Its major axis on the sky: a unit vector square to `direction`.
    pub major: [f64; 3],
    /// Central surface brightness, L☉ (V) per pc².
    pub sigma0: f64,
    /// Its magnitude is as observed, through the home galaxy's dust already.
    pub observed: bool,
}

impl SkyGalaxy {
    /// Surface brightness toward `dir` (a unit vector in the home galaxy's frame), L☉/pc².
    #[must_use]
    pub fn sigma(&self, dir: [f64; 3]) -> f64 {
        let c = dot(dir, self.direction);
        if c <= 0.0 {
            return 0.0;
        }
        let minor = cross(self.direction, self.major);
        let (x, y) = (dot(dir, self.major) / c, dot(dir, minor) / c);
        let r = (x * x + (y / self.axis_ratio).powi(2)).sqrt() / self.scale_deg.to_radians();
        if r > 3.0 * EDGE_SCALES {
            0.0
        } else {
            self.sigma0 * (-r).exp()
        }
    }

    /// Diameter along the major axis to its catalogued edge, degrees.
    #[must_use]
    pub fn size_deg(&self) -> f64 {
        2.0 * EDGE_SCALES * self.scale_deg
    }

    /// Dimmed by `mag` of the home galaxy's dust.
    #[must_use]
    pub fn dimmed(mut self, mag: f64) -> Self {
        if !self.observed {
            self.v += mag;
            self.sigma0 *= 10f64.powf(-0.4 * mag);
        }
        self
    }
}

/// The galaxies a system at `site` (its galaxy `home_model`, seen from `observer`) sees
/// brighter than `v_limit` before dust, brightest first.
#[must_use]
pub fn seen(site: &GalacticSite, system_seed: Seed, home_model: &Galaxy, observer: &Observer, v_limit: f64) -> Vec<SkyGalaxy> {
    let mut out = match Universe::of_site(site, system_seed) {
        None => real(observer),
        Some(u) => generated(&u, home_model, observer, v_limit),
    };
    out.retain(|g| g.v <= v_limit);
    out.sort_by(|a, b| a.v.total_cmp(&b.v));
    out
}

fn generated(u: &Universe, home_model: &Galaxy, observer: &Observer, v_limit: f64) -> Vec<SkyGalaxy> {
    let home = u.home(home_model);
    // The home galaxy's frame in the universe's: its axis is the north pole.
    let z = home.axis;
    let mut d = Draw::new(home.seed.child(b"frame"));
    let x = unit(cross(z, d.direction()));
    let frame = [x, cross(z, x), z];
    u.galaxies_seen(home.position_mpc, REACH_MPC, v_limit + 1.0).into_iter().filter(|g| g.id != home.id).filter_map(|g| place(&g, &home, &frame, observer)).collect()
}

fn place(g: &UniverseGalaxy, home: &UniverseGalaxy, frame: &[[f64; 3]; 3], observer: &Observer) -> Option<SkyGalaxy> {
    let rel = [0, 1, 2].map(|k| (g.position_mpc[k] - home.position_mpc[k]) * 1e6);
    let p = frame.map(|row| dot(row, rel));
    let v = [0, 1, 2].map(|k| p[k] - observer.position_pc[k]);
    let distance_pc = norm(v);
    if distance_pc < 1.0 {
        return None;
    }
    let direction = v.map(|c| c / distance_pc);
    let axis = frame.map(|row| dot(row, g.axis));
    let cos_i = dot(axis, direction).abs().min(1.0);
    let axis_ratio = (cos_i * cos_i + g.thickness.powi(2) * (1.0 - cos_i * cos_i)).sqrt();
    let on_sky = [0, 1, 2].map(|k| axis[k] - dot(axis, direction) * direction[k]);
    let major = if norm(on_sky) > 1e-6 { unit(cross(direction, on_sky)) } else { unit(cross(direction, if direction[2].abs() < 0.9 { [0.0, 0.0, 1.0] } else { [1.0, 0.0, 0.0] })) };
    let h_pc = g.half_light_kpc * 1000.0 / 1.678;
    let absolute_v = g.absolute_v();
    Some(SkyGalaxy {
        id: Some(g.id),
        name: g.name(),
        kind: g.kind(),
        direction,
        distance_pc,
        absolute_v,
        v: absolute_v + 5.0 * (distance_pc / 10.0).log10(),
        scale_deg: (h_pc / distance_pc).to_degrees(),
        axis_ratio,
        major,
        sigma0: g.luminosity_v() / (std::f64::consts::TAU * h_pc * h_pc * axis_ratio),
        observed: false,
    })
}

/// A real galaxy as the Sun sees it (J2000).
struct Real {
    name: &'static str,
    kind: &'static str,
    ra_deg: f64,
    dec_deg: f64,
    kpc: f64,
    v: f64,
    /// Major and minor axes, arcmin.
    major_arcmin: f64,
    minor_arcmin: f64,
    /// Position angle, degrees east of north.
    pa_deg: f64,
}

/// The real galaxies the Sun sees to about magnitude 7.5.
const REAL: [Real; 8] = [
    Real { name: "Large Magellanic Cloud", kind: "SBm", ra_deg: 80.894, dec_deg: -69.756, kpc: 49.6, v: 0.9, major_arcmin: 645.0, minor_arcmin: 550.0, pa_deg: 170.0 },
    Real { name: "Small Magellanic Cloud", kind: "dIrr", ra_deg: 13.187, dec_deg: -72.829, kpc: 62.4, v: 2.7, major_arcmin: 320.0, minor_arcmin: 185.0, pa_deg: 45.0 },
    Real { name: "Andromeda Galaxy", kind: "SAb", ra_deg: 10.685, dec_deg: 41.269, kpc: 765.0, v: 3.44, major_arcmin: 190.0, minor_arcmin: 60.0, pa_deg: 35.0 },
    Real { name: "Triangulum Galaxy", kind: "SAcd", ra_deg: 23.462, dec_deg: 30.660, kpc: 840.0, v: 5.72, major_arcmin: 70.8, minor_arcmin: 41.7, pa_deg: 23.0 },
    Real { name: "Centaurus A", kind: "S0 pec", ra_deg: 201.365, dec_deg: -43.019, kpc: 3800.0, v: 6.84, major_arcmin: 25.7, minor_arcmin: 20.0, pa_deg: 35.0 },
    Real { name: "Bode's Galaxy", kind: "SAab", ra_deg: 148.888, dec_deg: 69.065, kpc: 3600.0, v: 6.94, major_arcmin: 26.9, minor_arcmin: 14.1, pa_deg: 157.0 },
    Real { name: "Sculptor Galaxy", kind: "SABc", ra_deg: 11.888, dec_deg: -25.288, kpc: 3500.0, v: 7.1, major_arcmin: 27.5, minor_arcmin: 6.8, pa_deg: 52.0 },
    Real { name: "Southern Pinwheel", kind: "SABc", ra_deg: 204.254, dec_deg: -29.866, kpc: 4600.0, v: 7.54, major_arcmin: 12.9, minor_arcmin: 11.5, pa_deg: 45.0 },
];

fn real(observer: &Observer) -> Vec<SkyGalaxy> {
    let axes = observer.axes();
    // ICRS → the Sun's galactic axes → the galaxy's frame.
    let to_galaxy = |u: [f64; 3]| {
        let local = GALACTIC.map(|g| dot(g, u));
        [0, 1, 2].map(|k| local[0] * axes[0][k] + local[1] * axes[1][k] + local[2] * axes[2][k])
    };
    REAL.iter()
        .map(|&Real { name, kind, ra_deg: ra, dec_deg: dec, kpc, v, major_arcmin: major, minor_arcmin: minor, pa_deg: pa }| {
            let (a, dl) = (ra.to_radians(), dec.to_radians());
            let toward = [dl.cos() * a.cos(), dl.cos() * a.sin(), dl.sin()];
            let north = [-dl.sin() * a.cos(), -dl.sin() * a.sin(), dl.cos()];
            let east = [-a.sin(), a.cos(), 0.0];
            let (s, c) = pa.to_radians().sin_cos();
            let along = [0, 1, 2].map(|k| c * north[k] + s * east[k]);
            let distance_pc = kpc * 1000.0;
            let absolute_v = v - 5.0 * (distance_pc / 10.0).log10();
            let scale_deg = major / 60.0 / 2.0 / EDGE_SCALES;
            let axis_ratio = minor / major;
            let h_pc = scale_deg.to_radians() * distance_pc;
            let luminosity = 10f64.powf(-0.4 * (absolute_v - 4.83));
            SkyGalaxy {
                id: None,
                name: name.into(),
                kind: kind.into(),
                direction: to_galaxy(toward),
                distance_pc,
                absolute_v,
                v,
                scale_deg,
                axis_ratio,
                major: to_galaxy(along),
                sigma0: luminosity / (std::f64::consts::TAU * h_pc * h_pc * axis_ratio),
                observed: true,
            }
        })
        .collect()
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

fn unit(a: [f64; 3]) -> [f64; 3] {
    let n = norm(a).max(1e-300);
    a.map(|c| c / n)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sol;

    fn sun() -> (GalacticSite, Seed, Observer) {
        let system = sol::at(Seed(1), 2_451_545.0);
        let site = system.root.site.clone().unwrap();
        (site, system.root.seed, Observer::of_system(&system, 8.2))
    }

    /// The real galaxies come out where they are: the LMC in Dorado near the south
    /// galactic pole's side, Andromeda 21° below the plane toward l = 121°.
    #[test]
    fn the_sun_sees_the_real_neighbours() {
        let (site, seed, observer) = sun();
        let g = seen(&site, seed, &Galaxy::milky_way(), &observer, 6.5);
        assert_eq!(g.iter().map(|g| g.name.as_str()).collect::<Vec<_>>(), ["Large Magellanic Cloud", "Small Magellanic Cloud", "Andromeda Galaxy", "Triangulum Galaxy"]);
        let lb = |d: [f64; 3]| {
            let local = observer.axes().map(|ax| dot(ax, d));
            (local[1].atan2(local[0]).to_degrees().rem_euclid(360.0), local[2].asin().to_degrees())
        };
        let (l, b) = lb(g[2].direction);
        assert!((l - 121.17).abs() < 0.3 && (b + 21.57).abs() < 0.3, "Andromeda at l {l}, b {b}");
        let (l, b) = lb(g[0].direction);
        assert!((l - 280.47).abs() < 0.3 && (b + 32.89).abs() < 0.3, "LMC at l {l}, b {b}");
        // Andromeda's disk is ~22 mag/arcsec² on average inside its edge, brighter at the
        // centre: Σ₀ of hundreds of L☉/pc².
        assert!((100.0..2000.0).contains(&g[2].sigma0), "Andromeda Σ₀ {}", g[2].sigma0);
        assert!(g[2].sigma(g[2].direction) > 10.0 * g[2].sigma(rotate_toward(g[2].direction, g[2].major, 1.5)));
    }

    fn rotate_toward(d: [f64; 3], along: [f64; 3], deg: f64) -> [f64; 3] {
        let (s, c) = deg.to_radians().sin_cos();
        unit([0, 1, 2].map(|k| c * d[k] + s * along[k]))
    }

    /// A generated Milky Way: a naked-eye galaxy or a few (its satellites, a neighbour), the
    /// brightest seldom brighter than the LMC.
    #[test]
    fn generated_skies_hold_a_few_galaxies() {
        let mut counts = Vec::new();
        for k in 0..24 {
            let mut site = sun().0;
            site.universe = Some(geoforge_core::galactic::UniverseSite::Seeded(Seed(100 + k)));
            let g = seen(&site, Seed(k), &Galaxy::milky_way(), &sun().2, 6.5);
            assert!(g.iter().all(|g| g.v <= 6.5 && g.distance_pc > 1e3));
            counts.push(g.len());
        }
        let mean = counts.iter().sum::<usize>() as f64 / counts.len() as f64;
        assert!((0.5..8.0).contains(&mean), "naked-eye galaxies per sky: {counts:?}");
    }
}
