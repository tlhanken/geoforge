//! Plate seed placement by electrostatic repulsion on the sphere — a faithful port of
//! geoforge v1 (Stage 1.1). Charges with a power-law size distribution repel until they
//! settle; their positions become plate seeds. Arithmetic is kept identical to v1 so that
//! the same seed yields the same plates.

use geoforge_core::geom::UnitVec3;
use rand::rngs::StdRng;
use rand::Rng;

#[derive(Debug, Clone)]
pub struct Charge {
    pub position: UnitVec3,
    pub charge: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct ElectrostaticConfig {
    pub max_iterations: usize,
    pub convergence_threshold: f64,
    pub damping_factor: f64,
    pub time_step: f64,
    pub min_charge: f64,
    pub max_charge: f64,
}

impl Default for ElectrostaticConfig {
    fn default() -> Self {
        Self {
            max_iterations: 1000,
            convergence_threshold: 1e-6,
            damping_factor: 0.9,
            time_step: 0.01,
            min_charge: 0.01,
            max_charge: 100.0,
        }
    }
}

fn from_lat_lon_v1(lat: f64, lon: f64) -> UnitVec3 {
    UnitVec3::from_lat_lon(lat, lon)
}

/// Random charges, uniformly placed, with v1's power-law charge magnitudes.
pub fn random_charges(n: usize, rng: &mut StdRng, cfg: &ElectrostaticConfig) -> Vec<Charge> {
    let mut charges = Vec::with_capacity(n);
    for i in 0..n {
        let u = rng.r#gen::<f64>();
        let v = rng.r#gen::<f64>();
        let theta = 2.0 * std::f64::consts::PI * u;
        let phi = (2.0 * v - 1.0).acos();
        let lat = (std::f64::consts::PI / 2.0 - phi).to_degrees();
        let lon = theta.to_degrees();
        let lon = if lon > 180.0 { lon - 360.0 } else { lon };
        let position = from_lat_lon_v1(lat, lon);
        let charge = power_law_charge(i, n, rng, cfg);
        charges.push(Charge { position, charge });
    }
    charges
}

fn power_law_charge(index: usize, total: usize, rng: &mut StdRng, cfg: &ElectrostaticConfig) -> f64 {
    let superplate_count = 1;
    let major_plate_count = (total as f64 * 0.15).max(2.0) as usize;
    let medium_plate_count = (total as f64 * 0.25).max(3.0) as usize;
    let charge = if index < superplate_count {
        rng.gen_range(cfg.max_charge * 0.85..=cfg.max_charge)
    } else if index < major_plate_count {
        rng.gen_range(cfg.max_charge * 0.35..=cfg.max_charge * 0.65)
    } else if index < medium_plate_count {
        rng.gen_range(cfg.max_charge * 0.15..=cfg.max_charge * 0.35)
    } else {
        let remaining = total - medium_plate_count;
        let position_in_tail = (index - medium_plate_count) as f64 / remaining as f64;
        let base_charge = cfg.max_charge * 0.15 * (-5.0 * position_in_tail).exp();
        if position_in_tail > 0.8 {
            rng.gen_range(cfg.min_charge..=cfg.min_charge * 2.0)
        } else {
            rng.gen_range(cfg.min_charge..=base_charge.max(cfg.min_charge * 5.0))
        }
    };
    let variation = rng.gen_range(0.9..1.1);
    (charge * variation).clamp(cfg.min_charge, cfg.max_charge)
}

/// Damped repulsion until the largest force drops below the threshold or iterations run out.
pub fn simulate_equilibrium(charges: &mut [Charge], cfg: &ElectrostaticConfig) {
    let mut velocities = vec![(0.0, 0.0, 0.0); charges.len()];
    for _ in 0..cfg.max_iterations {
        let mut max_force: f64 = 0.0;
        let mut total = vec![(0.0, 0.0, 0.0); charges.len()];
        for i in 0..charges.len() {
            for j in 0..charges.len() {
                if i != j {
                    let f = force(&charges[i], &charges[j]);
                    total[i].0 += f.0;
                    total[i].1 += f.1;
                    total[i].2 += f.2;
                }
            }
            let m = (total[i].0.powi(2) + total[i].1.powi(2) + total[i].2.powi(2)).sqrt();
            max_force = max_force.max(m);
        }
        if max_force < cfg.convergence_threshold {
            break;
        }
        for i in 0..charges.len() {
            velocities[i].0 = velocities[i].0 * cfg.damping_factor + total[i].0 * cfg.time_step;
            velocities[i].1 = velocities[i].1 * cfg.damping_factor + total[i].1 * cfg.time_step;
            velocities[i].2 = velocities[i].2 * cfg.damping_factor + total[i].2 * cfg.time_step;
            let p = charges[i].position;
            let nx = p.x + velocities[i].0 * cfg.time_step;
            let ny = p.y + velocities[i].1 * cfg.time_step;
            let nz = p.z + velocities[i].2 * cfg.time_step;
            let norm = (nx * nx + ny * ny + nz * nz).sqrt();
            let (x, y, z) = (nx / norm, ny / norm, nz / norm);
            let lat = z.asin().to_degrees();
            let lon = y.atan2(x).to_degrees();
            charges[i].position = from_lat_lon_v1(lat, lon);
        }
    }
}

fn force(a: &Charge, b: &Charge) -> (f64, f64, f64) {
    let dx = a.position.x - b.position.x;
    let dy = a.position.y - b.position.y;
    let dz = a.position.z - b.position.z;
    let d2 = dx * dx + dy * dy + dz * dz;
    let d = d2.sqrt();
    if d < 1e-10 {
        return (0.0, 0.0, 0.0);
    }
    let m = (a.charge * b.charge) / d2;
    (m * (dx / d), m * (dy / d), m * (dz / d))
}

/// Seed points as v1 stored them (lat/lon round trip included, for bit-identical output).
#[must_use]
pub fn seed_points(charges: &[Charge]) -> Vec<UnitVec3> {
    charges
        .iter()
        .map(|c| {
            let ll = c.position.to_lat_lon();
            UnitVec3::from_lat_lon(ll.lat_deg, ll.lon_deg)
        })
        .collect()
}
