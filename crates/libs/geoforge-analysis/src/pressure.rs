//! Sea-level pressure, recovered from the stored surface winds by the balance surface air
//! keeps: pressure's push against the Coriolis turn and friction with the ground,
//! −∇p/ρ = f k̂×v + r·v (Ekman 1905's surface layer, as a damped geostrophic wind). Friction
//! makes winds cross the isobars toward low pressure, by atan(r/f): ~15° over sea and ~35°
//! over land at mid-latitudes, near 90° at the equator where f vanishes. The gradient field
//! this gives is not exactly a gradient, so pressure is its least-squares fit (a Poisson
//! solve on a 2° grid), with the mean set to the world's surface pressure.
//!
//! Derived, never stored: edit the winds and the pressure follows.

use geoforge_core::grid::{Cell, Grid};
use std::f64::consts::PI;

/// Surface friction rate over sea and land, 1/s.
const FRICTION_SEA: f64 = 2.5e-5;
const FRICTION_LAND: f64 = 7e-5;
/// Surface air density at 1 bar, kg/m³.
const AIR_DENSITY: f64 = 1.2;
/// Pressure added over land per kelvin colder than the latitude's mean, and taken away per
/// kelvin warmer, hPa.
const COLD_HPA_PER_K: f64 = 1.5;
const WARM_HPA_PER_K: f64 = 0.6;
/// Highest ground counted in the sea-level reduction, m.
const REDUCTION_MAX_M: f64 = 2000.0;
/// Gas constant of air, J/(kg·K).
const DRY_AIR_J_KG_K: f64 = 287.0;
/// Solve grid, degrees.
const STEP: f64 = 2.0;

/// What the pressure is recovered from, per cell: one season's surface winds (east, north;
/// m/s) and temperatures (°C), and the ground's height (m; land where above 0 on a world
/// with seas, as `land` says).
pub struct Air<'a> {
    pub wind: &'a [[f32; 2]],
    pub temperature: &'a [f64],
    pub elevation: &'a [f32],
    pub land: &'a [bool],
}

/// The body's figures the balance needs.
pub struct Planet {
    pub radius_km: f64,
    pub rotation_h: f64,
    /// Mean surface pressure, hPa.
    pub mean_hpa: f64,
    pub gravity_m_s2: f64,
}

/// Sea-level pressure, hPa, at every cell.
#[must_use]
pub fn sea_level_pressure(grid: &dyn Grid, planet: &Planet, air: &Air) -> Vec<f64> {
    let (wind, land) = (air.wind, air.land);
    // Temperatures brought down to sea level along the standard lapse rate (6.5 K/km at
    // Earth's gravity): a mountain is cold for its height, not a cold dome.
    let lapse = 6.5e-3 * planet.gravity_m_s2 / 9.81;
    let temperature: Vec<f64> = (0..grid.len()).map(|c| air.temperature[c] + lapse * f64::from(air.elevation[c]).max(0.0)).collect();
    let (radius_km, rotation_h, mean_hpa) = (planet.radius_km, planet.rotation_h, planet.mean_hpa);
    let (w, h) = ((360.0 / STEP) as usize, (180.0 / STEP) as usize);
    let lat_of = |j: usize| 90.0 - (j as f64 + 0.5) * STEP;
    // Mean wind and land share per coarse cell.
    let mut sum = vec![(0.0, 0.0, 0.0, 0.0); w * h];
    let mut warmth = vec![0.0; w * h];
    let mut height = vec![0.0; w * h];
    for c in 0..grid.len() {
        let ll = grid.center(c as Cell).to_lat_lon();
        let k = index(w, h, ll.lat_deg, ll.lon_deg);
        let a = grid.area_sr(c as Cell);
        let s = &mut sum[k];
        s.0 += a * f64::from(wind[c][0]);
        s.1 += a * f64::from(wind[c][1]);
        s.2 += a * f64::from(u8::from(land[c]));
        s.3 += a;
        warmth[k] += a * temperature[c];
        height[k] += a * f64::from(air.elevation[c]).max(0.0);
    }
    let omega = 2.0 * PI / (rotation_h.abs().max(1e-3) * 3600.0);
    let rho = AIR_DENSITY * mean_hpa / 1000.0;
    // The pressure gradient each cell's wind implies, Pa/m (east, north).
    let grad: Vec<(f64, f64)> = (0..w * h)
        .map(|k| {
            let s = sum[k];
            if s.3 <= 0.0 {
                return (0.0, 0.0);
            }
            let (u, v, land_share) = (s.0 / s.3, s.1 / s.3, s.2 / s.3);
            let f = 2.0 * omega * lat_of(k / w).to_radians().sin();
            let r = FRICTION_SEA + (FRICTION_LAND - FRICTION_SEA) * land_share;
            // −∇p/ρ = f k̂×v + r v, with k̂×(u, v) = (−v, u).
            (-rho * (-f * v + r * u), -rho * (f * u + r * v))
        })
        .collect();
    let mut p = fit(&grad, w, h, radius_km * 1e3);
    // Over land, cold air is a dense shallow dome and hot air a light one: a season's
    // departure from its latitude's mean (sea-level) temperature adds a high or low that the
    // winds alone underplay; cold domes count for more, as they hug the ground. On high
    // ground, reducing pressure to sea level through a colder column than usual adds
    // p·g·z/(R·T²) per kelvin more (the Siberian high's centre over Mongolia owes it much);
    // far above sea level the reduction means little, so heights count up to 2 km. The sea's
    // departures are small and the winds already carry its lows.
    for j in 0..h {
        let row: Vec<f64> = (0..w).map(|i| warmth[j * w + i] / sum[j * w + i].3.max(1e-30)).collect();
        let zonal = row.iter().sum::<f64>() / w as f64;
        for (i, &t) in row.iter().enumerate() {
            let k = j * w + i;
            let area = sum[k].3.max(1e-30);
            let (land_share, z) = (sum[k].2 / area, (height[k] / area).min(REDUCTION_MAX_M));
            let departure = t - zonal;
            let dome = if departure < 0.0 { COLD_HPA_PER_K } else { WARM_HPA_PER_K };
            let column_k = (t + 273.15).max(150.0);
            let reduction = mean_hpa * planet.gravity_m_s2 * z / (DRY_AIR_J_KG_K * column_k * column_k);
            p[k] -= (dome + reduction) * 100.0 * land_share * departure;
        }
    }
    // Mean over the sphere to the world's surface pressure.
    let weights: Vec<f64> = (0..h).map(|j| lat_of(j).to_radians().cos()).collect();
    let mean: f64 = (0..w * h).map(|k| p[k] * weights[k / w]).sum::<f64>() / (weights.iter().sum::<f64>() * w as f64);
    let coarse: Vec<f64> = p.iter().map(|x| mean_hpa + (x - mean) / 100.0).collect();
    (0..grid.len())
        .map(|c| {
            let ll = grid.center(c as Cell).to_lat_lon();
            bilinear(&coarse, w, h, ll.lat_deg, ll.lon_deg)
        })
        .collect()
}

fn index(w: usize, h: usize, lat: f64, lon: f64) -> usize {
    let j = (((90.0 - lat) / STEP) as usize).min(h - 1);
    let i = (((lon + 180.0).rem_euclid(360.0) / STEP) as usize) % w;
    j * w + i
}

fn bilinear(f: &[f64], w: usize, h: usize, lat: f64, lon: f64) -> f64 {
    let y = ((90.0 - lat) / STEP - 0.5).clamp(0.0, (h - 1) as f64);
    let x = (lon + 180.0) / STEP - 0.5;
    let (j0, fy) = (y.floor() as usize, y.fract());
    let j1 = (j0 + 1).min(h - 1);
    let (i0, fx) = (x.floor(), x - x.floor());
    let at = |j: usize, i: f64| f[j * w + (i as isize).rem_euclid(w as isize) as usize];
    (1.0 - fy) * ((1.0 - fx) * at(j0, i0) + fx * at(j0, i0 + 1.0)) + fy * ((1.0 - fx) * at(j1, i0) + fx * at(j1, i0 + 1.0))
}

/// The pressure (Pa, up to a constant) whose gradient best fits `grad` (Pa/m, east and north)
/// on the coarse grid, in the least-squares sense: Gauss–Seidel with over-relaxation on the
/// face-weighted normal equations.
fn fit(grad: &[(f64, f64)], w: usize, h: usize, radius_m: f64) -> Vec<f64> {
    let d_rad = STEP.to_radians();
    let cos_at = |lat: f64| lat.to_radians().cos().max(1e-6);
    let lat_of = |j: usize| 90.0 - (j as f64 + 0.5) * STEP;
    // Faces of each cell: (neighbour, conductance = face length / distance, expected rise).
    let mut faces: Vec<Vec<(usize, f64, f64)>> = vec![Vec::with_capacity(4); w * h];
    for j in 0..h {
        let lat = lat_of(j);
        let dx = radius_m * cos_at(lat) * d_rad;
        let dy = radius_m * d_rad;
        for i in 0..w {
            let k = j * w + i;
            let east = j * w + (i + 1) % w;
            let ge = 0.5 * (grad[k].0 + grad[east].0);
            let a = dy / dx;
            faces[k].push((east, a, ge * dx));
            faces[east].push((k, a, -ge * dx));
            if j + 1 < h {
                // The cell to the south.
                let south = (j + 1) * w + i;
                let gn = 0.5 * (grad[k].1 + grad[south].1);
                let face = radius_m * cos_at(lat - 0.5 * STEP) * d_rad;
                let a = face / dy;
                faces[k].push((south, a, -gn * dy));
                faces[south].push((k, a, gn * dy));
            }
        }
    }
    let mut p = vec![0.0; w * h];
    let omega = 1.9;
    for _ in 0..4000 {
        let mut change: f64 = 0.0;
        for k in 0..w * h {
            let (mut num, mut den) = (0.0, 0.0);
            for &(nb, a, rise) in &faces[k] {
                num += a * (p[nb] - rise);
                den += a;
            }
            let target = num / den;
            let d = omega * (target - p[k]);
            p[k] += d;
            change = change.max(d.abs());
        }
        if change < 1e-3 {
            break;
        }
    }
    p
}

/// Sea-level pressure in January (`july` false) or July from a world's stored winds,
/// temperatures and relief, hPa per cell; `None` while they are missing.
#[must_use]
pub fn seasonal(products: &geoforge_core::product::Products, grid: &dyn Grid, body: &geoforge_core::BodyContext, july: bool) -> Option<Vec<f64>> {
    let wind = if july { products.wind_jul() } else { products.wind_jan() }?;
    let (t, range) = (&products.temperature()?.values, &products.temperature_range()?.values);
    // July: the north at its warmest month, the south at its coldest.
    let season_t: Vec<f64> = (0..grid.len())
        .map(|c| {
            let north = grid.center(c as Cell).to_lat_lon().lat_deg >= 0.0;
            f64::from(t[c]) + 0.5 * f64::from(range[c]) * if north == july { 1.0 } else { -1.0 }
        })
        .collect();
    let elevation = &products.elevation().or(products.rough_elevation())?.values;
    let wet = body.interior.has_surface_liquid_water();
    let land: Vec<bool> = elevation.iter().map(|&e| !(wet && e < 0.0)).collect();
    let planet = Planet { radius_km: body.exterior.radius_km, rotation_h: body.exterior.rotation_h, mean_hpa: body.interior.surface_pressure_bar() * 1000.0, gravity_m_s2: body.gravity_m_s2() };
    Some(sea_level_pressure(grid, &planet, &Air { wind: &wind.values, temperature: &season_t, elevation, land: &land }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use geoforge_core::grid::GridSpec;

    /// Winds circling a point anticlockwise (northern hemisphere) and spiralling in: a low
    /// there, about as deep as the geostrophic balance says.
    #[test]
    fn a_cyclone_is_a_low() {
        let grid = GridSpec::Equirect { width: 360, height: 180 }.build().unwrap();
        let centre = geoforge_core::UnitVec3::from_lat_lon(45.0, 0.0);
        let wind: Vec<[f32; 2]> = (0..grid.len())
            .map(|c| {
                let p = grid.center(c as Cell);
                let ll = p.to_lat_lon();
                let d = p.angle_to(centre);
                if d > 0.35 {
                    return [0.0, 0.0];
                }
                // Direction to the centre, as east and north components.
                let (dlat, dlon) = (45.0 - ll.lat_deg, -ll.lon_deg * ll.lat_deg.to_radians().cos());
                let n = dlat.hypot(dlon).max(1e-9);
                let (to_e, to_n) = (dlon / n, dlat / n);
                let speed = 10.0 * (d / 0.1).min(1.0) * (1.0 - d / 0.35);
                // Anticlockwise: to the right of the inward direction, plus some inflow.
                let (tan_e, tan_n) = (to_n, -to_e);
                [(speed * (0.9 * tan_e + 0.3 * to_e)) as f32, (speed * (0.9 * tan_n + 0.3 * to_n)) as f32]
            })
            .collect();
        let land = vec![false; grid.len()];
        let planet = Planet { radius_km: 6371.0, rotation_h: 23.934, mean_hpa: 1013.25, gravity_m_s2: 9.81 };
        let p = sea_level_pressure(&grid, &planet, &Air { wind: &wind, temperature: &vec![15.0; grid.len()], elevation: &vec![-1000.0; grid.len()], land: &land });
        let at = |lat: f64, lon: f64| p[grid.locate(geoforge_core::UnitVec3::from_lat_lon(lat, lon)) as usize];
        let (centre_p, far) = (at(45.0, 0.0), at(-45.0, 180.0));
        assert!(far - centre_p > 5.0 && far - centre_p < 60.0, "low {centre_p}, far {far}");
        let mean: f64 = (0..grid.len()).map(|c| p[c] * grid.area_sr(c as Cell)).sum::<f64>() / (4.0 * PI);
        assert!((mean - 1013.25).abs() < 0.5, "{mean}");
    }
}
