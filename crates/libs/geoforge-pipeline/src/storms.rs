//! Storms across the map: seeded lows with their fronts, and tropical cyclones, born where
//! the climate makes them, steered by the season's winds, deepening and filling. The same
//! storms pass over every place, so neighbours share their weather.
//!
//! - **Mid-latitude lows** form along strong north–south temperature contrasts
//!   (baroclinic zones) at 25–70°, favouring the sea, about 2.5 a day per hemisphere on an
//!   Earth-sized world (the order of the lows tracked in reanalyses; Hoskins & Hodges
//!   2002), more where the contrast is stronger (winter) and fewer on slow rotators. They
//!   live 3–7 days, deepen to 10–35 hPa below their surroundings over 1000–1500 km, so a
//!   mid-latitude place is under a storm's rain on about a third of days, and drift east at
//!   10–15 m/s and
//!   slightly poleward: the steering of the upper westerlies, here 1.5 times the surface
//!   wind plus 5 m/s eastward. Each carries a warm front (to the east-southeast in the north,
//!   east-northeast in the south) and a trailing cold front (south-southwest, or north-
//!   northwest), with the warm sector between them; rain falls broadly ahead of the warm
//!   front, in a narrow band on the cold front and around the centre; the warm sector is
//!   mild and the air behind the cold front cold.
//! - **Tropical cyclones** form over sea above 26.5 °C at 5–25° (Earth's rate per area of
//!   such sea: ~85 a year), drift west and poleward with the trades, recurve into the
//!   westerlies, deepen to 15–70 hPa over ~250 km, and weaken fast over land or cool sea.
//! - **Mountains** steer storms: a low heading for high ground turns along it, choosing the
//!   way round that costs least (a 90° turn weighs as much as 800 m of rock ahead), and fills
//!   while over ground above 1500 m; lows form more readily in the lee of ranges, as the
//!   Colorado and Genoa lows do (genesis up to doubled 400 km downwind of 2 km of relief).
//! - **Winds** are gradient winds about each low (Coriolis, pressure gradient and
//!   curvature), cut by 30% for surface friction and turned 20° inward.
//!
//! A storm's track is worked out from its birth, which is seeded by its day, so the storms
//! on a date are the same however they are asked for.

use geoforge_core::grid::Grid;
use geoforge_core::product::Products;
use geoforge_core::seed::{unit, Seed};
use geoforge_core::{BodyContext, UnitVec3};
use serde::Serialize;

/// The longest a storm lives, days.
pub const MAX_LIFE_DAYS: f64 = 10.0;
/// Track step, days.
const STEP_DAYS: f64 = 0.25;
/// The storm climate's grid, degrees.
const COARSE_DEG: f64 = 2.0;
/// Closest two lows are born on the same day, km.
const SPACING_KM: f64 = 2000.0;
/// Air density near the ground, kg/m³.
const AIR_DENSITY: f64 = 1.2;
/// Earth's rotation rate, rad/s.
const EARTH_OMEGA: f64 = 7.292e-5;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum StormKind {
    /// A mid-latitude low with fronts.
    Low,
    TropicalCyclone,
}

/// Where a storm is at a moment, and how strong.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct Fix {
    pub lat: f64,
    pub lon: f64,
    /// Central pressure below its surroundings, hPa.
    pub depth_hpa: f64,
    /// Size of the pressure dip (e-folding), km.
    pub radius_km: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Storm {
    pub name: String,
    pub kind: StormKind,
    /// Earth days from the present.
    pub born_day: f64,
    /// Fixes every quarter day from birth.
    pub track: Vec<Fix>,
}

impl Storm {
    /// Where the storm is on day `t`, if alive.
    #[must_use]
    pub fn at(&self, t: f64) -> Option<Fix> {
        let x = (t - self.born_day) / STEP_DAYS;
        if x < 0.0 || x > (self.track.len() - 1) as f64 {
            return None;
        }
        let (i, f) = (x.floor() as usize, x.fract());
        let (a, b) = (self.track[i], self.track[(i + 1).min(self.track.len() - 1)]);
        let dlon = wrap(b.lon - a.lon);
        Some(Fix { lat: a.lat + f * (b.lat - a.lat), lon: wrap(a.lon + f * dlon), depth_hpa: a.depth_hpa + f * (b.depth_hpa - a.depth_hpa), radius_km: a.radius_km + f * (b.radius_km - a.radius_km) })
    }

    /// The warm and cold fronts at `fix`, as lines of (lat, lon) from the centre outward.
    #[must_use]
    pub fn fronts(&self, fix: Fix, planet_radius_km: f64) -> Option<[Vec<(f64, f64)>; 2]> {
        if self.kind != StormKind::Low || fix.depth_hpa < 4.0 {
            return None;
        }
        let line = |bearing: f64, length: f64| (0..=12).map(|k| offset(fix, bearing, length * fix.radius_km * f64::from(k) / 12.0, planet_radius_km)).collect();
        Some([line(WARM_FRONT_DEG, WARM_FRONT_LENGTH), line(COLD_FRONT_DEG, COLD_FRONT_LENGTH)])
    }
}

/// Bearings of the fronts from the centre, degrees clockwise from poleward, and their
/// lengths in storm radii.
const WARM_FRONT_DEG: f64 = 115.0;
const COLD_FRONT_DEG: f64 = 200.0;
const WARM_FRONT_LENGTH: f64 = 1.2;
const COLD_FRONT_LENGTH: f64 = 1.6;

/// What the storms do at a place at a moment.
#[derive(Debug, Clone, PartialEq, Default, Serialize)]
pub struct Effect {
    /// How strongly storms favour rain here, 0–1.
    pub rain: f64,
    /// Warm or cold air brought in, in units of the place's day-to-day spread.
    pub warmth: f64,
    /// Storm winds (east, north), m/s.
    pub wind: [f64; 2],
    /// Pressure change, hPa.
    pub pressure_hpa: f64,
    /// The storm that matters most here, if any.
    pub storm: Option<String>,
    pub kind: Option<StormKind>,
}

/// The climate storms are drawn from, on a coarse grid, and the planet's spin.
pub struct Atmosphere {
    rows: usize,
    cols: usize,
    mean_c: Vec<f64>,
    range_c: Vec<f64>,
    land: Vec<bool>,
    /// Highest ground in each coarse cell, m (0 at sea).
    height_m: Vec<f64>,
    wind_jan: Vec<[f64; 2]>,
    wind_jul: Vec<[f64; 2]>,
    pub radius_km: f64,
    /// Rotation rate, rad/s.
    pub omega: f64,
    year: f64,
    season: Box<dyn Fn(f64) -> f64>,
    seed: Seed,
}

/// A place's climate as storms see it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Local {
    /// Annual mean and warmest-to-coldest range, °C.
    pub mean_c: f64,
    pub range_c: f64,
    pub land: bool,
    /// The highest ground around, m (0 at sea).
    pub height_m: f64,
    /// January and July winds (east, north), m/s.
    pub wind_jan: [f64; 2],
    pub wind_jul: [f64; 2],
}

impl Atmosphere {
    /// From a climate given as a function of (lat, lon); `season(t)` is the northern season
    /// (−1 to 1) on Earth day `t`.
    pub fn from_fn(radius_km: f64, rotation_h: f64, year: f64, season: Box<dyn Fn(f64) -> f64>, seed: Seed, climate: impl Fn(f64, f64) -> Local) -> Self {
        let (rows, cols) = ((180.0 / COARSE_DEG) as usize, (360.0 / COARSE_DEG) as usize);
        let mut a = Self {
            rows,
            cols,
            mean_c: Vec::with_capacity(rows * cols),
            range_c: Vec::with_capacity(rows * cols),
            land: Vec::with_capacity(rows * cols),
            height_m: Vec::with_capacity(rows * cols),
            wind_jan: Vec::with_capacity(rows * cols),
            wind_jul: Vec::with_capacity(rows * cols),
            radius_km,
            omega: std::f64::consts::TAU / (rotation_h.abs().max(1.0) * 3600.0),
            year,
            season,
            seed,
        };
        for i in 0..rows * cols {
            let (lat, lon) = a.centre(i);
            let l = climate(lat, lon);
            a.mean_c.push(l.mean_c);
            a.range_c.push(l.range_c);
            a.land.push(l.land);
            a.height_m.push(l.height_m.max(0.0));
            a.wind_jan.push(l.wind_jan);
            a.wind_jul.push(l.wind_jul);
        }
        a
    }

    /// From a world's climate products; `None` if the climate is not generated.
    #[must_use]
    pub fn read(p: &Products, grid: &dyn Grid, body: &BodyContext, year: f64, season: Box<dyn Fn(f64) -> f64>, seed: Seed) -> Option<Self> {
        let (t, r, e) = (&p.temperature()?.values, &p.temperature_range()?.values, &p.elevation()?.values);
        let (jan, jul) = (p.wind_jan().map(|f| &f.values), p.wind_jul().map(|f| &f.values));
        let sea_world = body.interior.has_surface_liquid_water();
        let wind = |w: Option<&Vec<[f32; 2]>>, c: usize| w.map_or([0.0; 2], |v| v[c].map(f64::from));
        Some(Self::from_fn(body.exterior.radius_km, body.exterior.rotation_h, year, season, seed, |lat, lon| {
            let at = |la: f64, lo: f64| grid.locate(UnitVec3::from_lat_lon(la, lo)) as usize;
            let c = at(lat, lon);
            // The barrier a storm meets: the highest of a 4×4 sampling of the cell.
            let offsets = [-0.375, -0.125, 0.125, 0.375].map(|f| f * COARSE_DEG);
            let height_m = offsets.iter().flat_map(|&a| offsets.iter().map(move |&b| (a, b))).map(|(a, b)| f64::from(e[at((lat + a).clamp(-90.0, 90.0), lon + b)])).fold(0.0, f64::max);
            Local {
                mean_c: f64::from(t[c]),
                range_c: f64::from(r[c]),
                land: !sea_world || e[c] >= 0.0,
                height_m: if sea_world { height_m } else { 0.0 },
                wind_jan: wind(jan, c),
                wind_jul: wind(jul, c),
            }
        }))
    }

    fn centre(&self, i: usize) -> (f64, f64) {
        (90.0 - (((i / self.cols) as f64) + 0.5) * COARSE_DEG, -180.0 + (((i % self.cols) as f64) + 0.5) * COARSE_DEG)
    }

    fn cell(&self, lat: f64, lon: f64) -> usize {
        let row = (((90.0 - lat) / COARSE_DEG) as usize).min(self.rows - 1);
        let col = (((wrap(lon) + 180.0) / COARSE_DEG) as usize).min(self.cols - 1);
        row * self.cols + col
    }

    /// The northern season felt on day `t`: a month behind the star.
    fn felt(&self, t: f64) -> f64 {
        (self.season)(t - self.year / 12.0)
    }

    /// Air temperature at cell `i` in northern season `s`, °C.
    fn temperature(&self, i: usize, s: f64) -> f64 {
        let lat = self.centre(i).0;
        self.mean_c[i] + self.range_c[i] / 2.0 * s * (lat / 10.0).clamp(-1.0, 1.0)
    }

    /// The prevailing surface wind at (lat, lon) in northern season `s`.
    fn wind(&self, lat: f64, lon: f64, s: f64) -> [f64; 2] {
        let i = self.cell(lat, lon);
        let mix = 0.5 * (1.0 + s);
        [0, 1].map(|k| self.wind_jan[i][k] * (1.0 - mix) + self.wind_jul[i][k] * mix)
    }

    /// How much more readily lows form at cell `i` for the mountains upwind: up to double,
    /// 400 km downwind of 2 km of relief.
    fn lee(&self, i: usize, season: f64) -> f64 {
        let (lat, lon) = self.centre(i);
        let w = self.wind(lat, lon, season);
        let speed = w[0].hypot(w[1]);
        if speed < 0.5 {
            return 1.0;
        }
        let (up_lat, up_lon) = destination(lat, lon, (-w[0]).atan2(-w[1]).to_degrees(), 400.0, self.radius_km);
        let relief = (self.height(up_lat, up_lon) - self.height(lat, lon)).max(0.0);
        1.0 + (relief / 2000.0).min(1.0)
    }

    fn height(&self, lat: f64, lon: f64) -> f64 {
        self.height_m[self.cell(lat, lon)]
    }

    /// Turn a velocity away from high ground ahead: among headings up to 90° either side,
    /// the one with the least rock on the next 600 km, a 90° turn costing as much as 800 m.
    fn steer(&self, lat: f64, lon: f64, v: [f64; 2]) -> [f64; 2] {
        let speed = v[0].hypot(v[1]);
        if speed < 0.1 {
            return v;
        }
        let heading = v[0].atan2(v[1]).to_degrees();
        let barrier = |dev: f64| {
            let h = heading + dev;
            let (a, b) = (destination(lat, lon, h, 300.0, self.radius_km), destination(lat, lon, h, 600.0, self.radius_km));
            self.height(a.0, a.1).max(self.height(b.0, b.1)) / 1000.0 + dev.abs() / 90.0 * 0.8
        };
        let dev = [0.0, -30.0, 30.0, -60.0, 60.0, -90.0, 90.0].into_iter().min_by(|a, b| barrier(*a).total_cmp(&barrier(*b))).unwrap_or(0.0);
        let (s, c) = (heading + dev).to_radians().sin_cos();
        let slow = 1.0 - 0.3 * dev.abs() / 90.0;
        [speed * slow * s, speed * slow * c]
    }

    /// Storms alive at any time from `from_day` to `to_day`.
    #[must_use]
    pub fn storms(&self, from_day: f64, to_day: f64) -> Vec<Storm> {
        let first = (from_day - MAX_LIFE_DAYS).floor() as i64;
        (first..=to_day.floor() as i64).flat_map(|d| self.born_on(d)).filter(|s| s.born_day + (s.track.len() - 1) as f64 * STEP_DAYS >= from_day).collect()
    }

    /// The storms born on day `d`.
    fn born_on(&self, d: i64) -> Vec<Storm> {
        let seed = self.seed.child(b"storms").child_u64(b"day", d as u64);
        let mut n = 0u64;
        let mut u = || {
            n += 1;
            unit(seed.child_u64(b"draw", n))
        };
        let t = d as f64;
        let season = self.felt(t);
        let n_cells = self.rows * self.cols;
        let km = self.radius_km * COARSE_DEG.to_radians();
        let area = |i: usize| self.centre(i).0.to_radians().cos() * km * km;
        // Baroclinicity: north–south temperature difference per 1000 km.
        let gradient = |i: usize| {
            let (r, c) = (i / self.cols, i % self.cols);
            let (up, down) = (r.saturating_sub(1) * self.cols + c, (r + 1).min(self.rows - 1) * self.cols + c);
            (self.temperature(up, season) - self.temperature(down, season)).abs() / (km * (down / self.cols - up / self.cols) as f64) * 1000.0
        };
        let mut out = Vec::new();
        let rotation = (self.omega / EARTH_OMEGA).min(1.0);
        let size = (self.radius_km / 6371.0).powi(2);
        for north in [true, false] {
            let in_band = |i: usize| {
                let lat = self.centre(i).0;
                (lat > 0.0) == north && (25.0..=70.0).contains(&lat.abs())
            };
            let weights: Vec<f64> = (0..n_cells).map(|i| if in_band(i) { (gradient(i) - 2.0).max(0.0) * if self.land[i] { 1.0 } else { 1.5 } * self.lee(i, season) * area(i) } else { 0.0 }).collect();
            let band: f64 = (0..n_cells).filter(|&i| in_band(i)).map(area).sum();
            let mean_gradient = (0..n_cells).filter(|&i| in_band(i)).map(|i| gradient(i) * area(i)).sum::<f64>() / band.max(1.0);
            let rate = 2.5 * size * rotation * (mean_gradient / 5.0).clamp(0.0, 3.0);
            for _ in 0..poisson(rate, &mut u) {
                if let Some(i) = pick(&weights, u()) {
                    // Lows born the same day keep apart.
                    let (lat, lon) = self.centre(i);
                    if out.iter().all(|s: &Storm| great_circle_km(s.track[0].lat, s.track[0].lon, lat, lon, self.radius_km) > SPACING_KM) {
                        out.push(self.track(StormKind::Low, i, t + u(), &mut u));
                    }
                }
            }
        }
        // Tropical cyclones over warm sea.
        let warm: Vec<f64> = (0..n_cells)
            .map(|i| {
                let lat = self.centre(i).0.abs();
                let sst = self.temperature(i, season);
                if self.land[i] || !(5.0..=25.0).contains(&lat) || sst < 26.5 {
                    return 0.0;
                }
                (sst - 26.5 + 0.5) * area(i)
            })
            .collect();
        let warm_area: f64 = (0..n_cells).filter(|&i| warm[i] > 0.0).map(area).sum();
        for _ in 0..poisson(85.0 / 365.25 * warm_area / 6e7 * rotation, &mut u) {
            if let Some(i) = pick(&warm, u()) {
                out.push(self.track(StormKind::TropicalCyclone, i, t + u(), &mut u));
            }
        }
        for (k, s) in out.iter_mut().enumerate() {
            s.name = geoforge_names::name(seed.child_u64(b"name", k as u64).0);
        }
        out
    }

    /// A storm born at cell `i` at time `t`, followed through its life.
    fn track(&self, kind: StormKind, i: usize, t: f64, u: &mut impl FnMut() -> f64) -> Storm {
        let (mut lat, mut lon) = self.centre(i);
        lat += (u() - 0.5) * COARSE_DEG;
        lon += (u() - 0.5) * COARSE_DEG;
        let (life, peak, radius) = match kind {
            StormKind::Low => (3.0 + 4.0 * u(), 10.0 + 25.0 * u(), 1000.0 + 500.0 * u()),
            StormKind::TropicalCyclone => (4.0 + 6.0 * u(), 15.0 + 55.0 * u() * u(), 200.0 + 100.0 * u()),
        };
        let steps = (life / STEP_DAYS).ceil() as usize;
        let mut track = Vec::with_capacity(steps + 1);
        let mut weakening = 1.0;
        for s in 0..=steps {
            let age = s as f64 * STEP_DAYS;
            let now = t + age;
            let depth = peak * (std::f64::consts::PI * age / life).sin().max(0.0) * weakening;
            track.push(Fix { lat, lon, depth_hpa: depth, radius_km: radius });
            let season = self.felt(now);
            let w = self.wind(lat, lon, season);
            let pole = lat.signum();
            let v = match kind {
                StormKind::Low => [1.5 * w[0] + 5.0, 1.5 * w[1] + 1.5 * pole],
                // In the westerlies a tropical cyclone recurves and speeds up.
                StormKind::TropicalCyclone if lat.abs() > 28.0 => [1.2 * w[0] + 4.0, 1.2 * w[1] + 2.5 * pole],
                StormKind::TropicalCyclone => [0.8 * w[0], 0.8 * w[1] + 2.0 * pole],
            };
            let v = self.steer(lat, lon, v);
            let step_km = STEP_DAYS * 86.4;
            lat = (lat + (v[1] * step_km / self.radius_km).to_degrees()).clamp(-85.0, 85.0);
            lon = wrap(lon + (v[0] * step_km / (self.radius_km * lat.to_radians().cos().max(0.1))).to_degrees());
            // Tropical cyclones die over land or cool sea, lows away from the contrast that
            // feeds them.
            let c = self.cell(lat, lon);
            let starving = match kind {
                StormKind::TropicalCyclone => self.land[c] || self.temperature(c, season) < 25.0,
                StormKind::Low => lat.abs() < 25.0 || self.height_m[c] > 1500.0,
            };
            if starving {
                weakening *= 0.5f64.powf(STEP_DAYS / 0.75);
            }
        }
        Storm { name: String::new(), kind, born_day: t, track }
    }
}

/// The storms' effect at (`lat`, `lon`) on day `t`, on a planet `radius_km` turning at
/// `omega`.
#[must_use]
pub fn effect(storms: &[Storm], lat: f64, lon: f64, t: f64, radius_km: f64, omega: f64) -> Effect {
    let mut e = Effect::default();
    let (mut dry, mut strongest) = (1.0, 0.0);
    let hem = if lat >= 0.0 { 1.0 } else { -1.0 };
    for s in storms {
        let Some(fix) = s.at(t) else { continue };
        if fix.depth_hpa <= 0.5 {
            continue;
        }
        let (x, y) = local_km(fix, lat, lon, radius_km);
        let r = x.hypot(y).max(1.0);
        let q = r / fix.radius_km;
        if q > 3.0 {
            continue;
        }
        let gauss = (-q * q).exp();
        e.pressure_hpa -= fix.depth_hpa * gauss;
        // Gradient wind about the low, then friction and inflow.
        let (r_m, big_r) = (r * 1000.0, fix.radius_km * 1000.0);
        let dpdr = fix.depth_hpa * 100.0 * 2.0 * r_m / (big_r * big_r) * gauss;
        let f = 2.0 * omega * fix.lat.to_radians().sin().abs().max(5f64.to_radians().sin());
        let speed = 0.7 * (-f * r_m / 2.0 + ((f * r_m / 2.0).powi(2) + r_m / AIR_DENSITY * dpdr).sqrt());
        let turn = if fix.lat >= 0.0 { 1.0 } else { -1.0 };
        let (tx, ty) = (-y / r * turn, x / r * turn);
        let (inflow, along) = (20f64.to_radians().sin(), 20f64.to_radians().cos());
        let wind = [speed * (tx * along - x / r * inflow), speed * (ty * along - y / r * inflow)];
        e.wind = [e.wind[0] + wind[0], e.wind[1] + wind[1]];
        let strength = (fix.depth_hpa / 25.0).clamp(0.3, 1.4);
        let rain = match s.kind {
            StormKind::TropicalCyclone => (-(q / 1.5).powi(2)).exp() * (fix.depth_hpa / 30.0).min(1.0),
            StormKind::Low => {
                // Bearing from the centre, clockwise from poleward (mirrored in the south).
                let ys = y * hem;
                let bearing = x.atan2(ys).to_degrees().rem_euclid(360.0);
                let to_front = |b: f64, len: f64| segment_distance(x, ys, b, len * fix.radius_km);
                let (warm_front, cold_front) = (to_front(WARM_FRONT_DEG, WARM_FRONT_LENGTH), to_front(COLD_FRONT_DEG, COLD_FRONT_LENGTH));
                let warm_sector = (WARM_FRONT_DEG..COLD_FRONT_DEG).contains(&bearing) && q < 1.5;
                let behind = !warm_sector && (COLD_FRONT_DEG..360.0).contains(&bearing);
                let spread = (-(q / 1.5).powi(2)).exp() * strength;
                e.warmth += if warm_sector {
                    1.2 * spread
                } else if behind {
                    -spread
                } else {
                    -0.3 * spread
                };
                let ahead = if !(WARM_FRONT_DEG..=330.0).contains(&bearing) { (-(warm_front / 400.0).powi(2)).exp() * 0.75 } else { (-(warm_front / 150.0).powi(2)).exp() * 0.6 };
                let centre = 0.85 * (-(q / 0.6).powi(2)).exp();
                let cold = 0.85 * (-(cold_front / 200.0).powi(2)).exp();
                let sector = if warm_sector {
                    0.25
                } else if behind {
                    0.15 * gauss
                } else {
                    0.0
                };
                (centre.max(ahead).max(cold).max(sector) * strength).min(1.0)
            }
        };
        dry *= 1.0 - rain;
        let weight = speed + 20.0 * rain;
        if weight > strongest && (rain > 0.3 || speed > 12.0) {
            strongest = weight;
            e.storm = Some(s.name.clone());
            e.kind = Some(s.kind);
        }
    }
    e.rain = 1.0 - dry;
    e
}

/// East and north offsets (km) of (lat, lon) from a storm's centre.
fn local_km(fix: Fix, lat: f64, lon: f64, radius_km: f64) -> (f64, f64) {
    let mid = (0.5 * (lat + fix.lat)).to_radians().cos();
    (wrap(lon - fix.lon).to_radians() * mid * radius_km, (lat - fix.lat).to_radians() * radius_km)
}

/// A point `dist` km from a storm's centre along `bearing` (degrees clockwise from poleward).
fn offset(fix: Fix, bearing: f64, dist: f64, radius_km: f64) -> (f64, f64) {
    let hem = if fix.lat >= 0.0 { 1.0 } else { -1.0 };
    let (s, c) = bearing.to_radians().sin_cos();
    let lat = fix.lat + (dist * c * hem / radius_km).to_degrees();
    let lon = fix.lon + (dist * s / (radius_km * fix.lat.to_radians().cos().max(0.1))).to_degrees();
    (lat.clamp(-90.0, 90.0), wrap(lon))
}

/// Distance (km) from (x, y) to the segment from the origin along `bearing` for `len` km.
fn segment_distance(x: f64, y: f64, bearing: f64, len: f64) -> f64 {
    let (s, c) = bearing.to_radians().sin_cos();
    let along = (x * s + y * c).clamp(0.0, len);
    (x - along * s).hypot(y - along * c)
}

/// The point `dist` km from (lat, lon) along `bearing` (degrees clockwise from north).
fn destination(lat: f64, lon: f64, bearing: f64, dist: f64, radius_km: f64) -> (f64, f64) {
    let (p1, l1, d, b) = (lat.to_radians(), lon.to_radians(), dist / radius_km, bearing.to_radians());
    let p2 = (p1.sin() * d.cos() + p1.cos() * d.sin() * b.cos()).clamp(-1.0, 1.0).asin();
    let l2 = l1 + (b.sin() * d.sin() * p1.cos()).atan2(d.cos() - p1.sin() * p2.sin());
    (p2.to_degrees(), wrap(l2.to_degrees()))
}

fn great_circle_km(lat1: f64, lon1: f64, lat2: f64, lon2: f64, radius_km: f64) -> f64 {
    UnitVec3::from_lat_lon(lat1, lon1).angle_to(UnitVec3::from_lat_lon(lat2, lon2)) * radius_km
}

fn wrap(lon: f64) -> f64 {
    (lon + 180.0).rem_euclid(360.0) - 180.0
}

fn poisson(mean: f64, u: &mut impl FnMut() -> f64) -> usize {
    let (limit, mut p, mut k) = ((-mean).exp(), 1.0, 0);
    loop {
        p *= u();
        if p <= limit || k > 200 {
            return k;
        }
        k += 1;
    }
}

/// An index drawn with probability ∝ `weights`.
fn pick(weights: &[f64], u: f64) -> Option<usize> {
    let total: f64 = weights.iter().sum();
    if total <= 0.0 {
        return None;
    }
    let mut x = u * total;
    for (i, &w) in weights.iter().enumerate() {
        if x < w {
            return Some(i);
        }
        x -= w;
    }
    weights.iter().rposition(|&w| w > 0.0)
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    /// An Earth-like zonal climate: oceans everywhere but a band of land, westerlies at
    /// 35–65°, trades in the tropics.
    pub(crate) fn earthlike(seed: u64) -> Atmosphere {
        let season = Box::new(|t: f64| (std::f64::consts::TAU * t / 365.25).sin());
        Atmosphere::from_fn(6371.0, 24.0, 365.25, season, Seed(seed), |lat, lon| {
            let a = lat.abs();
            let wind = if a < 30.0 {
                [-6.0, -2.0 * lat.signum()]
            } else if a < 65.0 {
                [7.0, 0.5 * lat.signum()]
            } else {
                [-3.0, 0.0]
            };
            Local { mean_c: 30.0 - 0.0055 * a * a, range_c: 2.0 + a * 0.4, land: (0.0..40.0).contains(&lon) && a < 60.0, height_m: 0.0, wind_jan: wind, wind_jul: wind }
        })
    }

    /// The same with a 3000 m north–south range at 0–6°E from 35° to 65°N.
    fn with_range(seed: u64) -> Atmosphere {
        let base = earthlike(seed);
        let season = Box::new(|t: f64| (std::f64::consts::TAU * t / 365.25).sin());
        Atmosphere::from_fn(6371.0, 24.0, 365.25, season, Seed(seed), |lat, lon| {
            let i = base.cell(lat, lon);
            let high = (0.0..6.0).contains(&lon) && (35.0..65.0).contains(&lat);
            Local { mean_c: base.mean_c[i], range_c: base.range_c[i], land: base.land[i] || high, height_m: if high { 3000.0 } else { 0.0 }, wind_jan: base.wind_jan[i], wind_jul: base.wind_jul[i] }
        })
    }

    /// A range turns lows aside: far fewer cross it still deep, and more form in its lee.
    #[test]
    fn mountains_steer_storms() {
        let crossing = |atm: &Atmosphere| {
            let storms = atm.storms(0.0, 365.0);
            let crossed = storms.iter().filter(|s| s.kind == StormKind::Low && s.track.iter().any(|f| (1.0..5.0).contains(&f.lon) && (38.0..62.0).contains(&f.lat) && f.depth_hpa > 8.0)).count();
            let lee = storms.iter().filter(|s| s.kind == StormKind::Low && (6.0..12.0).contains(&s.track[0].lon) && (35.0..65.0).contains(&s.track[0].lat)).count();
            (crossed, lee)
        };
        let (flat, ridge) = (crossing(&earthlike(3)), crossing(&with_range(3)));
        assert!(ridge.0 * 3 < flat.0, "crossing deep: {} flat, {} with the range", flat.0, ridge.0);
        assert!(ridge.1 > flat.1, "born in the lee: {} flat, {} with the range", flat.1, ridge.1);
    }

    /// About 5 lows a day across both hemispheres and tens of tropical cyclones a year;
    /// lows drift east, tropical cyclones start west.
    #[test]
    fn storms_are_born_and_move() {
        let atm = earthlike(5);
        let storms = atm.storms(0.0, 365.0);
        let lows: Vec<_> = storms.iter().filter(|s| s.kind == StormKind::Low).collect();
        let tropical: Vec<_> = storms.iter().filter(|s| s.kind == StormKind::TropicalCyclone).collect();
        let per_day = lows.len() as f64 / 375.0;
        assert!((2.5..8.0).contains(&per_day), "{per_day} lows a day");
        assert!((30..200).contains(&tropical.len()), "{} tropical cyclones", tropical.len());
        let east = |s: &&Storm| wrap(s.track.last().unwrap().lon - s.track[0].lon);
        let mean_east = lows.iter().map(east).sum::<f64>() / lows.len() as f64;
        assert!(mean_east > 20.0, "lows drift {mean_east}° east");
        let first_day = |s: &&Storm| wrap(s.track[4].lon - s.track[0].lon);
        assert!(tropical.iter().map(first_day).sum::<f64>() < 0.0, "tropical cyclones start westward");
        assert!(lows.iter().all(|s| (23.5..71.5).contains(&s.track[0].lat.abs())));
    }

    /// Storms make neighbours' rain alike and distant places' rain unrelated; a deep low
    /// brings strong cyclonic winds.
    #[test]
    fn storms_are_shared_by_neighbours() {
        let atm = earthlike(8);
        let storms = atm.storms(0.0, 200.0);
        let series = |lat: f64, lon: f64| (0..800).map(|k| effect(&storms, lat, lon, f64::from(k) * 0.25, 6371.0, atm.omega).rain).collect::<Vec<_>>();
        let (a, near, far) = (series(50.0, -30.0), series(50.5, -29.0), series(50.0, 60.0));
        let near_r = correlation(&a, &near);
        let far_r = correlation(&a, &far);
        assert!(near_r > 0.7 && far_r < 0.3, "near {near_r}, far {far_r}");
        let fix = Fix { lat: 50.0, lon: 0.0, depth_hpa: 30.0, radius_km: 1000.0 };
        let deep = Storm { name: "Test".into(), kind: StormKind::Low, born_day: 0.0, track: vec![fix, fix] };
        // South of a northern low the wind blows from the west.
        let e = effect(std::slice::from_ref(&deep), 45.0, 0.0, 0.1, 6371.0, atm.omega);
        assert!(e.wind[0] > 8.0 && e.pressure_hpa < -15.0, "{e:?}");
    }

    pub(crate) fn correlation(a: &[f64], b: &[f64]) -> f64 {
        let n = a.len() as f64;
        let (ma, mb) = (a.iter().sum::<f64>() / n, b.iter().sum::<f64>() / n);
        let cov: f64 = a.iter().zip(b).map(|(x, y)| (x - ma) * (y - mb)).sum();
        let (va, vb) = (a.iter().map(|x| (x - ma).powi(2)).sum::<f64>(), b.iter().map(|y| (y - mb).powi(2)).sum::<f64>());
        cov / (va * vb).sqrt().max(1e-12)
    }
}
