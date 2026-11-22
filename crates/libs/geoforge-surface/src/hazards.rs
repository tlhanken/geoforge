//! Natural hazards: how often each strikes a place, from the stored tectonics, climate,
//! circulation and relief (derived views, never stored), plus the planet's magnetic field.
//!
//! Each hazard has λ, its expected events per year in the 10,000 km² region (100 × 100 km)
//! around a cell; the chance of at least one in T years is P = 1 − e^(−λT). Rare hazards are
//! shown as that chance over a window that suits them (a year for storms, a decade for
//! quakes and floods, a century for tsunamis); frequent ones (fog, dust, blizzards, auroras)
//! count days a year, shown as the share of days.
//!
//! - **Earthquakes** (magnitude 6 and up): each plate boundary's seismic moment rate,
//!   μ·length·depth·slip·coupling (Kostrov 1974), spread over a truncated Gutenberg–Richter
//!   law (b = 1) up to a corner magnitude by boundary class (Bird & Kagan 2004: subduction
//!   9.5, collision 8.5, continental transform 8, rifts 7.5, ridges and oceanic transforms 7),
//!   and over a band either side of it. Intraplate quakes are rare everywhere.
//! - **Tsunamis**: coasts within 600 km of a subduction zone.
//! - **Volcanic eruptions**: arcs, the active ends of hotspot tracks, rifts, ridges above sea.
//! - **Tropical cyclones**: born over sea warmer than 26.5 °C in the warm season, 5–25° from
//!   the equator, at Earth's rate (~85 a year over ~6×10⁷ km² of such sea); they run west
//!   with the trades, curve poleward and east past 25°, and die over sea under 18 °C or a few
//!   hundred kilometres inland (Artifexian, Log 43).
//! - **Tornadoes**: flat warm plains 25–50° from the equator within reach of both high ground
//!   (dry air off it) and warm sea (moist air): Tornado Alley's recipe.
//! - **Lightning** in Artifexian's four tiers (Log 42), from warm months, big land, rain and
//!   windward slopes; rain shadows one tier down.
//! - **Floods, wildfire, landslides, blizzards, dust storms, fog** by climate and relief, and
//!   **auroras** by magnetic latitude.

use crate::ground::Ground;
use crate::tectonics::boundaries::{analyze, BoundaryParams};
use geoforge_core::geom::UnitVec3;
use geoforge_core::grid::{cells_within, Cell, Grid};
use geoforge_core::product::Products;
use geoforge_core::surface::{BoundaryKind, CrustType};

/// A natural hazard.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Hazard {
    Earthquakes,
    Tsunamis,
    Volcanoes,
    Cyclones,
    Tornadoes,
    Floods,
    Wildfire,
    Landslides,
    Blizzards,
    DustStorms,
    Fog,
    Auroras,
}

impl Hazard {
    pub const ALL: [Self; 12] = [
        Self::Earthquakes,
        Self::Tsunamis,
        Self::Volcanoes,
        Self::Cyclones,
        Self::Tornadoes,
        Self::Floods,
        Self::Wildfire,
        Self::Landslides,
        Self::Blizzards,
        Self::DustStorms,
        Self::Fog,
        Self::Auroras,
    ];

    #[must_use]
    pub fn slug(self) -> &'static str {
        self.info().0
    }

    #[must_use]
    pub fn from_slug(s: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|h| h.slug() == s)
    }

    /// What one event is.
    #[must_use]
    pub fn event(self) -> &'static str {
        self.info().1
    }

    /// When in the year it comes.
    #[must_use]
    pub fn season(self) -> Season {
        self.info().2
    }

    /// How it is shown: the chance over a window, or the share of days.
    #[must_use]
    pub fn measure(self) -> Measure {
        self.info().3
    }

    /// It can strike at sea too (quakes, cyclones, auroras); otherwise land only.
    #[must_use]
    pub fn at_sea(self) -> bool {
        matches!(self, Self::Earthquakes | Self::Cyclones | Self::Auroras)
    }

    fn info(self) -> (&'static str, &'static str, Season, Measure) {
        use Measure::{Days, Window};
        match self {
            Self::Earthquakes => ("earthquakes", "an earthquake of magnitude 6 or more", Season::Any, Window(10.0)),
            Self::Tsunamis => ("tsunamis", "a damaging tsunami", Season::Any, Window(100.0)),
            Self::Volcanoes => ("volcanoes", "an eruption", Season::Any, Window(10.0)),
            Self::Cyclones => ("cyclones", "a tropical cyclone passing", Season::Warm, Window(1.0)),
            Self::Tornadoes => ("tornadoes", "a tornado", Season::Warm, Window(1.0)),
            Self::Floods => ("floods", "a river flood", Season::Warm, Window(10.0)),
            Self::Wildfire => ("wildfire", "a large wildfire", Season::Warm, Window(10.0)),
            Self::Landslides => ("landslides", "a large landslide", Season::Any, Window(10.0)),
            Self::Blizzards => ("blizzards", "a blizzard day", Season::Cold, Days),
            Self::DustStorms => ("dust_storms", "a sand or dust storm day", Season::Warm, Days),
            Self::Fog => ("fog", "a day of thick fog", Season::Any, Days),
            Self::Auroras => ("auroras", "a night with auroras overhead", Season::Any, Days),
        }
    }
}

/// How a hazard's map shows it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Measure {
    /// The chance of at least one in this many years.
    Window(f64),
    /// The share of the year's days with one.
    Days,
}

impl Measure {
    /// The map value for a rate λ a year.
    #[must_use]
    pub fn value(self, lambda: f64) -> f64 {
        match self {
            Self::Window(years) => 1.0 - (-lambda * years).exp(),
            Self::Days => (lambda / 365.25).min(1.0),
        }
    }
}

/// When in the year a hazard comes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Season {
    Any,
    /// The local summer half-year.
    Warm,
    /// The local winter half-year.
    Cold,
}

/// What hazards read beyond the ground: tectonics, circulation and the magnetic field.
pub struct Setting<'a> {
    pub products: &'a Products,
    /// The magnetic pole and the auroral ovals' magnetic latitude, if the world has a field.
    pub aurora: Option<(UnitVec3, f64)>,
}

/// Rigidity of the crust, Pa.
const SHEAR_MODULUS: f64 = 3e10;
/// Moment of a magnitude, N·m (Hanks & Kanamori 1979).
fn moment(m: f64) -> f64 {
    10f64.powf(1.5 * m + 9.1)
}
/// Gutenberg–Richter β = 2b/3 for b = 1.
const BETA: f64 = 2.0 / 3.0;
/// Magnitude 6+ quakes per year per 10,000 km² far from any boundary (stable continents
/// have a few a decade over ~10⁸ km²).
const INTRAPLATE: f64 = 3e-5;

/// A boundary class: corner magnitude, seismic depth (km), coupling, and half-width of the
/// band its quakes fall in (km).
struct Class {
    corner: f64,
    depth_km: f64,
    coupling: f64,
    band_km: f64,
}

fn class(kind: BoundaryKind, continental: bool, oceanic: bool) -> Class {
    match (kind, continental, oceanic) {
        (BoundaryKind::Convergent, _, true) => Class { corner: 9.5, depth_km: 100.0, coupling: 1.0, band_km: 200.0 },
        (BoundaryKind::Convergent, _, false) => Class { corner: 8.5, depth_km: 30.0, coupling: 0.4, band_km: 300.0 },
        (BoundaryKind::Transform, true, _) => Class { corner: 8.0, depth_km: 15.0, coupling: 0.8, band_km: 80.0 },
        (BoundaryKind::Transform, false, _) => Class { corner: 7.0, depth_km: 10.0, coupling: 0.2, band_km: 40.0 },
        (BoundaryKind::Divergent, true, _) => Class { corner: 7.5, depth_km: 15.0, coupling: 0.3, band_km: 100.0 },
        (BoundaryKind::Divergent, false, _) => Class { corner: 7.0, depth_km: 8.0, coupling: 0.1, band_km: 40.0 },
    }
}

/// Earthquakes: magnitude 6+ a year per 10,000 km², and the largest magnitude likely, per
/// cell, from the plates and crust in `products`.
#[must_use]
pub fn earthquakes(grid: &dyn Grid, radius_km: f64, products: &Products) -> (Vec<f64>, Vec<f64>) {
    let n = grid.len();
    let mut rate = vec![INTRAPLATE; n];
    let mut corner = vec![6.5; n];
    let (Some(plates), Some(crust)) = (products.plates(), products.crust_type()) else { return (rate, corner) };
    let area_km2 = |c: Cell| grid.area_sr(c) * radius_km * radius_km;
    for seg in analyze(grid, plates, radius_km, &BoundaryParams::default()) {
        let continental = seg.cells.iter().any(|&c| crust.values[c as usize] == CrustType::Continental);
        let oceanic = seg.cells.iter().any(|&c| crust.values[c as usize] == CrustType::Oceanic);
        let k = class(seg.kind, continental, oceanic);
        // Seismic moment rate, N·m a year, and the magnitude 6+ rate a truncated
        // Gutenberg–Richter law gives it.
        let moment_rate = SHEAR_MODULUS * seg.length_km * 1e3 * k.depth_km * 1e3 * seg.rate_cm_yr.abs() * 1e-2 * k.coupling;
        let per_year = moment_rate * (1.0 - BETA) / (BETA * moment(6.0).powf(BETA) * moment(k.corner).powf(1.0 - BETA));
        let band = cells_within(grid, &seg.cells, k.band_km / radius_km, |_| true);
        let area: f64 = band.iter().map(|&c| area_km2(c)).sum();
        if area <= 0.0 {
            continue;
        }
        let density = per_year / area * 1e4;
        for c in band {
            rate[c as usize] += density;
            corner[c as usize] = corner[c as usize].max(k.corner);
        }
    }
    (rate, corner)
}

/// Chance per year of a quake of magnitude `m` or more, given the magnitude 6+ rate and the
/// corner magnitude: Gutenberg–Richter, b = 1, nothing beyond the corner.
#[must_use]
pub fn quake_rate_at_least(rate6: f64, corner: f64, m: f64) -> f64 {
    if m > corner {
        0.0
    } else {
        rate6 * 10f64.powf(-(m - 6.0))
    }
}

/// λ for `hazard` at every cell: events a year in the 10,000 km² around it.
#[must_use]
pub fn rate(g: &Ground, s: &Setting, hazard: Hazard) -> Vec<f64> {
    let n = g.len();
    match hazard {
        Hazard::Earthquakes => earthquakes(g.grid, g.radius_km, s.products).0,
        Hazard::Tsunamis => tsunamis(g, s),
        Hazard::Volcanoes => volcanoes(g, s),
        Hazard::Cyclones => cyclones(g.grid, g.radius_km, &g.sea, g.temperature, g.temperature_range),
        Hazard::Tornadoes => tornadoes(g),
        Hazard::Auroras => auroras(g, s),
        Hazard::Floods => floods(g, s),
        _ => (0..n).map(|c| if g.sea[c] { 0.0 } else { local(g, s, hazard, c) }).collect(),
    }
}

/// Floodplains: low, flat ground within 30 km of a big river (Strahler 4 and up). Monsoon
/// rivers flood most years, snowmelt rivers often, others now and then.
fn floods(g: &Ground, s: &Setting) -> Vec<f64> {
    let rivers: Vec<Cell> = (0..g.len()).filter(|&c| g.river_order[c] >= 4).map(|c| c as Cell).collect();
    let near = g.distance_over_land_km(rivers, 30.0);
    let share = s.products.precipitation_warm_share();
    (0..g.len())
        .map(|c| {
            if g.sea[c] || !near[c].is_finite() || g.slope[c] > 2.0 {
                return 0.0;
            }
            if share.is_some_and(|w| w.values[c] > 0.7) {
                0.5
            } else if g.group(c) == 'D' {
                0.3
            } else {
                0.1
            }
        })
        .collect()
}

/// Hazards that depend on the cell's own climate and relief.
fn local(g: &Ground, s: &Setting, hazard: Hazard, c: usize) -> f64 {
    let p = f64::from(g.precipitation[c]);
    let code = g.code(c);
    match hazard {
        Hazard::Wildfire => match code {
            "Csa" | "Csb" | "Aw" => 0.3,
            "BSh" | "BSk" | "Cwa" | "Cwb" => 0.15,
            "Dfc" | "Dfd" | "Dwc" | "Dwd" if g.biome[c] == "boreal_forest" => 0.1,
            _ if g.group(c) == 'B' || g.group(c) == 'E' => 0.0,
            _ => 0.02,
        },
        Hazard::Landslides if g.slope[c] > 30.0 => 0.2 * (g.slope[c] / 30.0).min(3.0) * if p > 1000.0 { 2.0 } else { 1.0 },
        Hazard::Blizzards if matches!(g.group(c), 'D' | 'E') => {
            let wind = winter_wind(s.products, c, g.lat[c]);
            if wind > 5.0 {
                2.0 + 2.0 * (wind - 5.0)
            } else {
                0.0
            }
        }
        Hazard::DustStorms if g.slope[c] < 3.0 => match g.group(c) {
            'B' if code.starts_with("BW") => 15.0,
            'B' => 5.0,
            'C' | 'D' if p < 500.0 => 2.0,
            _ => 0.0,
        },
        Hazard::Fog => {
            let cold_water = s.products.sst_anomaly().is_some_and(|a| a.values[c] < -1.5);
            if g.coast(c) && cold_water {
                60.0
            } else if g.coast(c) && g.group(c) == 'C' && code.ends_with('b') {
                10.0
            } else {
                0.0
            }
        }
        _ => 0.0,
    }
}

/// Winter wind speed at a cell, m/s.
fn winter_wind(p: &Products, c: usize, lat: f64) -> f64 {
    let w = if lat >= 0.0 { p.wind_jan() } else { p.wind_jul() };
    w.map_or(0.0, |w| f64::from(w.values[c][0]).hypot(f64::from(w.values[c][1])))
}

fn tsunamis(g: &Ground, s: &Setting) -> Vec<f64> {
    let n = g.len();
    let (Some(plates), Some(crust)) = (s.products.plates(), s.products.crust_type()) else { return vec![0.0; n] };
    let trench: Vec<Cell> = analyze(g.grid, plates, g.radius_km, &BoundaryParams::default())
        .into_iter()
        .filter(|seg| seg.kind == BoundaryKind::Convergent && seg.cells.iter().any(|&c| crust.values[c as usize] == CrustType::Oceanic))
        .flat_map(|seg| seg.cells)
        .collect();
    let near = g.distance_over_land_or_sea_km(trench, 600.0);
    (0..n).map(|c| if g.coast(c) && g.elevation[c] < 50.0 && near[c].is_finite() { 0.02 * (1.0 - near[c] / 600.0) } else { 0.0 }).collect()
}

fn volcanoes(g: &Ground, s: &Setting) -> Vec<f64> {
    use geoforge_core::surface::GeologicProvince as P;
    let n = g.len();
    let mut rate: Vec<f64> = (0..n)
        .map(|c| match g.province[c] {
            Some(P::VolcanicArc) => 0.1,
            Some(P::ContinentalRift) => 0.03,
            Some(P::MidOceanRidge) if !g.sea[c] => 0.1,
            _ => 0.0,
        })
        .collect();
    // Hotspots erupt at their active end.
    if let Some(prov) = s.products.provinces() {
        let ends: Vec<Cell> = prov.table.rows.values().filter(|p| matches!(p.kind, P::HotspotTrack | P::OceanicHotspotTrack | P::ContinentalHotspotTrack)).filter_map(|p| p.anchor).map(|a| g.grid.locate(a)).collect();
        let near = g.distance_over_land_or_sea_km(ends, 200.0);
        for c in 0..n {
            if near[c].is_finite() {
                rate[c] = rate[c].max(0.1);
            }
        }
    }
    rate
}

/// Storms per year per km² of warm sea where they are born (Earth: ~85 a year over ~6×10⁷ km²).
const CYCLONE_GENESIS: f64 = 85.0 / 6e7;
/// Coarse grid for tracks, degrees.
const TRACK_STEP: f64 = 2.0;

/// Tropical cyclones passing a year per 10,000 km², from the sea's annual mean temperature and
/// range (its warm season is the mean plus half the range).
#[must_use]
pub fn cyclones(grid: &dyn Grid, radius_km: f64, sea: &[bool], temperature: &[f32], range: &[f32]) -> Vec<f64> {
    let (w, h) = ((360.0 / TRACK_STEP) as usize, (180.0 / TRACK_STEP) as usize);
    let at = |i: usize, j: usize| {
        let (lat, lon) = (90.0 - (j as f64 + 0.5) * TRACK_STEP, -180.0 + (i as f64 + 0.5) * TRACK_STEP);
        grid.locate(UnitVec3::from_lat_lon(lat, lon)) as usize
    };
    let warm_sea = |c: usize| f64::from(temperature[c]) + 0.5 * f64::from(range[c]);
    let mut passes = vec![0.0; w * h];
    for j in 0..h {
        let lat0 = 90.0 - (j as f64 + 0.5) * TRACK_STEP;
        if !(5.0..=25.0).contains(&lat0.abs()) {
            continue;
        }
        for i in 0..w {
            let c = at(i, j);
            if !sea[c] || warm_sea(c) < 26.5 {
                continue;
            }
            let cell_km2 = (TRACK_STEP.to_radians() * radius_km).powi(2) * lat0.to_radians().cos();
            let storms = CYCLONE_GENESIS * cell_km2;
            // West with the trades, drifting poleward; past 25°, poleward and east.
            let (mut lon, mut lat, mut inland) = (-180.0 + (i as f64 + 0.5) * TRACK_STEP, lat0, 0);
            for _ in 0..60 {
                let (fi, fj) = (((lon + 180.0).rem_euclid(360.0) / TRACK_STEP) as usize % w, ((90.0 - lat) / TRACK_STEP).clamp(0.0, h as f64 - 1.0) as usize);
                passes[fj * w + fi] += storms;
                let c = at(fi, fj);
                if sea[c] {
                    inland = 0;
                    if warm_sea(c) < 18.0 {
                        break;
                    }
                } else {
                    inland += 1;
                    if inland > 2 {
                        break;
                    }
                }
                let poleward = lat.signum();
                if lat.abs() < 25.0 {
                    lon -= TRACK_STEP;
                    lat += poleward * 0.5;
                } else {
                    lon += TRACK_STEP;
                    lat += poleward * 1.5;
                }
                if lat.abs() > 60.0 {
                    break;
                }
            }
        }
    }
    // Spread over a storm's damaging swath (~300 km) and the region (100 km): a 3×3 blur of
    // the coarse cells, read back bilinearly.
    let coarse = |i: isize, j: isize| passes[(j.clamp(0, h as isize - 1) as usize) * w + i.rem_euclid(w as isize) as usize];
    let blurred: Vec<f64> = (0..w * h).map(|k| {
        let (i, j) = ((k % w) as isize, (k / w) as isize);
        (-1..=1).flat_map(|dj| (-1..=1).map(move |di| (di, dj))).map(|(di, dj)| coarse(i + di, j + dj)).sum::<f64>() / 9.0
    }).collect();
    let smooth = |x: f64, y: f64| {
        let (i0, j0) = (x.floor(), y.floor());
        let (fx, fy) = (x - i0, y - j0);
        let at = |i: f64, j: f64| blurred[((j as isize).clamp(0, h as isize - 1) as usize) * w + (i as isize).rem_euclid(w as isize) as usize];
        (1.0 - fy) * ((1.0 - fx) * at(i0, j0) + fx * at(i0 + 1.0, j0)) + fy * ((1.0 - fx) * at(i0, j0 + 1.0) + fx * at(i0 + 1.0, j0 + 1.0))
    };
    (0..grid.len())
        .map(|c| {
            let ll = grid.center(c as Cell).to_lat_lon();
            smooth((ll.lon_deg + 180.0) / TRACK_STEP - 0.5, (90.0 - ll.lat_deg) / TRACK_STEP - 0.5)
        })
        .collect()
}

fn tornadoes(g: &Ground) -> Vec<f64> {
    let n = g.len();
    let high: Vec<Cell> = (0..n).filter(|&c| g.elevation[c] > 1500.0).map(|c| c as Cell).collect();
    let warm_sea: Vec<Cell> = (0..n).filter(|&c| g.sea[c] && f64::from(g.temperature[c]) > 24.0).map(|c| c as Cell).collect();
    let near_high = g.distance_over_land_or_sea_km(high, 1500.0);
    let near_sea = g.distance_over_land_or_sea_km(warm_sea, 1500.0);
    (0..n)
        .map(|c| {
            let plain = !g.sea[c] && g.slope[c] < 2.0 && g.elevation[c] < 1000.0 && (25.0..=50.0).contains(&g.lat[c].abs());
            let wet = (400.0..1300.0).contains(&f64::from(g.precipitation[c])) && g.warmest(c) >= 22.0;
            if plain && wet && near_high[c].is_finite() && near_sea[c].is_finite() {
                4.0
            } else if plain && g.warmest(c) >= 18.0 && g.group(c) != 'B' {
                0.2
            } else {
                0.0
            }
        })
        .collect()
}

fn auroras(g: &Ground, s: &Setting) -> Vec<f64> {
    let Some((pole, oval)) = s.aurora else { return vec![0.0; g.len()] };
    (0..g.len())
        .map(|c| {
            let magnetic_lat = 90.0 - g.grid.center(c as Cell).angle_to(pole).to_degrees();
            let off = (magnetic_lat.abs() - oval).abs();
            // Nights a year with auroras overhead: most nights under the oval, some within
            // 10° of it, only in great storms 10–20° away.
            if off < 5.0 {
                150.0
            } else if magnetic_lat.abs() > oval - 15.0 {
                10.0
            } else if magnetic_lat.abs() > oval - 25.0 {
                0.5
            } else {
                0.0
            }
        })
        .collect()
}

/// Lightning tiers (Artifexian, Log 42): flashes per km² a year under 1, 1–5, 5–15, over 15.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Lightning {
    Rare,
    Basic,
    Moderate,
    Severe,
}

impl Lightning {
    pub const ALL: [Self; 4] = [Self::Rare, Self::Basic, Self::Moderate, Self::Severe];

    #[must_use]
    pub fn info(self) -> (&'static str, &'static str, [u8; 3]) {
        match self {
            Self::Rare => ("rare", "Under 1 flash per km² a year", [235, 235, 240]),
            Self::Basic => ("basic", "1–5 flashes per km² a year", [250, 220, 120]),
            Self::Moderate => ("moderate", "5–15 flashes per km² a year", [240, 140, 50]),
            Self::Severe => ("severe", "Over 15 flashes per km² a year", [180, 30, 60]),
        }
    }
}

/// The lightning tier at every cell, `None` at sea.
#[must_use]
pub fn lightning(g: &Ground) -> Vec<Option<Lightning>> {
    let n = g.len();
    let inland = g.distance_over_land_or_sea_km((0..n).filter(|&c| g.sea[c]).map(|c| c as Cell), 300.0);
    (0..n)
        .map(|c| {
            if g.sea[c] {
                return None;
            }
            let warm = g.warmest(c) >= 18.0 && !matches!(g.group(c), 'B' | 'E');
            if !warm {
                return Some(Lightning::Rare);
            }
            let p = f64::from(g.precipitation[c]);
            let zones = 1 + usize::from(!inland[c].is_finite()) + usize::from(p > 1000.0) + usize::from(g.slope[c] > 10.0 && p > 800.0);
            let mut tier = match zones {
                1 => 1,
                2 => 2,
                _ => 3,
            };
            if p < 500.0 {
                tier -= 1;
            }
            Some(Lightning::ALL[tier])
        })
        .collect()
}

