//! Climate layer: mean annual temperature, annual temperature range (warmest minus coldest
//! month) and annual precipitation.
//!
//! The model is deliberately physical in structure and cheap to run:
//!
//! 1. **Zonal energy balance** (Budyko–Sellers) on the body's own atmosphere: absorbed
//!    starlight by latitude against gray-atmosphere longwave (greenhouse from composition
//!    and pressure) and poleward heat transport, with an ice-albedo feedback on wet worlds.
//! 2. **Seasons** from the body's obliquity, eccentricity and year: the annual harmonic of
//!    insolation sets the seasonal swing, damped by the air column and the ocean.
//! 3. **Circulation** from rotation: Hadley, Ferrel and polar cells (the Hadley cell widens
//!    on slowly rotating bodies), shifted seasonally with the ITCZ, which moves further over
//!    large heated land; plus monsoon flow into the summer heat lows.
//! 4. **Continentality**: a maritime-air tracer carried by the winds, so west coasts in the
//!    westerlies are mild and continental interiors have large seasonal swings.
//! 5. **Moisture**: water vapour carried along the winds on a coarse atmosphere grid;
//!    evaporation over the sea relaxes it toward saturation, it rains out where air rises
//!    (convergence zones, windward slopes) or cools, and not where it sinks or warms (the
//!    Hadley cells' sinking belts, the eastern flanks of the subtropical highs, the poleward
//!    side of summer heat lows, air off cold water, rain shadows). Sinking air needs to be
//!    more humid before it rains, and summer land under it bakes. Land recycles part of its
//!    rain.
//!    The global total is set by a water budget scaled with temperature and ocean area.
//!    Worlds without surface water get none.
//! 6. **Downscaling** to the world grid: lapse rate on the cell's own elevation (less on
//!    broad plateaus, which heat the air above them), coast distance on the fine grid, local
//!    orographic enhancement.
//!
//! Radiation, grids, transport and the circulation cells live in [`crate::atmosphere`],
//! shared with the circulation layer.

mod moisture;

pub use moisture::w_sat;

use crate::atmosphere::transport::carry;
use crate::atmosphere::grid::{ground_m, is_sea};
use crate::atmosphere::{Atmo, AtmoGrid, Radiation, Season};
use geoforge_core::grid::{distance_field, Cell};
use geoforge_core::product::{Generator, LayerId, Precondition, ProductData, ProductId, Products, RunContext};
use geoforge_core::shapes::Field;
use geoforge_core::{BodyContext, Error};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct ClimateParams {
    /// Atmosphere grid spacing, degrees.
    pub resolution_deg: f64,
    /// `None` = from gravity: moist adiabat on wet worlds, near-dry elsewhere.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub lapse_rate_c_per_km: Option<f64>,
    /// Extra greenhouse warming (or cooling, negative) on top of the atmosphere's own, °C.
    pub greenhouse_offset_c: f64,
    /// Earth's poleward heat transport coefficient C, W/m²/K; scaled with the air column.
    pub heat_transport: f64,
    /// Global mean precipitation, mm/yr. `None` = from temperature and ocean area.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub global_precip_mm: Option<f64>,
    /// How far maritime air reaches inland, km.
    pub maritime_reach_km: f64,
    /// Rain-out rate per atmosphere step (sets how fast air dries inland).
    pub rainout: f64,
}

impl Default for ClimateParams {
    fn default() -> Self {
        Self {
            resolution_deg: 1.0,
            lapse_rate_c_per_km: None,
            greenhouse_offset_c: 0.0,
            heat_transport: 3.2,
            global_precip_mm: None,
            maritime_reach_km: 2500.0,
            rainout: 0.06,
        }
    }
}

impl ClimateParams {
    /// Cooling with height, °C/km: g/cₚ times 0.66 for a moist adiabat or 0.8 for dry air
    /// (Earth 6.5, Mars 3.0, Venus 7.1, Titan 1.1).
    #[must_use]
    pub fn lapse_rate(&self, body: &BodyContext) -> f64 {
        let moist = if body.interior.has_surface_liquid_water() { 0.66 } else { 0.8 };
        self.lapse_rate_c_per_km.unwrap_or(moist * body.gravity_m_s2())
    }
}

pub struct Climate;

/// Everything the climate layer computes on the atmosphere grid.
pub struct Coarse {
    pub grid: AtmoGrid,
    /// Sea-level-equivalent mean temperature, °C.
    pub t_sl: Vec<f64>,
    /// Seasonal swing coefficient: range = swing × (0.13 + 0.87 × continentality).
    pub swing: Vec<f64>,
    pub maritime: Vec<f64>,
    pub precip: Vec<f64>,
    /// Precipitation of the northern-summer (July) and northern-winter (January) seasons.
    pub precip_jul: Vec<f64>,
    pub precip_jan: Vec<f64>,
    pub elev: Vec<f64>,
    /// Ground height smoothed over ~300 km (plateaus), m.
    pub broad_elev: Vec<f64>,
    pub global_mean_t: f64,
}

/// Broad high ground heats the air above it: over a wide plateau the air is warmer than free
/// air at the same height (~2 °C per km of plateau on Tibet and the Altiplano). Not over ice
/// or snow, which absorb too little sunlight.
const PLATEAU_WARMING_C_PER_KM: f64 = 2.0;

fn plateau_warming(t_sea_level: f64, broad_elevation_m: f64) -> f64 {
    PLATEAU_WARMING_C_PER_KM * ((t_sea_level + 5.0) / 15.0).clamp(0.0, 1.0) * broad_elevation_m.max(0.0) / 1000.0
}

/// How far inland air chilled over cold water keeps its chill, km.
const COLD_AIR_REACH_KM: f64 = 300.0;
/// Extra summer warmth of land under fully sinking air (the moisture model's view of a
/// desert summer), °C.
const DESERT_SUMMER_C: f64 = 8.0;

fn continentality(is_land: bool, maritime: f64, coast_km: f64) -> f64 {
    if is_land {
        let coastal = 0.55 * (-coast_km / 250.0).exp();
        (1.0 - maritime.max(coastal)).clamp(0.0, 1.0)
    } else {
        (0.5 * (1.0 - maritime)).clamp(0.0, 0.5)
    }
}

impl Climate {
    /// Run the atmosphere on the coarse grid. Winds and the sea-surface temperature anomaly
    /// come from the circulation layer when present.
    fn coarse(ctx: &RunContext<'_>, inputs: &Products, elev: &[f32], p: &ClimateParams) -> Coarse {
        let grid = ctx.grid;
        let body = ctx.body;
        let r = ctx.radius_km();
        let wet = body.interior.has_surface_liquid_water();
        let lapse = p.lapse_rate(body);
        let atmo = Atmo::from_world(grid, elev, p.resolution_deg, wet);
        let g = atmo.g;
        let obl_w = body.irradiance.obliquity_deg;
        // The rain belt and monsoons follow the circulation layer's defaults.
        let circulation = crate::circulation::CirculationParams::default();
        let season_winds: Vec<Vec<(f64, f64)>> = match (inputs.wind_jul(), inputs.wind_jan()) {
            (Some(jul), Some(jan)) => vec![
                atmo.vectors_from_world(grid, &jul.values, crate::circulation::WIND_REF_M_S),
                atmo.vectors_from_world(grid, &jan.values, crate::circulation::WIND_REF_M_S),
            ],
            _ => [1.0, -1.0]
                .into_iter()
                .map(|sign| crate::circulation::season_winds(&atmo, &Season::new(&atmo, body.exterior.rotation_h, obl_w, sign, circulation.monsoon), &circulation))
                .collect(),
        };
        let sst: Vec<f64> = match inputs.sst_anomaly() {
            Some(a) => atmo.sea_mean_from_world(grid, &a.values, elev),
            None => vec![0.0; g.len()],
        };

        // 1. Zonal energy balance.
        let lats: Vec<f64> = (0..g.h).map(|y| g.lat(y)).collect();
        let obl = body.irradiance.obliquity_deg;
        let radiation = Radiation::for_body(body, p.heat_transport, p.greenhouse_offset_c);
        let (tz, tmean) = radiation.zonal_temperature(&lats);
        let swing_z: Vec<f64> = lats.iter().map(|&l| radiation.seasonal_range(l)).collect();

        // 2–4. Continentality and current-warmed air from annual-mean winds.
        let wind_annual: Vec<(f64, f64)> = season_winds[0].iter().zip(&season_winds[1]).map(|(a, b)| (0.5 * (a.0 + b.0), 0.5 * (a.1 + b.1))).collect();
        let mar = carry(&atmo, &wind_annual, p.maritime_reach_km, r, &vec![1.0; g.len()]);
        // Air over warm currents is warmed deep and carries it far inland (Europe); air over
        // cold upwelling is chilled in a shallow marine layer that the first hills stop.
        let warm: Vec<f64> = sst.iter().map(|v| v.max(0.0)).collect();
        let cold: Vec<f64> = sst.iter().map(|v| v.min(0.0)).collect();
        let sst_air: Vec<f64> = carry(&atmo, &wind_annual, p.maritime_reach_km * 0.6, r, &warm)
            .iter()
            .zip(carry(&atmo, &wind_annual, COLD_AIR_REACH_KM, r, &cold))
            .map(|(w, c)| w + c)
            .collect();
        let broad = g.blur(&atmo.elev, 3);
        let mut t_sl = vec![0.0; g.len()];
        let mut swing = vec![0.0; g.len()];
        let monsoon = circulation.monsoon;
        let seasons_def = [1.0, -1.0].map(|sign| Season::new(&atmo, body.exterior.rotation_h, obl, sign, monsoon));
        let mut t_season = [vec![0.0; g.len()], vec![0.0; g.len()]];
        let mut t_mean_here = vec![0.0; g.len()];
        for i in 0..g.len() {
            let y = i / g.w;
            let lat = lats[y];
            let is_land = atmo.land[i] >= 0.5;
            let cont = continentality(is_land, mar[i], if is_land { 400.0 } else { 0.0 });
            // Interiors at higher latitude are colder than the zonal mean; mild west coasts warmer.
            let hi = (lat.to_radians().sin().abs() - 0.3).max(0.0);
            let t = tz[y] - 7.0 * (cont - 0.35) * hi * atmo.land[i].max(0.3) + sst_air[i];
            t_sl[i] = t;
            swing[i] = swing_z[y];
            let range = swing[i].abs() * (0.13 + 0.87 * cont);
            let half = 0.5 * range * swing[i].signum();
            let t_here = t - lapse * atmo.elev[i] / 1000.0 + plateau_warming(t, broad[i].min(atmo.elev[i]));
            t_mean_here[i] = t_here;
            t_season[0][i] = t_here + half; // northern summer
            t_season[1][i] = t_here - half; // northern winter
            // Summer land under sinking air bakes: clear skies, dry ground, adiabatic warming.
            for (k, season) in seasons_def.iter().enumerate() {
                if lat * season.sign > 0.0 {
                    let sink = ((1.0 - season.convergence(i % g.w, lat, atmo.land[i])) / 0.6).clamp(0.0, 1.0);
                    t_season[k][i] += DESERT_SUMMER_C * atmo.land[i] * sink;
                }
            }
        }

        // 5. Moisture, two solstice seasons; nothing to rain without surface water.
        let mut seasons = Vec::new();
        for (k, season) in seasons_def.iter().enumerate() {
            seasons.push(if wet { moisture::moisture(&atmo, &season_winds[k], season, monsoon, &t_season[k], &t_mean_here, &sst_air, p.rainout, r) } else { vec![0.0; g.len()] });
        }
        let precip: Vec<f64> = seasons[0].iter().zip(&seasons[1]).map(|(a, b)| 0.5 * (a + b)).collect();
        let precip_jan = seasons.pop().expect("two seasons");
        let precip_jul = seasons.pop().expect("two seasons");
        Coarse { grid: g, t_sl, swing, maritime: mar, precip, precip_jul, precip_jan, elev: atmo.elev, broad_elev: broad, global_mean_t: tmean }
    }
}

impl Generator for Climate {
    type Params = ClimateParams;
    const NAME: &'static str = "climate";
    const VERSION: u32 = 6;
    const LAYER: LayerId = LayerId::Climate;

    fn needs() -> &'static [ProductId] {
        &[ProductId::Elevation]
    }
    fn wants() -> &'static [ProductId] {
        &[ProductId::WindJan, ProductId::WindJul, ProductId::SstAnomaly]
    }
    fn applies(body: &BodyContext) -> Result<(), Precondition> {
        crate::preconditions::surface_air(body)
    }

    fn run(ctx: &RunContext<'_>, inputs: &Products, p: &ClimateParams) -> Result<Products, Error> {
        let mut out = Products::new();
        if inputs.temperature().is_some() && inputs.temperature_range().is_some() && inputs.precipitation().is_some() && inputs.precipitation_warm_share().is_some() {
            return Ok(out);
        }
        let elev = &inputs.elevation().ok_or(Error::MissingInput(ProductId::Elevation))?.values;
        let grid = ctx.grid;
        let r = ctx.radius_km();
        let co = Self::coarse(ctx, inputs, elev, p);
        let g = co.grid;
        let wet = ctx.body.interior.has_surface_liquid_water();
        let lapse = p.lapse_rate(ctx.body);
        let sea = |i: usize| is_sea(wet, f64::from(elev[i]));

        // Fine coast distance (km), land side and sea side.
        let coast: Vec<Cell> = grid
            .cells()
            .filter(|&c| {
                let s = sea(c as usize);
                grid.neighbors(c).iter().any(|n| sea(n as usize) != s)
            })
            .collect();
        let coast_km: Vec<f64> = distance_field(grid, coast, 1500.0 / r, |_| true).into_iter().map(|d| (d * r).min(1500.0)).collect();
        let e_smooth = g.blur(&co.elev, 1);

        let n = grid.len();
        // Stationary waves and regional weather regimes: smooth anomalies so climate zones
        // meander instead of following parallels (±~1 °C, ±~25% precipitation).
        let wave_t = crate::noise::TerrainNoise::new(ctx.seed(ProductId::Temperature).value());
        let wave_p = crate::noise::TerrainNoise::new(ctx.seed(ProductId::Precipitation).value());
        let fine: Vec<(f32, f32, f32)> = (0..n)
            .into_par_iter()
            .map(|i| {
                let pt = grid.center(i as Cell);
                let ll = pt.to_lat_lon();
                let e = f64::from(elev[i]);
                let is_land = !sea(i);
                let mar = g.sample(&co.maritime, ll.lat_deg, ll.lon_deg);
                let cont = continentality(is_land, mar, coast_km[i]);
                let t_sl = g.sample(&co.t_sl, ll.lat_deg, ll.lon_deg) + 3.5 * wave_t.fbm(pt, 3.0, 3);
                // Plateau warming needs the ground itself to be high, not just nearby ranges.
                let broad = g.sample(&co.broad_elev, ll.lat_deg, ll.lon_deg).min(e);
                let t = t_sl - lapse * ground_m(wet, e) / 1000.0 + if is_land { plateau_warming(t_sl, broad) } else { 0.0 };
                let swing = g.sample(&co.swing, ll.lat_deg, ll.lon_deg).abs();
                let range = swing * (0.13 + 0.87 * cont);
                let mut pr = g.sample(&co.precip, ll.lat_deg, ll.lon_deg).max(0.0) * (0.8 * wave_p.fbm(pt, 4.0, 3)).exp();
                if is_land {
                    // Local relief above the atmosphere grid's smooth surface catches more rain.
                    let rel = (e - g.sample(&e_smooth, ll.lat_deg, ll.lon_deg)) / 1000.0;
                    pr *= (1.0 + 0.45 * rel).clamp(0.55, 1.9);
                }
                (t as f32, range as f32, pr as f32)
            })
            .collect();
        // Dry land heats more (little evaporative cooling, clear skies): warmer means and
        // larger swings in deserts. Uses the normalised precipitation below.
        let arid = |pr_mm: f64| (1.0 - pr_mm / 900.0).clamp(0.0, 1.0);

        // Water budget: global mean precipitation from temperature and ocean area.
        let area: Vec<f64> = (0..n).map(|i| grid.area_sr(i as Cell)).collect();
        let total_a: f64 = area.iter().sum();
        let ocean_frac = (0..n).filter(|&i| sea(i)).map(|i| area[i]).sum::<f64>() / total_a;
        let raw_mean = fine.iter().zip(&area).map(|(f, a)| f64::from(f.2) * a).sum::<f64>() / total_a;
        let target = p.global_precip_mm.unwrap_or_else(|| 1000.0 * 1.03f64.powf(co.global_mean_t - 14.5) * (ocean_frac / 0.71).max(0.0).powf(0.6));
        let k = if raw_mean > 0.0 { target / raw_mean } else { 0.0 };

        let spec = grid.spec();
        let dry_warm = |i: usize| -> (f32, f32) {
            if sea(i) {
                return (0.0, 0.0);
            }
            let a = arid(f64::from(fine[i].2) * k);
            let lat = grid.center(i as Cell).to_lat_lon().lat_deg.abs();
            let sub = (1.0 - ((lat - 25.0) / 30.0).powi(2)).clamp(0.0, 1.0);
            ((5.0 * a * sub) as f32, (7.0 * a * sub) as f32)
        };
        if inputs.temperature().is_none() {
            out.insert(ProductData::Temperature(Field::from_values(spec, (0..n).map(|i| fine[i].0 + dry_warm(i).0).collect())));
        }
        if inputs.temperature_range().is_none() {
            out.insert(ProductData::TemperatureRange(Field::from_values(spec, (0..n).map(|i| fine[i].1 + dry_warm(i).1).collect())));
        }
        if inputs.precipitation_warm_share().is_none() {
            let share: Vec<f32> = (0..n)
                .into_par_iter()
                .map(|i| {
                    let ll = grid.center(i as Cell).to_lat_lon();
                    let jul = g.sample(&co.precip_jul, ll.lat_deg, ll.lon_deg).max(0.0);
                    let jan = g.sample(&co.precip_jan, ll.lat_deg, ll.lon_deg).max(0.0);
                    let warm = if ll.lat_deg >= 0.0 { jul } else { jan };
                    if jul + jan > 0.0 { (warm / (jul + jan)) as f32 } else { 0.5 }
                })
                .collect();
            out.insert(ProductData::PrecipitationWarmShare(Field::from_values(spec, share)));
        }
        if inputs.precipitation().is_none() {
            out.insert(ProductData::Precipitation(Field::from_values(spec, fine.iter().map(|f| (f64::from(f.2) * k) as f32).collect())));
        }
        Ok(out)
    }
}

/// Potential evapotranspiration, mm/yr, from mean annual temperature (Thornthwaite-like,
/// calibrated to ~1650 mm at 25 °C, ~720 mm at 10 °C and ~190 mm at 0 °C).
#[must_use]
pub fn pet_mm(t: f64) -> f64 {
    if t <= -5.0 { 0.0 } else { 28.0 * (t + 5.0).powf(1.2) }
}
