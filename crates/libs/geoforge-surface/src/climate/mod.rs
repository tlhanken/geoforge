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
//! 3. **Circulation** from rotation and size: a Hadley cell (wider on slow rotators) and as
//!    many eddy-driven jets as the Rhines scale allows (Earth one: Ferrel and polar cells),
//!    shifted seasonally with the ITCZ, which moves further over large heated land; plus
//!    monsoon flow into the summer heat lows. A tidally locked world is balanced and
//!    circulated around its substellar point instead, without seasons.
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

use crate::atmosphere::grid::{ground_m, is_sea};
use crate::atmosphere::transport::carry;
use crate::atmosphere::{Atmo, AtmoGrid, Radiation, Season, Substellar};
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
            heat_transport: geoforge_core::energy_balance::EARTH_HEAT_TRANSPORT,
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
    /// How maritime the air is in the local summer and winter half-years, each carried
    /// inland by that season's winds.
    pub maritime_summer: Vec<f64>,
    pub maritime_winter: Vec<f64>,
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

/// Freezing point of seawater, °C.
const SEA_FREEZE_C: f64 = -1.8;
/// Warm anomaly at which a current keeps a sea open all winter, °C.
const OPEN_SEA_ANOMALY_C: f64 = 2.0;
/// How far inland air chilled over cold water keeps its chill, km.
const COLD_AIR_REACH_KM: f64 = 300.0;
/// Extra summer warmth of land under fully sinking air (the moisture model's view of a
/// desert summer), °C.
const DESERT_SUMMER_C: f64 = 8.0;

/// Extra winter cold of a fully continental high-latitude interior under its winter high: clear,
/// calm, long nights over snow let the ground radiate its heat away, and the cold air pools
/// in a shallow inversion (Yakutsk's January is −38 °C).
const WINTER_HIGH_C: f64 = 5.0;
/// Cooling of an interior below its zonal mean per unit continentality above a coast's (× the
/// latitude factor), °C.
const INTERIOR_C: f64 = 7.0;
/// How far sea air reaches inland in summer, as a share of the winter reach.
const SUMMER_REACH_SHARE: f64 = 0.75;
/// How much of open water's warmth a frozen sea still gives the winter air (through leads and
/// thin ice).
const FROZEN_SEA_WARMTH: f64 = 0.2;

/// A place's warmest- and coldest-month sea-level temperatures, from its annual zonal and
/// sea-warmed temperature `t`, the zonal seasonal swing, and how continental its air is in
/// each half-year. Summer and winter are set apart: the summer swing follows the summer
/// winds' continentality, and the winter one the winter winds' (Europe's winter westerlies
/// bring the Atlantic in), and a high-latitude interior has the winter high's cold besides
/// (none at sea, nor on a world without seasons).
fn summer_winter(t: f64, swing: f64, lat: f64, land: f64, cont_summer: f64, cont_winter: f64, seasons: bool) -> (f64, f64) {
    let swing = swing.abs();
    let (hi, polar) = if seasons { ((lat.to_radians().sin().abs() - 0.3).max(0.0), crate::elevation::smoothstep(45.0, 65.0, lat.abs())) } else { (0.0, 0.0) };
    let interior = INTERIOR_C * (0.5 * (cont_summer + cont_winter) - 0.35) * hi * land.max(0.3);
    let summer = t + 0.5 * swing * (0.13 + 0.87 * cont_summer) - interior;
    let high = WINTER_HIGH_C * land * ((cont_winter - 0.5) / 0.5).clamp(0.0, 1.0) * polar;
    let winter = t - 0.5 * swing * (0.13 + 0.87 * cont_winter) - interior - high;
    (summer, winter)
}

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
        // The rain belt and monsoons follow the circulation layer's defaults.
        let circulation = crate::circulation::CirculationParams::default();
        let season_winds: Vec<Vec<(f64, f64)>> = match (inputs.wind_jul(), inputs.wind_jan()) {
            (Some(jul), Some(jan)) => vec![atmo.vectors_from_world(grid, &jul.values, crate::circulation::WIND_REF_M_S), atmo.vectors_from_world(grid, &jan.values, crate::circulation::WIND_REF_M_S)],
            _ => [1.0, -1.0].into_iter().map(|sign| crate::circulation::season_winds(&atmo, &Season::for_body(&atmo, body, sign, circulation.monsoon), &circulation)).collect(),
        };
        // The highs' rising and sinking is read from the winds over flat ground: air turned
        // round a range has not risen or sunk for it.
        let flat_winds: Vec<Vec<(f64, f64)>> =
            [1.0, -1.0].into_iter().map(|sign| crate::circulation::season_winds_flat(&atmo, &Season::for_body(&atmo, body, sign, circulation.monsoon), &circulation)).collect();
        let sst: Vec<f64> = match inputs.sst_anomaly() {
            Some(a) => atmo.sea_mean_from_world(grid, &a.values, elev),
            None => vec![0.0; g.len()],
        };

        // 1. Zonal energy balance.
        let lats: Vec<f64> = (0..g.h).map(|y| g.lat(y)).collect();
        let radiation = Radiation::for_body(body, p.heat_transport, p.greenhouse_offset_c);
        // A tidally locked world: by angle from the substellar point, with no seasons.
        let locked = Substellar::for_body(body);
        let (tz, tmean) = if locked.is_some() { radiation.substellar_temperature(&lats) } else { radiation.zonal_temperature(&lats) };
        let swing_z: Vec<f64> = lats.iter().map(|&l| if locked.is_some() { 0.0 } else { radiation.seasonal_range(l) }).collect();

        // 2–4. Current-warmed air from annual-mean winds; continentality from each season's.
        let wind_annual: Vec<(f64, f64)> = season_winds[0].iter().zip(&season_winds[1]).map(|(a, b)| (0.5 * (a.0 + b.0), 0.5 * (a.1 + b.1))).collect();
        // Air over warm currents is warmed deep and carries it far inland (Europe); air over
        // cold upwelling is chilled in a shallow marine layer that the first hills stop.
        // Open water cannot be colder than freezing: where a warm current keeps a high-latitude
        // sea open, it stands above a zonal mean pulled down by frozen land and ice by that
        // much more (the Norwegian Sea).
        let sst: Vec<f64> = sst
            .iter()
            .enumerate()
            .map(|(i, &a)| {
                let t0 = tz[i / g.w];
                let open = (a / OPEN_SEA_ANOMALY_C).clamp(0.0, 1.0);
                a + (SEA_FREEZE_C - t0).max(0.0) * open * (1.0 - atmo.land[i])
            })
            .collect();
        // How far the sea's air reaches inland, by each season's winds, read back by hemisphere
        // as the local summer's and winter's. A sea frozen in a season (all year on the polar
        // seas) warms or cools that season's air little: under ice it is no longer maritime.
        let open_sea = |k: usize| -> Vec<f64> {
            (0..g.len())
                .map(|i| {
                    let lat = g.lat(i / g.w);
                    let winter = (k == 1) == (lat >= 0.0);
                    let half = 0.5 * swing_z[i / g.w].abs() * 0.13;
                    let t_sea = tz[i / g.w] + sst[i] + if winter { -half } else { half };
                    if t_sea < SEA_FREEZE_C {
                        FROZEN_SEA_WARMTH
                    } else {
                        1.0
                    }
                })
                .collect()
        };
        // (The sea's own air keeps its open-water reckoning: the sea-ice view is set on it.)
        // Summer sun heats sea air over land within days, so it reaches inland only part as
        // far as in winter.
        let maritime = |k: usize, reach: f64| -> Vec<f64> {
            let (iced, open) = (carry(&atmo, &season_winds[k], reach, r, &open_sea(k)), carry(&atmo, &season_winds[k], reach, r, &vec![1.0; g.len()]));
            (0..g.len()).map(|i| atmo.land[i] * iced[i] + (1.0 - atmo.land[i]) * open[i]).collect()
        };
        let summer_reach = SUMMER_REACH_SHARE * p.maritime_reach_km;
        let (jul_s, jan_s, jul_w, jan_w) = (maritime(0, summer_reach), maritime(1, summer_reach), maritime(0, p.maritime_reach_km), maritime(1, p.maritime_reach_km));
        let (mut mar_summer, mut mar_winter) = (vec![0.0; g.len()], vec![0.0; g.len()]);
        for i in 0..g.len() {
            let north = g.lat(i / g.w) >= 0.0;
            (mar_summer[i], mar_winter[i]) = if north { (jul_s[i], jan_w[i]) } else { (jan_s[i], jul_w[i]) };
        }
        let warm: Vec<f64> = sst.iter().map(|v| v.max(0.0)).collect();
        let cold: Vec<f64> = sst.iter().map(|v| v.min(0.0)).collect();
        let sst_air: Vec<f64> = carry(&atmo, &wind_annual, p.maritime_reach_km * 0.6, r, &warm).iter().zip(carry(&atmo, &wind_annual, COLD_AIR_REACH_KM, r, &cold)).map(|(w, c)| w + c).collect();
        let broad = g.blur(&atmo.elev, 3);
        let mut t_sl = vec![0.0; g.len()];
        let mut swing = vec![0.0; g.len()];
        let monsoon = circulation.monsoon;
        let seasons_def = [1.0, -1.0].map(|sign| Season::for_body(&atmo, body, sign, monsoon));
        let mut t_season = [vec![0.0; g.len()], vec![0.0; g.len()]];
        let mut t_mean_here = vec![0.0; g.len()];
        for i in 0..g.len() {
            let y = i / g.w;
            let lat = lats[y];
            let is_land = atmo.land[i] >= 0.5;
            let coast = if is_land { 400.0 } else { 0.0 };
            let (cont_s, cont_w) = (continentality(is_land, mar_summer[i], coast), continentality(is_land, mar_winter[i], coast));
            let t0 = match &locked {
                Some(s) => interpolate(&lats, &tz, 90.0 - s.zenith(lat, g.lon(i % g.w))),
                None => tz[y],
            };
            let t = t0 + sst_air[i];
            t_sl[i] = t;
            swing[i] = swing_z[y];
            let (summer, winter) = summer_winter(t, swing[i], lat, atmo.land[i], cont_s, cont_w, locked.is_none());
            let mean = 0.5 * (summer + winter);
            let half = 0.5 * (summer - winter) * swing[i].signum();
            let t_here = mean - lapse * atmo.elev[i] / 1000.0 + plateau_warming(mean, broad[i].min(atmo.elev[i]));
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
            seasons.push(if wet { moisture::moisture(&atmo, &season_winds[k], &flat_winds[k], season, monsoon, &t_season[k], &t_mean_here, &sst_air, p.rainout, r) } else { vec![0.0; g.len()] });
        }
        let precip: Vec<f64> = seasons[0].iter().zip(&seasons[1]).map(|(a, b)| 0.5 * (a + b)).collect();
        let precip_jan = seasons.pop().expect("two seasons");
        let precip_jul = seasons.pop().expect("two seasons");
        Coarse { grid: g, t_sl, swing, maritime_summer: mar_summer, maritime_winter: mar_winter, precip, precip_jul, precip_jan, elev: atmo.elev, broad_elev: broad, global_mean_t: tmean }
    }
}

impl Generator for Climate {
    type Params = ClimateParams;
    const NAME: &'static str = "climate";
    const VERSION: u32 = 11;
    const LAYER: LayerId = LayerId::Climate;

    fn needs() -> &'static [ProductId] {
        &[ProductId::RoughElevation]
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
        let elev = &inputs.rough_elevation().ok_or(Error::MissingInput(ProductId::RoughElevation))?.values;
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
                let cont_s = continentality(is_land, g.sample(&co.maritime_summer, ll.lat_deg, ll.lon_deg), coast_km[i]);
                let cont_w = continentality(is_land, g.sample(&co.maritime_winter, ll.lat_deg, ll.lon_deg), coast_km[i]);
                let t_zonal = g.sample(&co.t_sl, ll.lat_deg, ll.lon_deg) + 3.5 * wave_t.fbm(pt, 3.0, 3);
                let swing = g.sample(&co.swing, ll.lat_deg, ll.lon_deg);
                let (summer, winter) = summer_winter(t_zonal, swing, ll.lat_deg, if is_land { 1.0 } else { 0.0 }, cont_s, cont_w, swing != 0.0);
                let (t_sl, range) = (0.5 * (summer + winter), summer - winter);
                // Plateau warming needs the ground itself to be high, not just nearby ranges.
                let broad = g.sample(&co.broad_elev, ll.lat_deg, ll.lon_deg).min(e);
                let t = t_sl - lapse * ground_m(wet, e) / 1000.0 + if is_land { plateau_warming(t_sl, broad) } else { 0.0 };
                let mut pr = g.sample(&co.precip, ll.lat_deg, ll.lon_deg).max(0.0) * (0.8 * wave_p.fbm(pt, 4.0, 3)).exp();
                if is_land {
                    // Local relief above the atmosphere grid's smooth surface catches more rain,
                    // up to half as much again (a plateau's top behind its rim is not wetter
                    // than the slope that faces the wind).
                    let rel = (e - g.sample(&e_smooth, ll.lat_deg, ll.lon_deg)) / 1000.0;
                    pr *= (1.0 + 0.3 * rel).clamp(0.6, 1.5);
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
                    if jul + jan > 0.0 {
                        (warm / (jul + jan)) as f32
                    } else {
                        0.5
                    }
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

/// Temperature on the final heightmap. The climate ran on the rough heightmap, so ground the
/// erosion layer has cut lower since is warmer by the lapse rate (Earth's moist 6.5 °C/km;
/// a `lapse_rate_c_per_km` override in the climate's parameters is not seen here).
#[must_use]
pub fn on_ground(temperature: &[f32], rough: Option<&[f32]>, elevation: &[f32], body: &BodyContext) -> Vec<f32> {
    let Some(rough) = rough else { return temperature.to_vec() };
    let lapse = ClimateParams::default().lapse_rate(body) / 1000.0;
    temperature.iter().zip(rough.iter().zip(elevation)).map(|(&t, (&r, &e))| if e >= 0.0 && r >= 0.0 { t + (lapse * f64::from(r - e)) as f32 } else { t }).collect()
}

/// `f` at `x` by linear interpolation between samples `f[k]` at `xs[k]` (any order).
fn interpolate(xs: &[f64], f: &[f64], x: f64) -> f64 {
    let mut idx: Vec<usize> = (0..xs.len()).collect();
    idx.sort_by(|&a, &b| xs[a].total_cmp(&xs[b]));
    let (first, last) = (idx[0], idx[idx.len() - 1]);
    if x <= xs[first] {
        return f[first];
    }
    if x >= xs[last] {
        return f[last];
    }
    let k = idx.windows(2).find(|w| x <= xs[w[1]]).expect("x inside the samples");
    let (a, b) = (k[0], k[1]);
    f[a] + (f[b] - f[a]) * (x - xs[a]) / (xs[b] - xs[a])
}

/// Potential evapotranspiration, mm/yr, from mean annual temperature (Thornthwaite-like,
/// calibrated to ~1650 mm at 25 °C, ~720 mm at 10 °C and ~190 mm at 0 °C).
#[must_use]
pub fn pet_mm(t: f64) -> f64 {
    if t <= -5.0 {
        0.0
    } else {
        28.0 * (t + 5.0).powf(1.2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A continental interior at 62°N: a much colder winter than the sea's, and a summer only
    /// somewhat warmer (the winter high).
    #[test]
    fn a_high_latitude_interior_has_its_winter_high() {
        let (sea_s, sea_w) = summer_winter(-2.0, 46.0, 62.0, 0.0, 0.05, 0.05, true);
        let (land_s, land_w) = summer_winter(-2.0, 46.0, 62.0, 1.0, 0.9, 0.9, true);
        assert!(land_s > sea_s + 5.0 && sea_w - land_w > land_s - sea_s + 5.0, "sea {sea_s:.1}/{sea_w:.1}, land {land_s:.1}/{land_w:.1}");
        // No winter high at 30°N: the range is the swing's continental share alone.
        let (s, w) = summer_winter(15.0, 25.0, 30.0, 1.0, 0.9, 0.9, true);
        assert!((s - w - 25.0 * (0.13 + 0.87 * 0.9)).abs() < 1e-9, "{s} {w}");
    }

    /// Maritime winters and continental summers (Europe under its westerlies) make a warmer
    /// year than continental winters and maritime summers.
    #[test]
    fn maritime_winters_warm_the_year() {
        let (es, ew) = summer_winter(5.0, 43.0, 55.0, 1.0, 0.7, 0.4, true);
        let (ms, mw) = summer_winter(5.0, 43.0, 55.0, 1.0, 0.4, 0.7, true);
        assert!(es + ew > ms + mw + 3.0 && es - ew < ms - mw);
    }

    /// Without seasons nothing swings.
    #[test]
    fn a_locked_world_has_no_seasons() {
        assert_eq!(summer_winter(10.0, 0.0, 60.0, 1.0, 0.9, 0.9, false), (10.0, 10.0));
    }
}
