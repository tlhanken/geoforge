//! Elevation of a rigid lid: the crust's highland–lowland step, each province's landform,
//! craters by surface age, and roughness. Heights scale with gravity (see
//! [`relief_scale`]). Without surface water elevation is relative to the mean radius; with
//! it, to the sea level that holds the water.

use super::landforms::{basin, basin_depth_m, corona, crater, crater_depth_m, craters, dome, relief_scale, shield, Canvas, BASIN_MIN_KM};
use crate::elevation::sea_level_for_volume;
use crate::noise::TerrainNoise;
use crate::preconditions;
use geoforge_core::grid::{distance_field, Cell};
use geoforge_core::product::{Generator, LayerId, Precondition, ProductData, ProductId, Products, RunContext};
use geoforge_core::shapes::Field;
use geoforge_core::surface::{CrustType, GeologicProvince as G};
use geoforge_core::{BodyContext, Error};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// Heights of loads the lithosphere holds up (volcanoes, rises, rifts) are at Earth gravity
/// and scale with the body's; the highland–lowland step is isostatic and does not.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct LidElevationParams {
    /// Step from lowlands up to highlands, m. `None` = from the crust: thin lowland crust
    /// floats lower by about a ninth of the crust's thickness (Mars ≈ 5.5 km).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub dichotomy_m: Option<f64>,
    /// Width of the step, km.
    pub dichotomy_width_km: f64,
    /// Height of a full-strength volcanic rise, m.
    pub rise_m: f64,
    /// Height of a full-strength shield volcano, m (Mauna Loa from its base ≈ 9 km).
    pub shield_m: f64,
    pub rift_depth_m: f64,
    pub corona_m: f64,
    /// Smallest crater drawn, in grid cells across.
    pub crater_min_cells: f64,
    /// Amplitude of rough terrain, m.
    pub roughness_m: f64,
}

impl Default for LidElevationParams {
    fn default() -> Self {
        Self { dichotomy_m: None, dichotomy_width_km: 400.0, rise_m: 3000.0, shield_m: 8500.0, rift_depth_m: 2800.0, corona_m: 1000.0, crater_min_cells: 3.0, roughness_m: 400.0 }
    }
}

/// Highland–lowland step per metre of crust: a lowland crust ~60% thinner, floating on a
/// mantle ~12% denser than the crust, plus the highlands' thicker cover.
const ISOSTATIC_STEP: f64 = 0.11;

/// No volcano taller than this, m, however weak the gravity.
const MAX_SHIELD_M: f64 = 26_000.0;

fn smoothstep(x: f64) -> f64 {
    let t = x.clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

pub struct LidElevation;

impl Generator for LidElevation {
    type Params = LidElevationParams;
    const NAME: &'static str = "lid_elevation";
    const VERSION: u32 = 1;
    const LAYER: LayerId = LayerId::Elevation;

    fn needs() -> &'static [ProductId] {
        &[ProductId::CrustType, ProductId::CrustAge, ProductId::Provinces]
    }
    fn applies(body: &BodyContext) -> Result<(), Precondition> {
        preconditions::rigid_lid(body)
    }

    fn run(ctx: &RunContext<'_>, inputs: &Products, p: &LidElevationParams) -> Result<Products, Error> {
        let mut out = Products::new();
        if inputs.rough_elevation().is_some() {
            return Ok(out);
        }
        // Randomness is drawn from the `elevation` seed, as before the heightmap was split.
        let grid = ctx.grid;
        let n = grid.len();
        let r = ctx.radius_km();
        let scale = relief_scale(ctx.body);
        let crust = &inputs.crust_type().ok_or(Error::MissingInput(ProductId::CrustType))?.values;
        let crust_age = &inputs.crust_age().ok_or(Error::MissingInput(ProductId::CrustAge))?.values;
        let provinces = inputs.provinces().ok_or(Error::MissingInput(ProductId::Provinces))?;
        let high = |c: usize| crust[c] == CrustType::Continental;

        // 1. The highland–lowland step, smoothed across its width.
        let edge: Vec<Cell> = grid.cells().filter(|&c| grid.neighbors(c).iter().any(|nb| high(nb as usize) != high(c as usize))).collect();
        let w_rad = p.dichotomy_width_km / r;
        let d = distance_field(grid, edge, w_rad, |_| true);
        let crust_m = ctx.body.interior.surface_shell().map_or(30_000.0, |s| s.thickness_km * 1000.0);
        let step = p.dichotomy_m.unwrap_or(ISOSTATIC_STEP * crust_m);
        let base: Vec<f64> = (0..n)
            .map(|c| {
                let side = if high(c) { 1.0 } else { -1.0 };
                let t = 0.5 + 0.5 * side * (d[c] / w_rad).min(1.0);
                -0.6 * step + step * smoothstep(t)
            })
            .collect();
        let mut canvas = Canvas { grid, radius_km: r, height_m: base };

        // 2. Landforms, in the provinces' painting order (old under young).
        for (key, prov) in &provinces.table.rows {
            let (Some(at), radius) = (prov.anchor, prov.width_km / 2.0) else { continue };
            match prov.kind {
                G::ImpactBasin => canvas.stamp(at, radius, 1.8, basin_depth_m(prov.width_km) * scale.powf(0.2), basin),
                G::VolcanicRise => canvas.stamp(at, radius, 1.0, p.rise_m * scale * prov.intensity, dome),
                G::ShieldVolcano => canvas.stamp(at, radius, 1.0, (p.shield_m * scale * prov.intensity).min(MAX_SHIELD_M), shield),
                G::Corona => canvas.stamp(at, radius, 1.1, p.corona_m * scale.sqrt() * prov.intensity, corona),
                G::ContinentalRift => {
                    let own = key.0 as u32;
                    let outside: Vec<Cell> =
                        grid.cells().filter(|&c| provinces.owner.values[c as usize] != own && grid.neighbors(c).iter().any(|nb| provinces.owner.values[nb as usize] == own)).collect();
                    let half = (prov.width_km / 2.0) / r;
                    let inside = distance_field(grid, outside, half, |c| provinces.owner.values[c as usize] == own);
                    for (c, di) in inside.iter().enumerate() {
                        if provinces.owner.values[c] == own && di.is_finite() {
                            canvas.height_m[c] -= p.rift_depth_m * scale * prov.intensity * smoothstep(di / half);
                        }
                    }
                }
                _ => {}
            }
        }

        // 3. Craters, as dense as the surface is old: volcanic ground is as young as its province.
        let surface_age: Vec<f32> = (0..n)
            .map(|c| match provinces.row_at(c as Cell) {
                Some(pv) if matches!(pv.kind, G::VolcanicRise | G::ShieldVolcano | G::Corona | G::ImpactBasin) => pv.age_ma as f32,
                _ => crust_age[c],
            })
            .collect();
        let d_min = p.crater_min_cells * grid.spacing_rad() * r;
        for imp in craters(grid, &surface_age, r, d_min, BASIN_MIN_KM, ctx.seed(ProductId::Elevation).child(b"craters").value()) {
            let peak = imp.diameter_km > 30.0;
            canvas.stamp(imp.at, imp.diameter_km / 2.0, 2.5, crater_depth_m(imp.diameter_km) * scale.powf(0.3), |x| crater(x, peak));
        }
        let mut e = canvas.height_m;

        // 4. Roughness: rugged highlands and tesserae, smoother plains.
        let noise = TerrainNoise::new(ctx.seed(ProductId::Elevation).value());
        let rough: Vec<f64> = (0..n as Cell)
            .into_par_iter()
            .map(|c| {
                let pt = grid.center(c);
                let kind = provinces.row_at(c).map(|pv| pv.kind);
                let base = p.roughness_m * scale * noise.fbm(pt, 6.0, 4);
                match kind {
                    Some(G::Tessera) => base + 1.5 * p.roughness_m * scale * (noise.ridged(pt, 14.0, 3) - 0.5),
                    Some(G::VolcanicPlains) => 0.4 * base,
                    _ if high(c as usize) => base,
                    _ => 0.5 * base,
                }
            })
            .collect();
        for (z, dz) in e.iter_mut().zip(&rough) {
            *z += dz;
        }

        // 5. Datum: sea level holding the surface water, or the mean radius.
        let area: Vec<f64> = grid.cells().map(|c| grid.area_sr(c)).collect();
        let gel_m = ctx.body.interior.surface_ocean_gel_km() * 1000.0;
        let datum = if ctx.body.interior.has_surface_liquid_water() && gel_m > 0.0 {
            sea_level_for_volume(&e, &area, gel_m)
        } else {
            e.iter().zip(&area).map(|(z, a)| z * a).sum::<f64>() / area.iter().sum::<f64>()
        };
        out.insert(ProductData::RoughElevation(Field::from_values(grid.spec(), e.iter().map(|z| (z - datum) as f32).collect())));
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lid::{LidGeology, LidGeologyParams, LidTectonics, LidTectonicsParams};
    use geoforge_core::grid::{Grid, GridSpec};
    use geoforge_core::key::Key;
    use geoforge_core::{Archetype, Seed};

    fn relief(a: Archetype) -> (f64, f64, f64) {
        let body = a.context(Key(1), a.slug());
        let grid = GridSpec::Equirect { width: 512, height: 256 }.build().unwrap();
        let ctx = RunContext { root_seed: Seed(9), body: &body, grid: &grid };
        let mut all = LidTectonics::run(&ctx, &Products::new(), &LidTectonicsParams::default()).unwrap();
        all.extend(LidGeology::run(&ctx, &all, &LidGeologyParams::default()).unwrap());
        all.extend(LidElevation::run(&ctx, &all, &LidElevationParams::default()).unwrap());
        let e = &all.rough_elevation().unwrap().values;
        let lo = f64::from(e.iter().copied().fold(f32::INFINITY, f32::min));
        let hi = f64::from(e.iter().copied().fold(f32::NEG_INFINITY, f32::max));
        let mean = e.iter().zip(grid.cells()).map(|(z, c)| f64::from(*z) * grid.area_sr(c)).sum::<f64>() / (4.0 * std::f64::consts::PI);
        (lo, hi, mean)
    }

    #[test]
    fn relief_matches_real_bodies() {
        // Mars: Hellas floor ≈ −8 km, Olympus Mons ≈ +22 km (relative to the mean).
        let (lo, hi, mean) = relief(Archetype::Mars);
        assert!(mean.abs() < 1.0, "datum is the mean radius");
        assert!((-12_000.0..-5_000.0).contains(&lo), "Mars lowest {lo}");
        assert!((12_000.0..27_000.0).contains(&hi), "Mars highest {hi}");
        // Venus: lower, smoother relief; Maat Mons ≈ 8 km above the plains.
        let (lo, hi, _) = relief(Archetype::Venus);
        assert!(hi - lo < 16_000.0 && hi > 4000.0, "Venus {lo}..{hi}");
        // The Moon: deep basins, no great volcanoes.
        let (lo, hi, _) = relief(Archetype::Moon);
        assert!(lo < -6000.0 && hi < 12_000.0, "Moon {lo}..{hi}");
    }
}
