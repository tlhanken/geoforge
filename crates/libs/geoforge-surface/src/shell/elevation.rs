//! Elevation of an ice shell: low relief, because ice flows. Europa's is a few hundred metres
//! over most of the surface: ridges 100–300 m high, chaos a little raised with blocks tilted
//! in it, bands a little above the plains, its largest craters shallow and flattened (Schenk
//! 2002; Nimmo & Schenk 2008). A thick shell holds more (Ganymede's grooves a few hundred
//! metres deep, regional relief of a kilometre or two), and its oldest basins are palimpsests,
//! nearly flat. Heights are the shell's own, not scaled from Earth's. Elevation is relative
//! to the mean radius.

use crate::lid::landforms::{basin, basin_depth_m, crater, crater_depth_m, craters, Canvas, BASIN_MIN_KM};
use crate::noise::TerrainNoise;
use crate::preconditions;
use crate::shell::tectonics::{IceShellTectonicsParams, ShellStyle};
use geoforge_core::grid::Cell;
use geoforge_core::product::{Generator, LayerId, Precondition, ProductData, ProductId, Products, RunContext};
use geoforge_core::shapes::Field;
use geoforge_core::surface::GeologicProvince as G;
use geoforge_core::{BodyContext, Error};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct IceShellElevationParams {
    /// Long-wavelength relief of a thin shell, m; a thick one has twice it.
    pub regional_m: f64,
    /// Roughness of the plains, m.
    pub roughness_m: f64,
    /// How far chaos stands above the plains, m, and the relief of its blocks.
    pub chaos_rise_m: f64,
    pub chaos_blocks_m: f64,
    /// How far a band stands above the plains, m.
    pub band_m: f64,
    /// Height of a ridge linea, m.
    pub ridge_m: f64,
    /// Relief of grooved terrain, m.
    pub groove_m: f64,
    /// Craters wider than this (km) on a thin shell are shallower than rock's, as ice relaxes.
    pub relaxation_km: f64,
    /// Smallest crater drawn, in grid cells across.
    pub crater_min_cells: f64,
}

impl Default for IceShellElevationParams {
    fn default() -> Self {
        Self { regional_m: 400.0, roughness_m: 120.0, chaos_rise_m: 100.0, chaos_blocks_m: 250.0, band_m: 100.0, ridge_m: 200.0, groove_m: 400.0, relaxation_km: 8.0, crater_min_cells: 3.0 }
    }
}

/// An old thick shell's basins keep this share of a rocky basin's depth.
const PALIMPSEST_DEPTH: f64 = 0.2;

pub struct IceShellElevation;

impl Generator for IceShellElevation {
    type Params = IceShellElevationParams;
    const NAME: &'static str = "ice_shell_elevation";
    const VERSION: u32 = 1;
    const LAYER: LayerId = LayerId::Elevation;

    fn needs() -> &'static [ProductId] {
        &[ProductId::CrustAge, ProductId::Provinces]
    }
    fn applies(body: &BodyContext) -> Result<(), Precondition> {
        preconditions::ice_shell(body)
    }

    fn run(ctx: &RunContext<'_>, inputs: &Products, p: &IceShellElevationParams) -> Result<Products, Error> {
        let mut out = Products::new();
        if inputs.rough_elevation().is_some() {
            return Ok(out);
        }
        let grid = ctx.grid;
        let n = grid.len();
        let r = ctx.radius_km();
        let style = ShellStyle::of(ctx.body, &IceShellTectonicsParams::default());
        let crust_age = &inputs.crust_age().ok_or(Error::MissingInput(ProductId::CrustAge))?.values;
        let provinces = inputs.provinces().ok_or(Error::MissingInput(ProductId::Provinces))?;
        let noise = TerrainNoise::new(ctx.seed(ProductId::Elevation).value());

        // 1. Each province's own ground, over the shell's long-wavelength relief.
        let regional = p.regional_m * (1.0 + style.thick);
        let ground: Vec<f64> = (0..n as Cell)
            .into_par_iter()
            .map(|c| {
                let pt = grid.center(c);
                let base = regional * noise.fbm(pt, 2.0, 3);
                let rough = p.roughness_m * noise.fbm(pt, 24.0, 3);
                let ridged = |amp: f64, f: f64| amp * 2.0 * (noise.ridged(pt, f, 3) - 0.5);
                let own = match provinces.row_at(c).map(|pv| (pv.kind, pv.intensity)) {
                    Some((G::ChaosTerrain, _)) => p.chaos_rise_m + ridged(p.chaos_blocks_m, 30.0),
                    Some((G::DilationalBand, k)) => p.band_m * k + 0.3 * rough,
                    Some((G::RidgeLinea, k)) => p.ridge_m * k + 0.5 * rough,
                    Some((G::GroovedTerrain, _)) => ridged(p.groove_m, 60.0),
                    Some((G::CrateredIce, _)) => 1.5 * rough,
                    _ => rough + ridged(0.5 * p.roughness_m, 40.0),
                };
                base + own
            })
            .collect();
        let mut canvas = Canvas { grid, radius_km: r, height_m: ground };

        // 2. Basins, flattened on old ice.
        for prov in provinces.table.rows.values().filter(|pv| pv.kind == G::ImpactBasin) {
            if let Some(at) = prov.anchor {
                canvas.stamp(at, prov.width_km / 2.0, 1.8, PALIMPSEST_DEPTH * basin_depth_m(prov.width_km), basin);
            }
        }

        // 3. Craters, as dense as the ground is old, shallow beyond the relaxation size.
        let surface_age: Vec<f32> =
            (0..n).map(|c| provinces.row_at(c as Cell).filter(|pv| pv.kind == G::RidgeLinea || pv.kind == G::DilationalBand).map_or(crust_age[c], |pv| pv.age_ma as f32)).collect();
        let d_min = p.crater_min_cells * grid.spacing_rad() * r;
        for imp in craters(grid, &surface_age, r, d_min, BASIN_MIN_KM, ctx.seed(ProductId::Elevation).child(b"craters").value()) {
            let relaxed = (p.relaxation_km / imp.diameter_km).min(1.0).powf(0.6);
            canvas.stamp(imp.at, imp.diameter_km / 2.0, 2.5, crater_depth_m(imp.diameter_km) * relaxed, |x| crater(x, imp.diameter_km > 30.0));
        }

        // 4. Datum: the mean radius.
        let e = canvas.height_m;
        let area: Vec<f64> = grid.cells().map(|c| grid.area_sr(c)).collect();
        let datum = e.iter().zip(&area).map(|(z, a)| z * a).sum::<f64>() / area.iter().sum::<f64>();
        out.insert(ProductData::RoughElevation(Field::from_values(grid.spec(), e.iter().map(|z| (z - datum) as f32).collect())));
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::shell::geology::tests::run;
    use geoforge_core::grid::{Grid, GridSpec};
    use geoforge_core::key::Key;
    use geoforge_core::{Archetype, Seed};

    /// Lowest, highest and the share of the surface within ±500 m.
    fn relief(a: Archetype) -> (f64, f64, f64) {
        let body = a.context(Key(1), a.slug());
        let grid = GridSpec::Equirect { width: 720, height: 360 }.build().unwrap();
        let ctx = RunContext { root_seed: Seed(5), body: &body, grid: &grid };
        let all = run(a);
        let e = IceShellElevation::run(&ctx, &all, &IceShellElevationParams::default()).unwrap().rough_elevation().unwrap().values.clone();
        let lo = f64::from(e.iter().copied().fold(f32::INFINITY, f32::min));
        let hi = f64::from(e.iter().copied().fold(f32::NEG_INFINITY, f32::max));
        let low = grid.cells().filter(|&c| e[c as usize].abs() < 500.0).map(|c| grid.area_sr(c)).sum::<f64>() / (4.0 * std::f64::consts::PI);
        (lo, hi, low)
    }

    #[test]
    fn ice_shells_have_low_relief() {
        // Europa: most of it within a few hundred metres, nothing as high as rock's mountains.
        let (lo, hi, within) = relief(Archetype::Europa);
        assert!(lo > -2500.0 && hi < 2000.0 && within > 0.75, "Europa {lo}..{hi}, {within} within ±500 m");
        // A thick, old shell holds more.
        let (tlo, thi, twithin) = relief(Archetype::Titan);
        assert!(thi - tlo > hi - lo && twithin < within && thi < 6000.0, "thick {tlo}..{thi}");
    }
}
