//! Tectonics of an ice shell: one shell, old plains and the ice renewed from the ocean below,
//! and the surface's age.
//!
//! A thin shell (Europa's ~20 km) is young all over, ~60 Ma by its few craters (Bierhaus et
//! al. 2009; Leonard, Patthoff & Senske 2024), and about a third of it renewed more recently:
//! many small patches where the shell broke into chaos, and bands (assumed: a quarter chaos,
//! a twentieth bands; ridged plains cover the rest, the most widespread unit on the global
//! maps). A thick shell (Ganymede's) is old: dark cratered terrain ~4 Ga in large polygons,
//! and two thirds of it bright grooved terrain ~2 Ga (Patterson et al. 2010; Zahnle et al.
//! 2003). Renewed ice is stored as oceanic (thin) crust and old plains as continental.

use crate::noise::TerrainNoise;
use crate::preconditions;
use geoforge_core::grid::Cell;
use geoforge_core::key::Key;
use geoforge_core::product::{Generator, LayerId, Precondition, ProductData, ProductId, Products, RunContext};
use geoforge_core::shapes::{Event, EventKind, Field, Partition, Place, Table, Timeline};
use geoforge_core::surface::{CrustType, Plate};
use geoforge_core::{BodyContext, Error, UnitVec3};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// Every value `None` = chosen from the shell's thickness and the system's age.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct IceShellTectonicsParams {
    /// Share of the surface renewed: chaos and bands, or grooved terrain.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub renewed_fraction: Option<f64>,
    /// Age of the old plains, Ma.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub plains_age_ma: Option<f64>,
    /// Age of the renewed ice, Ma.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub renewed_age_ma: Option<f64>,
}

/// The shell's style with every default filled in.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShellStyle {
    /// 0 for a thin shell like Europa's, 1 for a thick one like Ganymede's.
    pub thick: f64,
    pub renewed_fraction: f64,
    pub plains_age_ma: f64,
    pub renewed_age_ma: f64,
}

/// Shells thinner than this are drawn as thin ones, thicker than the next as thick (km).
const THIN_KM: f64 = 25.0;
const THICK_KM: f64 = 60.0;

impl ShellStyle {
    #[must_use]
    pub fn of(body: &BodyContext, p: &IceShellTectonicsParams) -> Self {
        let km = body.interior.surface_shell().map_or(THICK_KM, |s| s.thickness_km);
        let t = ((km - THIN_KM) / (THICK_KM - THIN_KM)).clamp(0.0, 1.0);
        let thick = t * t * (3.0 - 2.0 * t);
        let oldest = (body.epoch.system_age_gyr * 1000.0 - 250.0).max(100.0);
        let log_mix = |thin: f64, thick_v: f64| (thin.ln() + thick * (thick_v.ln() - thin.ln())).exp();
        let plains = log_mix(60.0, oldest.min(4000.0));
        let renewed = log_mix(25.0, (0.5 * oldest).min(2000.0));
        Self {
            thick,
            renewed_fraction: p.renewed_fraction.unwrap_or(0.30 + 0.35 * thick).clamp(0.0, 1.0),
            plains_age_ma: p.plains_age_ma.unwrap_or(plains),
            renewed_age_ma: p.renewed_age_ma.unwrap_or(renewed).min(p.plains_age_ma.unwrap_or(plains)),
        }
    }
}

/// The value below which a share `q` of the total area lies.
fn area_quantile(values: &[f64], area: &[f64], q: f64) -> f64 {
    let mut idx: Vec<usize> = (0..values.len()).collect();
    idx.sort_by(|a, b| values[*a].total_cmp(&values[*b]));
    let target = q.clamp(0.0, 1.0) * area.iter().sum::<f64>();
    let mut acc = 0.0;
    for i in idx {
        acc += area[i];
        if acc >= target {
            return values[i];
        }
    }
    f64::INFINITY
}

pub struct IceShellTectonics;

impl Generator for IceShellTectonics {
    type Params = IceShellTectonicsParams;
    const NAME: &'static str = "ice_shell_tectonics";
    const VERSION: u32 = 1;
    const LAYER: LayerId = LayerId::Tectonics;

    fn needs() -> &'static [ProductId] {
        &[]
    }
    fn applies(body: &BodyContext) -> Result<(), Precondition> {
        preconditions::ice_shell(body)
    }

    fn run(ctx: &RunContext<'_>, inputs: &Products, p: &IceShellTectonicsParams) -> Result<Products, Error> {
        let grid = ctx.grid;
        let n = grid.len();
        let style = ShellStyle::of(ctx.body, p);
        let seed = ctx.seed(ProductId::CrustType);
        let pole = UnitVec3::from_lat_lon(0.0, 0.0);
        let mut out = Products::new();

        if inputs.plates().is_none() {
            let mut table = Table::new();
            table.insert(Key(1), Plate { name: "Shell".into(), euler_pole: pole, rate_deg_per_myr: 0.0, seed_point: pole });
            out.insert(ProductData::Plates(Partition { owner: Field::filled(grid.spec(), 1), table }));
        }

        let crust: Vec<CrustType> = match inputs.crust_type() {
            Some(c) => c.values.clone(),
            None => {
                // Renewed ice in patches: many small ones on a thin shell (chaos), large
                // polygons between old dark blocks on a thick one.
                let noise = TerrainNoise::new(seed.value());
                let (f1, f2) = (6.0 - 4.5 * style.thick, 16.0 - 10.0 * style.thick);
                let renewal: Vec<f64> = (0..n as Cell).into_par_iter().map(|c| noise.fbm(grid.center(c), f1, 4) + 0.35 * noise.fbm(grid.center(c), f2, 3)).collect();
                let area: Vec<f64> = grid.cells().map(|c| grid.area_sr(c)).collect();
                let cut = area_quantile(&renewal, &area, 1.0 - style.renewed_fraction);
                let crust: Vec<CrustType> = renewal.iter().map(|&v| if v >= cut { CrustType::Oceanic } else { CrustType::Continental }).collect();
                out.insert(ProductData::CrustType(Field::from_values(grid.spec(), crust.clone())));
                crust
            }
        };

        if inputs.crust_age().is_none() {
            let noise = TerrainNoise::new(ctx.seed(ProductId::CrustAge).value());
            let age: Vec<f32> = (0..n as Cell)
                .into_par_iter()
                .map(|c| {
                    let u = (0.5 + 0.5 * noise.fbm(grid.center(c), 3.0, 3)).clamp(0.0, 1.0);
                    let a = match crust[c as usize] {
                        CrustType::Continental => style.plains_age_ma * (0.8 + 0.2 * u),
                        CrustType::Oceanic => style.renewed_age_ma * (0.3 + 0.7 * u),
                    };
                    a as f32
                })
                .collect();
            out.insert(ProductData::CrustAge(Field::from_values(grid.spec(), age)));
        }

        if inputs.tectonic_events().is_none() {
            let mut events = Table::new();
            let note = if style.thick > 0.5 { "bright grooved ice floods between the old dark terrain" } else { "the shell breaks into chaos and opens into bands" };
            events.insert(
                Key(1),
                Event { t_ma: style.renewed_age_ma, kind: EventKind::ShellRenewal, at: Place::Point { at: pole }, magnitude: style.renewed_fraction, causes: Vec::new(), note: note.into() },
            );
            out.insert(ProductData::TectonicEvents(Timeline { events }));
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use geoforge_core::grid::{Grid, GridSpec};
    use geoforge_core::Archetype;

    #[test]
    fn europa_is_young_and_a_third_renewed() {
        let body = Archetype::Europa.context(Key(1), "Europa");
        let style = ShellStyle::of(&body, &IceShellTectonicsParams::default());
        assert!(style.thick < 0.05 && (50.0..70.0).contains(&style.plains_age_ma), "{style:?}");
        let grid = GridSpec::Equirect { width: 360, height: 180 }.build().unwrap();
        let ctx = RunContext { root_seed: geoforge_core::Seed(2), body: &body, grid: &grid };
        let out = IceShellTectonics::run(&ctx, &Products::new(), &IceShellTectonicsParams::default()).unwrap();
        let crust = &out.crust_type().unwrap().values;
        let renewed = grid.cells().filter(|&c| crust[c as usize] == CrustType::Oceanic).map(|c| grid.area_sr(c)).sum::<f64>() / (4.0 * std::f64::consts::PI);
        assert!((renewed - 0.30).abs() < 0.03, "renewed {renewed}");
        assert!(out.crust_age().unwrap().values.iter().all(|a| *a <= 60.0));
        // Titan's thick shell is drawn the Ganymede way: old, two thirds renewed long ago.
        let titan = ShellStyle::of(&Archetype::Titan.context(Key(1), "Titan"), &IceShellTectonicsParams::default());
        assert!(titan.thick > 0.95 && titan.plains_age_ma > 3000.0 && (0.6..0.7).contains(&titan.renewed_fraction), "{titan:?}");
    }
}
