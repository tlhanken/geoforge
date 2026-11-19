//! Tectonics of a rigid lid: one plate, a crust split into old highlands and younger
//! lowlands, and the ground's age.
//!
//! Stagnant lids (Mars, Mercury, the Moon) gather their lowlands in one hemisphere (Mars's
//! northern plains, the Moon's nearside maria); an episodic lid (Venus) renewed most of its
//! surface at once, leaving scattered older highlands. Lowlands are stored as oceanic
//! (thin) crust and highlands as continental (thick) crust.

use crate::noise::TerrainNoise;
use crate::preconditions;
use geoforge_core::body::{TectonicRegime, EARTH_MASS_KG};
use geoforge_core::grid::Cell;
use geoforge_core::key::Key;
use geoforge_core::product::{Generator, LayerId, Precondition, ProductData, ProductId, Products, RunContext};
use geoforge_core::shapes::{Event, EventKind, Field, Partition, Place, Table, Timeline};
use geoforge_core::surface::{CrustType, Plate};
use geoforge_core::{BodyContext, Error, UnitVec3};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// Every value `None` = chosen from the body's regime and age.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct LidTectonicsParams {
    /// Share of the surface that is lowland.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub lowland_fraction: Option<f64>,
    /// How strongly the lowlands gather in one hemisphere: 1 = a single cap (Mars),
    /// 0 = scattered (Venus).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub dichotomy: Option<f64>,
    /// Age of the highland surface, Ma.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub highland_age_ma: Option<f64>,
    /// Age of the lowland surface, Ma.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub lowland_age_ma: Option<f64>,
}

/// The parameters with every default filled in for a body.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LidStyle {
    pub lowland_fraction: f64,
    pub dichotomy: f64,
    pub highland_age_ma: f64,
    pub lowland_age_ma: f64,
}

impl LidTectonicsParams {
    #[must_use]
    pub fn style(&self, body: &BodyContext) -> LidStyle {
        let oldest = body.epoch.system_age_gyr * 1000.0 - 250.0;
        // Small bodies cool before flooding much: lava plains cover ~17% of the Moon,
        // ~27% of Mercury, ~35% of Mars.
        let x = (body.exterior.mass_kg / EARTH_MASS_KG / 0.1).clamp(0.0, 1.0);
        let stagnant_lowland = 0.17 + 0.18 * x * x * (3.0 - 2.0 * x);
        let (lowland, dichotomy, high, low) = match body.interior.regime {
            // Venus: tesserae up to ~1 Ga, plains renewed ~500 Ma ago.
            TectonicRegime::Episodic => (0.8, 0.1, 900.0, 450.0),
            // Dead since it formed: no lava plains, all old ground.
            TectonicRegime::Inactive => (0.0, 0.0, oldest + 100.0, oldest + 100.0),
            _ => (stagnant_lowland, 0.75, oldest, 0.82 * oldest),
        };
        LidStyle {
            lowland_fraction: self.lowland_fraction.unwrap_or(lowland).clamp(0.0, 1.0),
            dichotomy: self.dichotomy.unwrap_or(dichotomy).clamp(0.0, 1.0),
            highland_age_ma: self.highland_age_ma.unwrap_or(high),
            lowland_age_ma: self.lowland_age_ma.unwrap_or(low),
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

pub struct LidTectonics;

impl Generator for LidTectonics {
    type Params = LidTectonicsParams;
    const NAME: &'static str = "lid_tectonics";
    const VERSION: u32 = 1;
    const LAYER: LayerId = LayerId::Tectonics;

    fn needs() -> &'static [ProductId] {
        &[]
    }
    fn applies(body: &BodyContext) -> Result<(), Precondition> {
        preconditions::rigid_lid(body)
    }

    fn run(ctx: &RunContext<'_>, inputs: &Products, p: &LidTectonicsParams) -> Result<Products, Error> {
        let grid = ctx.grid;
        let n = grid.len();
        let style = p.style(ctx.body);
        let seed = ctx.seed(ProductId::CrustType);
        let pole = UnitVec3::from_lat_lon(geoforge_core::seed::range(seed.child(b"lat"), -60.0, 60.0), geoforge_core::seed::range(seed.child(b"lon"), -180.0, 180.0));
        let mut out = Products::new();

        if inputs.plates().is_none() {
            let mut table = Table::new();
            table.insert(Key(1), Plate { name: "Lid".into(), euler_pole: pole, rate_deg_per_myr: 0.0, seed_point: pole });
            out.insert(ProductData::Plates(Partition { owner: Field::filled(grid.spec(), 1), table }));
        }

        let crust: Vec<CrustType> = match inputs.crust_type() {
            Some(c) => c.values.clone(),
            None => {
                // Lowness: the dichotomy cap plus warped noise, cut at the lowland share.
                let noise = TerrainNoise::new(seed.value());
                let low: Vec<f64> = (0..n as Cell)
                    .into_par_iter()
                    .map(|c| {
                        let pt = grid.center(c);
                        style.dichotomy * pt.dot(pole) + (1.0 - style.dichotomy) * 1.5 * noise.fbm(pt, 1.5, 4) + 0.25 * noise.fbm(pt, 4.0, 3)
                    })
                    .collect();
                let area: Vec<f64> = grid.cells().map(|c| grid.area_sr(c)).collect();
                let cut = area_quantile(&low, &area, 1.0 - style.lowland_fraction);
                let crust: Vec<CrustType> = low.iter().map(|&v| if v >= cut { CrustType::Oceanic } else { CrustType::Continental }).collect();
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
                        CrustType::Continental => style.highland_age_ma * (0.93 + 0.07 * u),
                        CrustType::Oceanic => (style.lowland_age_ma * (0.85 + 0.25 * u)).min(style.highland_age_ma),
                    };
                    a as f32
                })
                .collect();
            out.insert(ProductData::CrustAge(Field::from_values(grid.spec(), age)));
        }

        if inputs.tectonic_events().is_none() {
            let (kind, t_ma, note) = match ctx.body.interior.regime {
                TectonicRegime::Episodic => (EventKind::GlobalResurfacing, style.lowland_age_ma, "the lid founders and lava renews the plains"),
                _ => (EventKind::CrustalDichotomy, style.highland_age_ma - 100.0, "lowlands form in one hemisphere"),
            };
            let mut events = Table::new();
            events.insert(Key(1), Event { t_ma, kind, at: Place::Point { at: pole }, magnitude: style.lowland_fraction, causes: Vec::new(), note: note.into() });
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
    fn mars_has_a_hemispheric_dichotomy() {
        let body = Archetype::Mars.context(Key(1), "Mars");
        let grid = GridSpec::Equirect { width: 180, height: 90 }.build().unwrap();
        let ctx = RunContext { root_seed: geoforge_core::Seed(3), body: &body, grid: &grid };
        let out = LidTectonics::run(&ctx, &Products::new(), &LidTectonicsParams::default()).unwrap();
        let crust = &out.crust_type().unwrap().values;
        let area = |f: &dyn Fn(usize) -> bool| grid.cells().filter(|&c| f(c as usize)).map(|c| grid.area_sr(c)).sum::<f64>() / (4.0 * std::f64::consts::PI);
        let low = area(&|i| crust[i] == CrustType::Oceanic);
        assert!((low - 0.35).abs() < 0.03, "lowland share {low}");
        // Lowlands cluster around the dichotomy pole.
        let pole = out.plates().unwrap().table.rows[&Key(1)].seed_point;
        let near = area(&|i| crust[i] == CrustType::Oceanic && grid.center(i as Cell).dot(pole) > 0.0);
        assert!(near > 0.8 * low, "{near} of {low} near the pole");
        let age = &out.crust_age().unwrap().values;
        let mean = |t: CrustType| crust.iter().zip(age).filter(|(c, _)| **c == t).map(|(_, a)| f64::from(*a)).sum::<f64>() / crust.iter().filter(|c| **c == t).count() as f64;
        assert!(mean(CrustType::Continental) > 4000.0 && mean(CrustType::Oceanic) < mean(CrustType::Continental));
    }
}
