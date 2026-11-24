//! The hydrology layer of a world whose rain is methane (Titan): rivers, lakes and seas, and
//! their watersheds, by the same drainage as water's ([`crate::hydrology`]).
//!
//! Methane's cycle is moved by the cold trap. The air carries methane from where the ground is
//! warmest to where it is coldest and it rains out there: Titan's tropics, 2–3 K warmer than
//! its poles, lose more by evaporation than they get, and its poles gain (Mitchell 2008;
//! Schneider et al. 2012), so the seas and nearly all the lakes lie poleward of 55° and the low
//! latitudes are dry, with dune fields and dry riverbeds (Hayes 2016). The rain and the demand
//! of evaporation follow the ground's temperature against the world's mean, scaled by how far
//! it ranges: R·C^(−z) and D·C^z, z = (T − T̄)/half the range, with C = 3, R = 20 and D = 60
//! mm/yr of liquid methane (the order of Titan's mean rain in GCMs, Lora & Mitchell 2015), so
//! ground colder than the mean by half its half-range gets as much as it gives back, and only
//! the coldest is wet. With no sea to drain to, every basin ends in a lake or a dry pan.
//! Methane lakes neither freeze nor silt up in the model, so every closed depression that
//! fills keeps its lake (Titan's sit in sharp-rimmed karstic basins), until the world's liquid
//! runs out: its lakes hold a layer 0.85 m deep over the world (Titan's ~70,000 km³), filled
//! coldest first. That leaves Titan ~1% lakes, all poleward of 55° (Titan: 1.1%, Hayes 2016).
//! Not modelled: the seasons' shifts of the rain (the summer pole's storms), the orbit's
//! north–south asymmetry that keeps Titan's seas in the north (Aharonson et al. 2009), and
//! ethane and nitrogen in the liquid.

use crate::hydrology::{drain, HydrologyParams, Water};
use geoforge_core::product::{Generator, LayerId, Precondition, ProductId, Products, RunContext};
use geoforge_core::{BodyContext, Error};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct MethaneCycleParams {
    /// Rain over ground at the world's mean temperature, mm/yr of liquid methane.
    pub rain_mm: f64,
    /// Evaporative demand over ground at the world's mean temperature, mm/yr.
    pub demand_mm: f64,
    /// How many times the mean rain the coldest ground gets, and how many times less demand
    /// (and the warmest the reverse).
    pub contrast: f64,
    /// Liquid on the surface, as a layer over the whole world, m: Titan's seas and lakes hold
    /// ~70,000 km³ (Hayes 2016), 0.85 m.
    pub inventory_m: f64,
    /// The drainage, as water's.
    pub network: HydrologyParams,
}

impl Default for MethaneCycleParams {
    fn default() -> Self {
        // Titan's rivers are small (Vid Flumina carries ~10–100 m³/s in GCM storms).
        Self {
            rain_mm: 20.0,
            demand_mm: 60.0,
            contrast: 3.0,
            inventory_m: 0.85,
            network: HydrologyParams { min_river_m3_s: 5.0, min_lake_depth_m: 10.0, arid_channel_loss_per_1000km: 0.9, ..HydrologyParams::default() },
        }
    }
}

/// Rain and evaporative demand (mm/yr) of ground at `t` with the world's mean `mean` and half
/// range `half_range` (same units).
#[must_use]
pub fn rain_and_demand(t: f64, mean: f64, half_range: f64, p: &MethaneCycleParams) -> (f64, f64) {
    let z = ((t - mean) / half_range.max(0.5)).clamp(-1.0, 1.0);
    (p.rain_mm * p.contrast.powf(-z), p.demand_mm * p.contrast.powf(z))
}

pub struct MethaneCycle;

impl Generator for MethaneCycle {
    type Params = MethaneCycleParams;
    const NAME: &'static str = "methane_cycle";
    const VERSION: u32 = 1;
    const LAYER: LayerId = LayerId::Hydrology;

    fn needs() -> &'static [ProductId] {
        &[ProductId::Elevation, ProductId::Temperature]
    }
    fn wants() -> &'static [ProductId] {
        &[ProductId::Provinces, ProductId::RoughElevation]
    }
    fn applies(body: &BodyContext) -> Result<(), Precondition> {
        crate::preconditions::surface_methane(body)
    }

    fn run(ctx: &RunContext<'_>, inputs: &Products, p: &MethaneCycleParams) -> Result<Products, Error> {
        if inputs.rivers().is_some() && inputs.lakes().is_some() && inputs.watersheds().is_some() {
            return Ok(Products::new());
        }
        let elev = &inputs.elevation().ok_or(Error::MissingInput(ProductId::Elevation))?.values;
        let temp = &inputs.temperature().ok_or(Error::MissingInput(ProductId::Temperature))?.values;
        let temp = crate::climate::on_ground(temp, inputs.rough_elevation().map(|f| f.values.as_slice()), elev, ctx.body);
        let grid = ctx.grid;
        let area: Vec<f64> = (0..grid.len()).map(|c| grid.area_sr(c as geoforge_core::grid::Cell)).collect();
        let mean = temp.iter().zip(&area).map(|(&t, a)| f64::from(t) * a).sum::<f64>() / area.iter().sum::<f64>();
        let (lo, hi) = temp.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &t| (lo.min(f64::from(t)), hi.max(f64::from(t))));
        let (precip, pet): (Vec<f64>, Vec<f64>) = temp.iter().map(|&t| rain_and_demand(f64::from(t), mean, 0.5 * (hi - lo), p)).unzip();
        let r = ctx.radius_km();
        let inventory_km3 = p.inventory_m / 1000.0 * 4.0 * std::f64::consts::PI * r * r;
        Ok(drain(ctx, inputs, &Water { elev, precip: &precip, pet: &pet, lake_temp: None, sea: false, inventory_km3: Some(inventory_km3) }, &p.network))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_cold_trap_wets_the_coldest_ground() {
        let p = MethaneCycleParams::default();
        let (pole, eq, mid) = (rain_and_demand(-182.0, -180.0, 2.0, &p), rain_and_demand(-178.0, -180.0, 2.0, &p), rain_and_demand(-180.0, -180.0, 2.0, &p));
        assert!((pole.0 / pole.1 - 3.0).abs() < 1e-9 && (eq.0 / eq.1 - 1.0 / 27.0).abs() < 1e-9);
        assert!((mid.0 - 20.0).abs() < 1e-9 && (mid.1 - 60.0).abs() < 1e-9);
        // Ground half a half-range colder than the mean gives back what it gets.
        let edge = rain_and_demand(-181.0, -180.0, 2.0, &p);
        assert!((edge.0 - edge.1).abs() < 1e-9);
        // A world with hardly any range is treated as half a degree wide.
        assert!(rain_and_demand(-180.1, -180.0, 0.0, &p).0 < rain_and_demand(-180.0, -180.0, 0.0, &p).0 * 1.3);
    }
}
