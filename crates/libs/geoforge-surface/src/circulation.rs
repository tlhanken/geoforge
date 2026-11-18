//! Circulation layer: seasonal winds, ocean currents and sea-surface temperature anomaly.
//! (Placeholder until the layer is implemented.)

use geoforge_core::product::{Generator, LayerId, Precondition, ProductId, Products, RunContext};
use geoforge_core::{BodyContext, Error};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct CirculationParams {}

pub struct Circulation;

impl Generator for Circulation {
    type Params = CirculationParams;
    const LAYER: LayerId = LayerId::Circulation;
    fn needs() -> &'static [ProductId] {
        &[ProductId::Elevation]
    }
    fn applies(_body: &BodyContext) -> Result<(), Precondition> {
        Err(Precondition("not implemented yet".into()))
    }
    fn run(_ctx: &RunContext<'_>, _inputs: &Products, _p: &CirculationParams) -> Result<Products, Error> {
        Ok(Products::new())
    }
}
