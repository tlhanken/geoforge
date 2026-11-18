//! Static dispatch over the layers: one place that knows every generator.

use geoforge_core::product::{Generator, LayerId, Precondition, ProductId, Products, RunContext};
use geoforge_core::BodyContext;
use geoforge_surface::{Biomes, Climate, Elevation, Geology, Hydrology, Tectonics};

macro_rules! dispatch {
    ($l:expr, $g:ident => $body:expr) => {
        match $l {
            LayerId::Tectonics => {
                type $g = Tectonics;
                $body
            }
            LayerId::Geology => {
                type $g = Geology;
                $body
            }
            LayerId::Elevation => {
                type $g = Elevation;
                $body
            }
            LayerId::Climate => {
                type $g = Climate;
                $body
            }
            LayerId::Hydrology => {
                type $g = Hydrology;
                $body
            }
            LayerId::Biomes => {
                type $g = Biomes;
                $body
            }
        }
    };
}

/// Bump a layer's version when its algorithm changes output for the same inputs; stored
/// products of that layer (and everything downstream) then go stale.
#[must_use]
pub fn version(l: LayerId) -> u32 {
    match l {
        LayerId::Tectonics => 1,
        LayerId::Geology => 1,
        LayerId::Elevation => 1,
        LayerId::Climate => 1,
        LayerId::Hydrology => 1,
        LayerId::Biomes => 1,
    }
}

#[must_use]
pub fn needs(l: LayerId) -> &'static [ProductId] {
    dispatch!(l, G => G::needs())
}

#[must_use]
pub fn wants(l: LayerId) -> &'static [ProductId] {
    dispatch!(l, G => G::wants())
}

/// # Errors
/// Why the layer does not apply to the body.
pub fn applies(l: LayerId, body: &BodyContext) -> Result<(), Precondition> {
    dispatch!(l, G => G::applies(body))
}

#[must_use]
pub fn default_params(l: LayerId) -> serde_json::Value {
    strip_nulls(dispatch!(l, G => serde_json::to_value(<G as Generator>::Params::default()).expect("params serialize")))
}

/// TOML has no null: drop null members (they mean "default").
#[must_use]
pub fn strip_nulls(v: serde_json::Value) -> serde_json::Value {
    match v {
        serde_json::Value::Object(m) => serde_json::Value::Object(m.into_iter().filter(|(_, v)| !v.is_null()).map(|(k, v)| (k, strip_nulls(v))).collect()),
        serde_json::Value::Array(a) => serde_json::Value::Array(a.into_iter().map(strip_nulls).collect()),
        o => o,
    }
}

/// Stored parameters with defaults filled in, normalised.
///
/// # Errors
/// If the stored value does not parse as the layer's parameters.
pub fn resolve_params(l: LayerId, stored: Option<&serde_json::Value>) -> crate::Result<serde_json::Value> {
    dispatch!(l, G => {
        let p: <G as Generator>::Params = match stored {
            Some(v) => serde_json::from_value(v.clone()).map_err(|e| crate::PipelineError::Invalid(format!("[params.{}]: {e}", l.name())))?,
            None => Default::default(),
        };
        Ok(serde_json::to_value(p)?)
    })
}

/// # Errors
/// From the generator.
pub fn run(l: LayerId, ctx: &RunContext<'_>, inputs: &Products, params: &serde_json::Value) -> crate::Result<Products> {
    dispatch!(l, G => {
        let p: <G as Generator>::Params = serde_json::from_value(params.clone())?;
        Ok(G::run(ctx, inputs, &p)?)
    })
}
