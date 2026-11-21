//! Every generator the pipeline knows, and how each layer picks one for a body.
//!
//! A layer is a slot. Its implementations are listed in order of preference; the first whose
//! preconditions hold for the body fills it. Implementations of one layer produce the same
//! products, so later layers never care which one ran.

use geoforge_core::product::{Generator, LayerId, Precondition, ProductId, Products, RunContext};
use geoforge_core::BodyContext;
use geoforge_surface as surface;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Implementation {
    Tectonics,
    LidTectonics,
    Geology,
    LidGeology,
    Elevation,
    LidElevation,
    Circulation,
    Climate,
    Erosion,
    Hydrology,
    Biomes,
    GiantBands,
}

/// Run `$body` with `$g` bound to the generator type of implementation `$i`.
macro_rules! dispatch {
    ($i:expr, $g:ident => $body:expr) => {
        match $i {
            Implementation::Tectonics => {
                type $g = surface::Tectonics;
                $body
            }
            Implementation::LidTectonics => {
                type $g = surface::LidTectonics;
                $body
            }
            Implementation::Geology => {
                type $g = surface::Geology;
                $body
            }
            Implementation::LidGeology => {
                type $g = surface::LidGeology;
                $body
            }
            Implementation::Elevation => {
                type $g = surface::Elevation;
                $body
            }
            Implementation::LidElevation => {
                type $g = surface::LidElevation;
                $body
            }
            Implementation::Circulation => {
                type $g = surface::Circulation;
                $body
            }
            Implementation::Climate => {
                type $g = surface::Climate;
                $body
            }
            Implementation::Erosion => {
                type $g = surface::Erosion;
                $body
            }
            Implementation::Hydrology => {
                type $g = surface::Hydrology;
                $body
            }
            Implementation::Biomes => {
                type $g = surface::Biomes;
                $body
            }
            Implementation::GiantBands => {
                type $g = surface::GiantBands;
                $body
            }
        }
    };
}

impl Implementation {
    /// In order of preference within each layer.
    pub const ALL: [Self; 12] = [
        Self::Tectonics,
        Self::LidTectonics,
        Self::Geology,
        Self::LidGeology,
        Self::Elevation,
        Self::LidElevation,
        Self::Circulation,
        Self::Climate,
        Self::Erosion,
        Self::Hydrology,
        Self::Biomes,
        Self::GiantBands,
    ];

    #[must_use]
    pub fn name(self) -> &'static str {
        dispatch!(self, G => G::NAME)
    }
    #[must_use]
    pub fn version(self) -> u32 {
        dispatch!(self, G => G::VERSION)
    }
    #[must_use]
    pub fn layer(self) -> LayerId {
        dispatch!(self, G => G::LAYER)
    }
    #[must_use]
    pub fn needs(self) -> &'static [ProductId] {
        dispatch!(self, G => G::needs())
    }
    #[must_use]
    pub fn wants(self) -> &'static [ProductId] {
        dispatch!(self, G => G::wants())
    }
    /// # Errors
    /// Why this implementation does not apply to the body.
    pub fn applies(self, body: &BodyContext) -> Result<(), Precondition> {
        dispatch!(self, G => G::applies(body))
    }
    #[must_use]
    pub fn default_params(self) -> serde_json::Value {
        strip_nulls(dispatch!(self, G => serde_json::to_value(<G as Generator>::Params::default()).expect("params serialize")))
    }
    /// Stored parameters with defaults filled in, normalised.
    ///
    /// # Errors
    /// If the stored value does not parse as this implementation's parameters.
    pub fn resolve_params(self, stored: Option<&serde_json::Value>) -> crate::Result<serde_json::Value> {
        dispatch!(self, G => {
            let p: <G as Generator>::Params = match stored {
                Some(v) => serde_json::from_value(v.clone()).map_err(|e| crate::PipelineError::Invalid(format!("[params.{}]: {e}", self.name())))?,
                None => Default::default(),
            };
            Ok(serde_json::to_value(p)?)
        })
    }
    /// # Errors
    /// From the generator.
    pub fn run(self, ctx: &RunContext<'_>, inputs: &Products, params: &serde_json::Value) -> crate::Result<Products> {
        dispatch!(self, G => {
            let p: <G as Generator>::Params = serde_json::from_value(params.clone())?;
            Ok(G::run(ctx, inputs, &p)?)
        })
    }
}

/// The implementations of `layer`, in order of preference.
pub fn candidates(layer: LayerId) -> impl Iterator<Item = Implementation> {
    Implementation::ALL.into_iter().filter(move |i| i.layer() == layer)
}

/// The implementation that fills `layer` on this body.
///
/// # Errors
/// Why none applies, one reason per candidate.
pub fn select(layer: LayerId, body: &BodyContext) -> Result<Implementation, Precondition> {
    let mut why = Vec::new();
    for i in candidates(layer) {
        match i.applies(body) {
            Ok(()) => return Ok(i),
            Err(p) => why.push(format!("{} {}", i.name(), p.0)),
        }
    }
    Err(Precondition(why.join("; ")))
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

#[cfg(test)]
mod tests {
    use super::*;
    use geoforge_core::{Archetype, Key};

    #[test]
    fn names_are_unique_and_match_serde() {
        let mut seen = std::collections::BTreeSet::new();
        for i in Implementation::ALL {
            assert!(seen.insert(i.name()), "duplicate name {}", i.name());
            assert_eq!(serde_json::to_value(i).unwrap(), serde_json::Value::from(i.name()));
        }
    }

    #[test]
    fn every_layer_has_an_implementation() {
        for l in LayerId::ALL {
            assert!(candidates(l).next().is_some(), "{} has no implementation", l.name());
        }
    }

    #[test]
    fn selection_follows_the_body() {
        let earth = Archetype::Earth.context(Key(1), "Earth");
        for l in LayerId::ALL {
            assert_eq!(select(l, &earth).is_ok(), l != LayerId::Bands, "{} on Earth", l.name());
        }
        let mars = Archetype::Mars.context(Key(1), "Mars");
        assert_eq!(select(LayerId::Tectonics, &mars), Ok(Implementation::LidTectonics));
        assert_eq!(select(LayerId::Elevation, &mars), Ok(Implementation::LidElevation));
        assert_eq!(select(LayerId::Climate, &mars), Ok(Implementation::Climate));
        assert!(select(LayerId::Hydrology, &mars).is_err(), "no surface water");
        let giant = Archetype::GasGiant.context(Key(1), "Jove");
        let why = select(LayerId::Tectonics, &giant).unwrap_err();
        assert!(why.0.contains("tectonics"), "{why}");
        assert_eq!(select(LayerId::Bands, &giant), Ok(Implementation::GiantBands));
    }
}
