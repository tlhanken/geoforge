//! Pipeline stage 2 — supercluster from universe layer.

#![deny(clippy::todo, clippy::unimplemented)]

use geoforge_types::seed::{draw, range, unit};
use geoforge_types::coordinates::GalacticPoint;
use geoforge_types::cosmology::{Supercluster, Universe};

/// Generate supercluster framing using only [`Universe`] data.
#[must_use]
pub fn generate(universe: &Universe) -> Supercluster {
    let s = universe.root.child(b"supercluster");
    let label = format!(
        "Supercluster-{:04}",
        (unit(draw(s, b"label")) * 10_000.0) as u32
    );
    let extent_ly = universe.config.load_radius_ly * range(draw(s, b"extent"), 4.0, 12.0);
    let origin_ly = GalacticPoint::new(
        range(draw(s, b"ox"), -extent_ly * 0.1, extent_ly * 0.1),
        range(draw(s, b"oy"), -extent_ly * 0.05, extent_ly * 0.05),
        range(draw(s, b"oz"), -extent_ly * 0.1, extent_ly * 0.1),
    );
    let focus_ly = GalacticPoint::origin();

    Supercluster {
        root: universe.root,
        config: universe.config.clone(),
        label,
        extent_ly,
        origin_ly,
        focus_ly,
    }
}
