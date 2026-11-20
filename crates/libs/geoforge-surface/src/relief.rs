//! How tall topography a body holds.
//!
//! Rock at the foot of a mountain or plateau bears a load of `ρ g h`. When that exceeds what
//! the crust can carry, it flows or breaks until the load fits, so the tallest relief scales
//! as `1/g` (Mars at 0.38 g holds Olympus Mons, 2.5× Mauna Loa). Landform heights are set
//! for Earth and multiplied by [`relief_scale`]; widths stay in km.

use geoforge_core::BodyContext;

/// Relief multiplier from gravity, Earth = 1: weaker gravity holds taller topography.
#[must_use]
pub fn relief_scale(body: &BodyContext) -> f64 {
    (9.81 / body.gravity_m_s2()).clamp(0.5, 3.0)
}
