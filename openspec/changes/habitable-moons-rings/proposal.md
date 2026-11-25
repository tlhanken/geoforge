# habitable-moons-rings

## Why

The default world of a universe was chosen for its water, not its warmth: four of five seeds
gave frozen worlds (seed 12345: −4 °C, 96% of its land ice or tundra). There was no way to
find systems worth a look in the viewer. Ice giants had too few moons and all moon systems
were evenly shared; nothing was captured. Generated giants had no rings, and the night sky
had no zodiacal light.

## What Changes

- **A climate estimate before generating** (`geoforge_cosmos::worlds::climate_estimate`):
  the latitude-band energy balance the climate layer is built on, its mean taken 2.5 K
  colder for the ground standing above the sea. It follows the generated climate within ~2 K
  (seven worlds, −6 to +20 °C). `habitability()` grades a world: temperate (mean 5–25 °C,
  three quarters of the surface above freezing), habitable (open water somewhere, short of
  the moist greenhouse), or none.
- **A temperate default world**: the default star rule asks first for a Sun-like star whose
  most Earth-like body is temperate (then wet, as before). The Earth-like score weighs
  habitability. Saved universes keep their focus.
- **Habitable systems flagged in the viewer**: the star view circles systems with a world of
  open water (green temperate, blue other) and counts them; the system view marks bodies
  (✓, ~); a body's panel gives its climate estimate. Worked out in parallel, cached.
  `geoforge show` prints the estimate.
- **Moons** (new `geoforge_cosmos::moons`, shared by both system builders): gas giants 3–6
  regular moons, ice giants 3–5; mass shared over 2.5 decades (one may dominate); 1–3 small
  inner moons; a third of ice giants a captured, backward giant moon (Triton) leaving one
  small moon; every giant's irregular moons as a population; captured specks for rocky
  planets just inside an asteroid belt (Phobos, Deimos).
- **Rings**: every giant a faint dust ring; a quarter of those under 150 K bright ice
  (optical depth 0.2–1.6, inside the Roche limit), drawn from their own stream.
  `Population::optical_depth`. The viewer draws them, banded, tilted with the planet.
- **Zodiacal light** on the star chart, by the system's belts against the Solar System's
  (`zodiacal_dust`, `Sky::zodiacal`), with the gegenschein.

## Impact

- `geoforge-core`: `PowerLaw::mean_square_km2`, `Population::optical_depth`.
- `geoforge-cosmos`: `worlds` (estimate, habitability), `resolve` (default rule), `galactic`
  (score, moons, rings, focus fallback to belts round stars), `moons` (new), `system`.
- `geoforge-almanac`: `Sky::zodiacal_dust`, `Sky::zodiacal`, `zodiacal_dust`.
- `geoforge-store`: `render_sky` takes the zodiacal grid.
- `geoforge-viewer`: `/api/habitable`, `Viewer::habitability`, rings and climate in scenes.
- New universes resolve to other default bodies; generated systems gain moons and rings
  (systems built before keep theirs as stored).
