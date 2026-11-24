# universe-and-figures

## Why

The night sky is drawn from the world's own galaxy, but that galaxy sits alone: no
Andromeda, no Magellanic Clouds, nothing beyond it. The cosmos should be a universe made
from a seed and zoomed into, its galaxies generated where they are looked at, with any
world's galaxy given a place in it. And a world with no constellations of its own shows a
sky with no figures at all, though people everywhere have drawn them from the brightest
stars.

## What Changes

- A universe from a seed: dark-matter halos by mass, placed in nested cubes (one size per
  half-decade of mass, as the galaxy's stars are) and clustered in a lognormal density
  field; in each halo a central galaxy and its satellites, with stellar masses, types and
  sizes from the halo's mass and its surroundings. Every galaxy has a stable identity that
  seeds it, so any observer sees the same galaxy and a query generates only the cubes it
  reaches.
- A world's galaxy has a place in a universe: the generated galaxy like it (the same broad
  class, within 0.25 dex of its stellar mass) nearest the universe's origin; its
  satellites and neighbours are that galaxy's. The Solar System's universe is the real one.
- Other galaxies in the night sky: every galaxy brighter than magnitude 8 is listed (name,
  type, magnitude, distance, size, position) and drawn into the glow by its surface
  brightness through the home galaxy's dust. The Sun sees the real Local Group and the
  bright nearby galaxies.
- Automatic constellations: a world without its own gets figures drawn from its brightest
  stars, which it can save and then rename or edit.
- CLI: `geoforge universe` lists the galaxies around a world's own; `geoforge constellation
  auto` saves the automatic figures. The report and `geoforge sky` list the galaxies and
  draw the automatic figures.

Not in this change: visiting another galaxy's stars (the zoom from a universe galaxy into
its sectors, which needs a galaxy model scaled to the generated galaxy), superclusters'
filaments beyond the lognormal field, and a world holding more than one system.

## Capabilities

### New Capabilities
- `universe`: halos and galaxies generated in nested cubes; a world's galaxy placed in it.

### Modified Capabilities
- `night-sky`: other galaxies in the sky and its glow; automatic constellations.

## Impact

- `geoforge-core`: `GalaxyId`, `UniverseSite`; `GalacticSite.universe` (optional; old
  systems get a universe seeded from their galaxy's seed).
- `geoforge-cosmos`: new `universe` (`halos`, `members`, `sky`) and `constellations`
  modules; `sky::Neighbourhood` keeps the galactic pole's brightness and the extinction out
  of the galaxy along each ray; `sol` marks the real universe.
- `geoforge-almanac`: `Sky.galaxies`, the galaxies in `Sky::glow`, automatic constellations.
- `geoforge-pipeline`: the report's sky lists the galaxies and draws automatic figures.
- CLI: `universe`, `constellation auto`; `sky` lists galaxies.
- Worlds made before this change: Earth worlds (`--archetype earth`) must be made again to
  see the real Local Group; until then they see a generated universe.
