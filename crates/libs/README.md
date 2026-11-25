# Geoforge libraries

| Crate | Role | Depends on |
|---|---|---|
| `geoforge-core` | Data and the body's own physics: keys, seeds, geometry, the `Grid` trait and `EquirectGrid`, the five shapes (`Table`, `Field`, `Partition`, `Network`, `Timeline`), body and system types, derived quantities, the gray-atmosphere radiation model, magnetism, surface entity types, product ids, sources, the `Generator` trait | serde |
| `geoforge-names` | Pronounceable names from a seed | — |
| `geoforge-cosmos` | Galaxies and their populations, stellar evolution, solar systems and what each body is made of, the `BodyContext` handoff, the stellar neighbourhood of a system | core, names |
| `geoforge-surface` | The surface layers as `Generator`s: tectonics, geology (rock catalog), elevation, erosion, rigid-lid worlds, the shared atmosphere model, circulation, climate, hydrology, biomes (biome catalog), giants' bands; preconditions that pick between them | core |
| `geoforge-analysis` | What a surface means: Köppen classes, soils, timber and ore deposits (with the giant-deposit generator), natural hazards, sea ice, reefs and productivity, sea-level pressure | core, surface |
| `geoforge-store` | The world directory: manifest, lossless products with content hashes, edit views (PNG + world file + sidecar, scaled, polar caps), patches, preview renders, the derived cache | core |
| `geoforge-almanac` | The world in time: Kepler positions and calendars, eclipses, the night sky, the storm field, daily weather and its phenomena, hazards near a place | core, cosmos, analysis, store, names |
| `geoforge-pipeline` | The generator registry, plan, provenance, incremental generation, export/import/reset, derived views, place reports, realism stats | all of the above |

Programs, side by side on the pipeline: `crates/bins/geoforge-cli` (`geoforge`) and
`crates/bins/geoforge-viewer` (`geoforge-viewer`, the viewer in a browser).

Rules:

- `geoforge-core` has no generation algorithms and no I/O; its physics is pure functions of
  a body that both cosmology and the surface use.
- A generator never branches on the kind of body: a different kind of body gets a different
  generator for the same layer, chosen by its preconditions.
- Surface layers see a `BodyContext` and earlier products, never cosmology.
- Algorithms ask a `Grid` for positions, areas and neighbours; they never use pixel x/y, so a
  Goldberg grid can be added without touching them.
- Only `geoforge-store` reads or writes files.

See `docs/GETTING_STARTED.md` to use them through the `geoforge` command.
