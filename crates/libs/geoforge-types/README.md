# geoforge-types

Core domain types for Geoforge v2: seeds, coordinates, cosmology, planetary parameters, tectonics, geology, and pipeline layer IDs.

**Data-only** — no stage generation, no I/O.

## Modules

- **`seed`** — [`Seed`](src/seed.rs), deterministic `draw` / `unit` / `range`
- **`cosmology`** — universe → solar system types, [`CosmologyLayers`](src/cosmology/layers.rs), seed paths
- **`planetary`** — [`PlanetaryParams`](src/planetary/params.rs), [`PlanetaryBodies`](src/planetary/body.rs) (stage 5)
- **`planetary`**, **`tectonics`**, **`geology`**, **`pipeline`**

Generation lives in `stage*` crates; export in `geoforge-export`; CLIs in `crates/bins/`.
