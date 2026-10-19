# v2-architecture

## Why

Geoforge's surface pipeline grew organically: geology as struct/vector overlays, fuzzy stage
boundaries, recomputed intermediates, and no way to hand-author or edit a layer. v2 replaces
this with a strictly linear, per-pixel pipeline of named layers, each of which can be
generated from a seed or injected from an editable file, over a shared world container.

## What Changes

- Introduces a canonical 16-layer pipeline: six cosmology layers (an "interior of the
  container above" descent) plus ten surface layers on the focus body's grid.
- Defines the world as a directory: a manifest, a path-namespaced cosmology tree, and
  per-body layer files. Rasters are PNG + JSON sidecar; structs are JSON/TOML.
- Defines the CLI: per-layer and whole-pipeline generation, layer insert/export, driven by a
  seed or by injected data.
- Establishes determinism as a hard contract over `(seed, injected layers, code version)`.
- Reserves `culture` and `religion` layers as unspecified placeholders.

## Capabilities

- `cli` — the `geoforge` command surface
- `world-container` — on-disk world layout and manifest
- `determinism` — reproducibility contract
- `layer-applicability` — declared preconditions, skip semantics, not-applicable markers
- `layer-universe`, `layer-superclusters`, `layer-galaxies`, `layer-solar-systems`,
  `layer-solar-system`, `layer-planetary-body` — cosmology layers
- `layer-tectonics`, `layer-geology`, `layer-elevation`, `layer-climate`,
  `layer-hydrology`, `layer-biomes`, `layer-resources`, `layer-hazards`,
  `layer-settlements`, `layer-kingdoms` — surface layers
- `layer-culture`, `layer-religion` — reserved placeholders (no design detail)

## Impact

- New crate layout (`geoforge-types` as the shared-types sink, logic crates, and a separate
  CLI crate). Supersedes the current `stageN_*` / `geoforge-types` attempts, which become
  reference material only.
- Supersedes `PipelineLayerId` (a subset; renamed and extended to the 16-layer roster).
- Build migration to Bloomery is tracked as a separate change, not here.
