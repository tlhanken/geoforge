# Geoforge v2 libraries

| Crate | Role | Depends on |
|---|---|---|
| `geoforge-core` | Data and the body's own physics: keys, seeds, geometry, the `Grid` trait and `EquirectGrid`, the five shapes (`Table`, `Field`, `Partition`, `Network`, `Timeline`), body and system types (`body/`: system, matter, interior with `TectonicRegime`, context, `Archetype` presets), derived quantities and the gray-atmosphere `radiation` model, surface entity types (rock and biome catalog rows), product ids, sources, `Generator` trait | serde |
| `geoforge-cosmos` | Thin cosmology: one solar system per seed (`system`), what each body is made of (`worlds`), and the `BodyContext` handoff for any body, belt members included (`context`) | core |
| `geoforge-surface` | The surface layers as `Generator`s: plate tectonics, geology (with the rock catalog), elevation (with erosion); `lid/` generators for rigid-lid worlds; the shared `atmosphere/` model; circulation, climate, hydrology, biomes (with the biome catalog); `preconditions` that pick between them | core |
| `geoforge-store` | The world directory: manifest, lossless product storage with content hashes, edit views (PNG + world file + sidecar), patches, preview renders | core |
| `geoforge-pipeline` | The generator `registry` (layer slots), plan, provenance, incremental generation, export/import/reset, realism stats | all of the above |

Binaries: `crates/bins/geoforge-cli` (`geoforge-v2`) and `crates/bins/geoforge` (v1, kept as
the reference for parity tests).

Rules:

- `geoforge-core` has no generation algorithms and no I/O; its physics is pure functions of
  a body (derived quantities, radiative balance) that both cosmology and the surface use.
- A generator never branches on the kind of body: a different kind of body gets a different
  generator for the same layer, chosen by its preconditions.
- Surface layers see a `BodyContext` and earlier products, never cosmology.
- Algorithms ask a `Grid` for positions, areas and neighbours; they never use pixel x/y, so a
  Goldberg grid can be added without touching them.
- Only `geoforge-store` reads or writes files.

See `docs/V2_PIPELINE.md` for the user guide.
