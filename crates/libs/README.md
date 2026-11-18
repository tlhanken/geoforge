# Geoforge v2 libraries

| Crate | Role | Depends on |
|---|---|---|
| `geoforge-core` | Data only: keys, seeds, geometry, the `Grid` trait and `EquirectGrid`, the five shapes (`Table`, `Field`, `Partition`, `Network`, `Timeline`), body and system types, surface entity types, product ids, sources, `Generator` trait | serde |
| `geoforge-cosmos` | Thin cosmology: one solar system per seed, and the `BodyContext` handoff for its focus world | core |
| `geoforge-surface` | The surface layers as `Generator`s: tectonics, geology (with lithology), elevation (with erosion), circulation, climate, hydrology, biomes | core |
| `geoforge-store` | The world directory: manifest, lossless product storage with content hashes, edit views (PNG + world file + sidecar), patches, preview renders | core |
| `geoforge-pipeline` | Plan, provenance, incremental generation, export/import/reset, realism stats | all of the above |

Binaries: `crates/bins/geoforge-cli` (`geoforge-v2`) and `crates/bins/geoforge` (v1, kept as
the reference for parity tests).

Rules:

- `geoforge-core` has no algorithms and no I/O.
- Surface layers see a `BodyContext` and earlier products, never cosmology.
- Algorithms ask a `Grid` for positions, areas and neighbours; they never use pixel x/y, so a
  Goldberg grid can be added without touching them.
- Only `geoforge-store` reads or writes files.

See `docs/V2_PIPELINE.md` for the user guide.
