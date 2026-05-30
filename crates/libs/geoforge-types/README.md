# geoforge-types

Core domain types for Geoforge v2: seeds, coordinates, planetary parameters, tectonic and geologic enums, and pipeline layer identifiers.

This crate is **data-only** — no raster grids, no generators, no I/O. Downstream crates (`geoforge-grid`, `geoforge-tectonics`, etc.) depend on these types.

## Modules

| Module | Purpose |
|--------|---------|
| `seed` | Deterministic `Seed` derivation for hierarchical generation |
| `coordinates` | `LatLon`, `PixelCoord`, `MapExtent`, orbital and Euclidean points |
| `planetary` | `PlanetaryParams` presets (Earth, Mars, Venus) |
| `tectonics` | Crust, plates, boundaries, `TectonicLayerData` |
| `geology` | `GeologicProvince`, tectonic context, `GeologyLayerData` |
| `pipeline` | `PipelineLayerId` — ordered pipeline stages |
| `cosmology` | Optional large-scale stubs (stars, galaxies) for Stage 0 |

## Design rules

1. **Per-pixel layers** are described here as typed fields (e.g. `CrustType`); actual `TerrainMap` storage lives in `geoforge-grid`.
2. **Plate character ≠ crust composition** — `PlateType` is per-plate metadata; `CrustType` is per-pixel.
3. **Serde** on all public types for snapshots and tests.
4. **No `todo!()`** — every public API is implemented or explicitly documented as future.
