# geoforge-types

Core domain types for Geoforge v2: seeds, coordinates, cosmology, planetary parameters, tectonics, geology, and pipeline layer IDs.

**Data-only** — no generators, no rasters, no I/O.

## Cosmology (Stage 0)

- **`GalacticPoint`** — absolute positions in light-years for rendering
- **`RegionId`** — cylindrical (disk/spiral/ring) or spherical (bubble/elliptical) cells
- **`StellarSystemMarker`** — galaxy-map LOD (position + phenotype, no interiors)
- **`SolarSystem`** — zoom LOD (stars, planet slots, belts)
- **`CosmologyPreset`** — `Minimal` (1 system), `Rich` (local region), `Expansive`

Regenerate everything from **`Seed`** + preset; no JSON persistence required.

## Planetary pipeline (Stages 1+)

See module docs for `tectonics`, `geology`, `planetary`, `pipeline`.
