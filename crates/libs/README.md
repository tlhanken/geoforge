# Geoforge v2 libraries

## Layout

| Crate | Role |
|-------|------|
| `geoforge-types` | Layer structs, enums, seeds (no generation) |
| `stage1_universe` … `stage5_planetary_body` | One `generate()` per pipeline step |
| `geoforge-export` | Run stages 1–5, query helpers, JSON/text reports |
| `Materials` | Legacy materials data (not in workspace yet) |

Empty directories `stage6_tectonics`, `stage7_geologic_domains`, `stage8_heightmap` are placeholders for the planetary surface pipeline.

## Naming

Folder `stageN_*` follows **pipeline layer order** (see `PipelineLayerId`), not the old v1 “Stage 1 = tectonics” README wording.

| Folder | `PipelineLayerId` |
|--------|-------------------|
| `stage1_universe` | `universe` |
| `stage2_supercluster` | `supercluster` |
| `stage3_galaxy` | `galaxy` |
| `stage4_solar_system` | `solar_system` |
| `stage5_planetary_body` | `planetary_body` |
| `stage6_tectonics` | `tectonics` (future) |

## Dependency rule

Each stage crate depends only on `geoforge-types` and (for stage 4→5) the prior stage crate when shared helpers are needed (`stage4` → `stage3::marker_position`).

`geoforge-export` may depend on all stage crates to orchestrate and export; bins depend on `geoforge-export` + `geoforge-types`.
