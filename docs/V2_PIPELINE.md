# Geoforge v2 pipeline

## Review workflow

```bash
./scripts/v2-generate.sh
```

Then open:

- `outputs/cosmology_seed42.txt` — readable summary
- `outputs/cosmology_seed42.json` — full structured data for tools

## Layer order (pre-tectonics)

| Layer ID | Crate | Output type |
|----------|--------|-------------|
| `universe` | `stage1_universe` | `Universe` |
| `supercluster` | `stage2_supercluster` | `Supercluster` |
| `galaxy` | `stage3_galaxy` | `GalaxyCatalog` |
| `solar_system` | `stage4_solar_system` | `SolarSystem` |
| `planetary_body` | `stage5_planetary_body` | `PlanetaryBodies` |
| `tectonics` | `stage6_tectonics` | *(planned)* |

Folder names `stageN_*` match this table, not the v1 README where “Stage 1” meant tectonics.

## Crate layout

- **`geoforge-types`** — structs only (`CosmologyLayers`, `PlanetaryParams`, …)
- **`stage1` … `stage5`** — `generate()` per layer
- **`geoforge-export`** — `run_pre_tectonic_layers`, reports, file I/O
- **`geoforge_v2`** — CLI

## Presets

- `minimal` — one galaxy, one primary system (fast review)
- `rich` — more markers near focus
- `expansive` — multiple galaxies

```bash
cargo run -p geoforge_v2 -- generate --seed 42 --cosmology-scale rich -d outputs
```
