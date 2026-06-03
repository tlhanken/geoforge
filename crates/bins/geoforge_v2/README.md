# geoforge_v2

CLI for Geoforge v2 cosmology (and future pipeline stages).

## Generate and inspect

By default, each run writes **JSON** under `outputs/` (no SQLite):

```bash
cargo run -p geoforge_v2 -- generate --seed 42
# → outputs/cosmology_seed42.json
```

### Export formats

| Flag | Files |
|------|--------|
| `--export json` (default) | `cosmology_seed<N>.json` |
| `--export text` | `cosmology_seed<N>.txt` |
| `--export both` | JSON + text |
| `--no-export` | stdout only |

```bash
cargo run -p geoforge_v2 -- generate --seed 42 --export text -d ./my_outputs
cargo run -p geoforge_v2 -- generate --seed 42 --cosmology-scale rich --export both -v
```

Open the JSON in any editor or `jq` — it includes galaxies, all loaded stellar markers (positions in ly), and the primary `SolarSystem` (stars, planet slots, belts).

## Other flags

- `--cosmology-scale minimal|rich|expansive`
- `--to-layer solar_system` (default)
- `-v` verbose stdout
