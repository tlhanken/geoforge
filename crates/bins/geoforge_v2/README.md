# geoforge_v2

CLI for the v2 **pre-tectonic** pipeline (universe → planetary body). Tectonics and surface stages are not in this binary yet.

## Quick review (recommended)

From the repository root:

```bash
chmod +x scripts/v2-generate.sh
./scripts/v2-generate.sh
```

Writes **`outputs/cosmology_seed42.json`** and **`outputs/cosmology_seed42.txt`** (seed 42, full pipeline through `planetary_body`).

Custom seed or directory:

```bash
./scripts/v2-generate.sh 99 ./outputs both
./scripts/v2-generate.sh 42 ./inspect text -v
```

## Manual `cargo` command

Always put `--` before `generate` so Cargo does not swallow flags:

```bash
cargo build -p geoforge_v2

cargo run -p geoforge_v2 -- generate \
  --seed 42 \
  --format both \
  --output-directory outputs \
  -v
```

## Defaults

| Flag | Default |
|------|---------|
| `--seed` | `0x60F06E55EED` if omitted |
| `--to-layer` | `planetary_body` |
| `--from-layer` | `universe` |
| `--format` | `both` (JSON + text) |
| `-d` / `--output-directory` | `outputs` |
| `--cosmology-scale` | `minimal` |

## Output files

| File | Contents |
|------|----------|
| `outputs/cosmology_seed<N>.json` | Full report (galaxies, markers, solar system, planetary bodies) |
| `outputs/cosmology_seed<N>.txt` | Same data, human-readable |

The `outputs/` directory is gitignored; generated files are for local review only.

## Tests

```bash
cargo test -p geoforge_v2 --test cli_smoke
cargo test -p geoforge-export
```

See also [docs/V2_PIPELINE.md](../../docs/V2_PIPELINE.md) and [crates/libs/README.md](../../libs/README.md).
