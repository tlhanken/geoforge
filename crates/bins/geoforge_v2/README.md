# geoforge_v2

CLI for Geoforge v2 cosmology (and future pipeline stages).

**Requires branch `cursor/v2-types-foundation-b88c` (or later merged main)** — older `cleanup_and_nix` only had a **placeholder** CLI.

## Important: SQLite was never implemented

On the old stub CLI, the terminal printed something like `Output format: Sqlite` because the README listed `-f sqlite` as a planned format. **No SQLite file was ever written** — it was echo-only scaffolding, not real export.

Current CLI writes real files: **JSON** (default), **text**, or **both**.

## Generate and inspect

```bash
# From repo root — rebuild so you pick up the latest CLI
cargo build -p geoforge_v2

cargo run -p geoforge_v2 -- generate --seed 42
# → outputs/cosmology_seed42.json
```

### Export flags (any of these work)

| Flag | Files |
|------|--------|
| `-f json` or `--export-format json` (default) | `cosmology_seed<N>.json` |
| `-f text` or `--export-format text` | `cosmology_seed<N>.txt` |
| `-f both` | JSON + text |
| `--no-export` | stdout only |

```bash
cargo run -p geoforge_v2 -- generate --seed 42 -f text -d ./inspect
cargo run -p geoforge_v2 -- generate --seed 42 --export-format both -d ./inspect
```

If you see `unexpected argument '--export'`, you are on an **old binary** — run `git pull`, checkout the v2 branch above, and `cargo build -p geoforge_v2` again. Then run `cargo run -p geoforge_v2 -- generate --help` and confirm you see `--export-format` / `-f`.

### Verify the CLI version

```bash
cargo run -p geoforge_v2 -- generate --help | grep -E 'export-format|sqlite'
```

You should see `export-format` and **no** sqlite.

## Other flags

- `--cosmology-scale minimal|rich|expansive`
- `--to-layer solar_system` (default)
- `-v` verbose stdout
