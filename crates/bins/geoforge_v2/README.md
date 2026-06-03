# geoforge_v2

CLI for Geoforge v2 cosmology. **Branch:** `cursor/v2-types-foundation-b88c` (or later).

## Copy-paste commands that work

From the **repository root**:

```bash
cargo build -p geoforge_v2

# Text report in ./inspect (recommended for reading)
cargo run -p geoforge_v2 -- generate --seed 42 --format text -d ./inspect

# JSON (default format)
cargo run -p geoforge_v2 -- generate --seed 42 -d ./outputs

# Both files
cargo run -p geoforge_v2 -- generate --seed 42 --format both -d ./inspect
```

Or use the wrapper script (same behavior):

```bash
chmod +x scripts/v2-generate.sh
./scripts/v2-generate.sh 42 ./inspect text
```

## Critical: the `--` before `generate`

`cargo` will steal flags if they appear **before** `--`:

```bash
# WRONG — Cargo error: unexpected argument '--export'
cargo run -p geoforge_v2 --export text -- generate --seed 42

# RIGHT — flags after `--` go to geoforge_v2
cargo run -p geoforge_v2 -- generate --seed 42 --format text -d ./inspect
```

## Export flags (on the binary, after `--`)

| Flag | Output |
|------|--------|
| `--format json` (default) | `cosmology_seed<N>.json` |
| `--format text` | `cosmology_seed<N>.txt` |
| `--format both` | JSON + text |
| `--export-format text` | alias of `--format` |
| `--export text` | alias of `--format` |
| `--no-export` | terminal only |

**Do not use** `-f sqlite` — SQLite was never implemented (old stub only printed the word).

## Verify your binary

```bash
cargo run -p geoforge_v2 -- generate --help | grep format
```

You should see `--format` in the help. Startup should print:

`geoforge_v2 0.2.0 — cosmology export enabled`

## Run smoke tests

```bash
cargo test -p geoforge_v2 --test cli_smoke
```
