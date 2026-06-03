# Running geoforge_v2 locally

Cloud agent work is on branch **`cursor/v2-types-foundation-b88c`**. It is **not merged into `main` yet**.

## If commands fail locally

| Symptom | Cause | Fix |
|---------|--------|-----|
| `package geoforge_v2 not found` | You are on **`main`** (single crate at repo root) | Checkout v2 branch (below) |
| `unexpected argument '--export'` | Flag passed to **cargo** before `--` | Put flags after `--` (see below) |
| `unexpected argument '--export'` on generate | **Old** geoforge_v2 stub (no export) | Checkout v2 branch + rebuild |
| `edition2024` required | On **`cleanup_and_nix`** with old Cargo | Use v2 branch + `rustup update stable` |
| `Output format: Sqlite` but no file | Old stub only **printed** sqlite | Checkout v2 branch |

## One-time setup

```bash
cd /path/to/geoforge
git fetch origin
git checkout cursor/v2-types-foundation-b88c
rustup update stable   # needs Rust 1.75+; 1.83+ recommended
./scripts/doctor.sh    # must print "All checks passed"
```

## Working command

```bash
cargo run -p geoforge_v2 -- generate --seed 42 --format text -d ./inspect
#                            ^^ required — separates cargo args from program args
```

Inspect: `./inspect/cosmology_seed42.txt`

## Yes: pull from cloud / GitHub

Your **local** repo must have the same branch as the cloud VM:

```bash
git pull origin cursor/v2-types-foundation-b88c
```

Or open PR [#4](https://github.com/tlhanken/geoforge/pull/4) and merge when ready.
