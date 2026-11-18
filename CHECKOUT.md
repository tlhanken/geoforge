# Running geoforge v2 locally

```bash
git fetch origin
git checkout v2-implementation
cargo build --release -p geoforge-cli      # binary: target/release/geoforge-v2
./scripts/v2-generate.sh 42                 # world + renders + stats under ./outputs/world-42
```

| Symptom | Cause | Fix |
|---|---|---|
| `package geoforge-cli not found` | Wrong branch | `git checkout v2-implementation` |
| `package geoforge_v2 not found` | The old pre-tectonic CLI was replaced by `geoforge-cli` | Use `geoforge-v2` (see `docs/V2_PIPELINE.md`) |
| Something else | — | `./scripts/doctor.sh` |

User guide: [docs/V2_PIPELINE.md](docs/V2_PIPELINE.md). Realism check: [docs/V2_REALISM.md](docs/V2_REALISM.md).
