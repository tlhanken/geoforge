# geoforge_v2

CLI entry point for Geoforge v2. Generation logic will live in stage crates; this binary orchestrates layers defined in `geoforge-types`.

```bash
cargo run -p geoforge_v2 -- generate --seed 42 --from-layer tectonics --to-layer geology
```
