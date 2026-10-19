# geoforge-export

Serialize and write pipeline outputs. **No generation logic** beyond orchestrating stage crates.

## Modules

- **`pipeline`** — `run_cosmology_layers`, query helpers (`galaxy`, `markers_near_focus`, …)
- **`report`** — `CosmologyReport`, JSON/text export

Layer data lives in **`geoforge-types`** (`CosmologyLayers`, `PlanetaryBodies`, …).

## Example

```rust
use geoforge_export::{run_cosmology_layers, build_cosmology_report, write_cosmology_report, ExportFormat};
use geoforge_types::{PipelineLayerId, Seed};
use geoforge_types::cosmology::CosmologyPreset;

let layers = run_cosmology_layers(Seed::new(42), CosmologyPreset::Minimal, PipelineLayerId::PlanetaryBody);
let report = build_cosmology_report(&layers, Seed::new(42), CosmologyPreset::Minimal, PipelineLayerId::PlanetaryBody);
write_cosmology_report(std::path::Path::new("./inspect"), &layers, &report, ExportFormat::Both)?;
```
