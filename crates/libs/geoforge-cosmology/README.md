# geoforge-cosmology

Seed-driven cosmology generation for Geoforge v2.

## Usage

```rust
use geoforge_cosmology::CosmologyContext;
use geoforge_types::cosmology::CosmologyPreset;
use geoforge_types::Seed;

let ctx = CosmologyContext::new(Seed::new(42), CosmologyPreset::Minimal);

// Local region markers (no full galaxy array)
let markers = ctx.markers_near_focus(0);

// Zoom: full solar system regenerated from reference
let system = ctx.solar_system(ctx.primary_ref());
```

## Design

- **Region streaming** — only load sectors intersecting `focus ± load_radius_ly`
- **Morphology-aware partitions** — cylindrical rings×wedges×layers vs spherical shells
- **Per-system subseeds** — `root → galaxy → region → system`
- **Phenotype** on map; **stars + orbit slots** on zoom (no planetary surfaces)
