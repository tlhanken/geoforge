# Geoforge v2

A world is a directory. Seeds generate every layer; any layer can be exported to files that
mapping or paint software opens, edited or drawn from scratch, and imported back. Everything
downstream then regenerates from your version.

```text
cosmology (thin: one solar system)  →  BodyContext  →  surface layers  →  (refinement, later)
                                  tectonics → geology → elevation → circulation → climate → hydrology → biomes
```

## Quick start

```bash
cargo build --release -p geoforge-cli          # binary: target/release/geoforge-v2

geoforge-v2 new worlds/demo --seed 42          # a solar system and its focus world
geoforge-v2 generate worlds/demo --render      # ~40 s at 2048×1024; previews in worlds/demo/renders
geoforge-v2 stats worlds/demo                  # summary numbers next to Earth's
```

`new --earth-like` uses an Earth-like reference body instead of a generated system (useful
for calibration). `new --body my_world.toml` uses a body you wrote (a `BodyContext` in TOML
or JSON; `bodies/<name>/body.json` of any world is a template).

## Commands

| Command | What it does |
|---|---|
| `new <dir> --seed N [--width 2048] [--earth-like \| --body file]` | Create a world |
| `generate <dir> [--render] [--from LAYER] [--until LAYER] [--force]` | Build what is missing or stale |
| `plan <dir>` | What `generate` would run, and why |
| `status <dir>` | Source (generate / patched / authored) and state of every product |
| `export <dir> [products\|layers\|all] [--to DIR]` | Write edit views (default `<dir>/edit`) |
| `import <dir> <product> [--from DIR] [--authored]` | Bring an edit back as a patch, or as authored data |
| `reset <dir> <product>` | Drop a product's edits |
| `render <dir> [products] [--to DIR]` | Preview PNGs |
| `show <dir>` | The solar system and the focus body |
| `seed <dir> <product> [--owner KEY]` | The seed a product is generated from |
| `stats <dir>` | Realism summary with Earth reference values |

## The world directory

```text
world.toml                       seed, grid, layer parameters, product sources, provenance
system.json                      the solar system (absent for a hand-written body)
bodies/<body>/body.json          BodyContext: the only thing the surface layers see
bodies/<body>/products/          generated products (cache: always rebuildable)
bodies/<body>/edits/             your authored data and patches (never overwritten)
cache/baselines/                 what each edit view showed when exported (for diffs)
cache/bases/                     generated product under each patch
edit/                            default edit-view folder
renders/                         preview images
```

Layer parameters are written out in `world.toml` with their defaults, so they can be edited.
Changing one makes exactly that layer and everything after it stale.

## Products

| Layer | Products | Shape |
|---|---|---|
| tectonics | `plates`, `crust_type`, `crust_age`, `tectonic_events` | partition, field, field, timeline |
| geology | `provinces`, `geology_events`, `lithology` (14 rock classes), `rock_age` (Ma) | partition, timeline, field, field |
| elevation | `elevation` (m above sea level) | field |
| circulation | `wind_jan`, `wind_jul` (m/s), `ocean_currents` (m/s), `sst_anomaly` (°C) | vector fields, field |
| climate | `temperature` (°C), `temperature_range` (warmest − coldest month, °C), `precipitation` (mm/yr), `precipitation_warm_share` (share in the warm half-year) | fields |
| hydrology | `rivers`, `lakes`, `watersheds` | network, partition, partition |
| biomes | `biome` | field |

Derived values (plate boundaries, slope, sea level, flow directions) are recomputed, not stored.

## Editing a layer

```bash
geoforge-v2 export worlds/demo elevation       # edit/elevation.png (+ .pgw, .prj, .edit.json)
#   … paint in GIMP, Krita, QGIS, Photoshop …
geoforge-v2 import worlds/demo elevation       # a patch: only the pixels you changed
geoforge-v2 generate worlds/demo --render      # climate, rivers and biomes follow your edit
```

Edit views:

| Product kind | View |
|---|---|
| Continuous field | 16-bit greyscale PNG; `value = offset + scale × pixel` (range in the sidecar; elevation −12 000 … 12 000 m) |
| Vector field (winds, currents) | 16-bit RGB PNG: red = eastward, green = northward, each `offset + scale × channel` |
| Category field (crust type, lithology, biome) | RGB PNG in the category palette; off-palette pixels (anti-aliasing) snap to the nearest colour |
| Partition (plates, provinces, lakes, watersheds) | RGB PNG, one colour per key, black = none, plus `<p>.table.json` (key, colour, row). Add a region by adding a row with a new key and colour |
| Timeline, network | The JSON document |

Every raster has a `.pgw` world file and a `.prj` sphere definition (the body's radius), so GIS
tools place it on the globe. The grid is equirectangular with pixel centres at
`lon = −180 + (x + ½)·360/W`, `lat = 90 − (y + ½)·180/H`.

**Patch or authored.** By default an import is diffed against the view as it was exported and
stored as a patch: deltas for continuous fields, new values for categories and owners. If the
layer is later regenerated underneath (a parameter changed, an upstream edit), the patch is
re-applied on top and `status` shows it as *rebased*; edits that no longer fit (a cell outside
the grid, a key with no row) are counted as *orphaned*. `--authored` replaces the product
outright: draw a whole map and the generator for that product is skipped. A new plate key
without a row gets a seed point at its centroid and a seeded motion; a new lake gets its
surface at the highest cell.

## Provenance

Each stored product records a hash of everything it was built from: the layer's version, the
root seed, the body, the grid, the layer's parameters, its products' sources and edit files,
and the provenance of every product the layer reads. `generate` reuses a product whose recorded
hash matches and rebuilds the rest, so an edit to elevation reruns elevation, climate,
hydrology and biomes, and nothing else. Bump a layer's version in
`geoforge-pipeline/src/layers.rs` when its algorithm changes output; until then use
`generate --from <layer>` while iterating on code.

## Seeds

A product's seed is `root_seed → product name → owner key`, never an ancestor, an index or a
name. Adding a layer, renaming a body or editing an unrelated product never changes another
product's randomness. `seed <dir> <product>` prints it.

## Models (short)

- **Tectonics**: v1's electrostatic plate seeding (power-law charges), now with each cell
  assigned by distance ÷ charge^0.22 so a seed's charge sets its plate's size
  (`size_weighting = 0` gives v1's plates exactly); domain warp and island removal;
  Euler-pole motion, per-cell boundary classification, crust age from ridge distance.
  Continents are embedded in plates with ocean floor to their passive margins; the largest
  plate and about 30% of the rest are oceanic.
- **Geology**: v1's 20 province types (collision belts now swell and pinch along their
  length), plus lithology: each province kind has a rock mix (craton: basement, greenstone,
  granite; platform: sandstone, limestone, salt; orogen: metamorphic core, granite,
  deformed sediment, ophiolite; arc: volcanics and plutons…) picked by smooth noise, with
  rock ages from the province's history.
- **Elevation**: GDH1 ocean depth from crust age, continental freeboard with plateaus,
  province relief, shelves and slopes planed to sea level, sea level from the ocean's volume,
  then stream-power erosion (Braun & Willett implicit scheme) with erodibility from
  lithology and uplift in active belts: rivers cut dendritic valleys while ranges stay high.
- **Circulation**: January and July winds (Hadley/Ferrel/polar cells scaled by rotation,
  seasonal ITCZ, monsoon, subtropical highs over each ocean basin); wind-driven ocean gyres
  (Sverdrup interior, Stommel western boundary currents, Ekman drift, equatorial currents,
  circumpolar current); SST anomaly from currents crossing the temperature gradient, carried
  downstream, plus coastal upwelling.
- **Climate**: zonal energy balance with ice–albedo feedback and a carbonate–silicate
  thermostat; seasons from obliquity and eccentricity; continentality and current-warmed air
  carried by the winds; water vapour carried by each season's winds, raining out where air
  rises or cools (storm tracks follow the seasons), recycled over land; global
  precipitation from a water budget; lapse rate and relief on the fine grid.
- **Hydrology**: Budyko (Fu) runoff; priority-flood drainage with jittered routing across
  flats; lakes from each depression's water balance (open where overflowing and kept open by
  rifting, depth or glaciation; endorheic and shrunk to balance where not); arid channel
  losses; Strahler-ordered rivers; watersheds by terminal outlet.
- **Biomes**: Köppen thresholds (10 °C tree line, 18 °C tropics, aridity threshold that
  depends on temperature and rain seasonality, dry-summer Mediterranean) and Whittaker's
  rainfall steps, alpine by altitude, wetlands on flat floodplains.

## Realism status

See `docs/V2_REALISM.md` for the numbers across seeds and the known gaps.
