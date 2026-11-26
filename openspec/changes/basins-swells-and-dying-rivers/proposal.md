# basins-swells-and-dying-rivers

## Why

On seven Earth-like worlds the five largest drainage basins took 28% of the land (Earth 14%),
the largest was 2–6× the Amazon's share, peaks stood 10–11.8 km high at Earth's gravity, and
only 7% of land drained to no sea (Earth 18%). Run on Earth's own one-degree map, the same
hydrology gives Earth's basins (Amazon 6.2 Mkm², Congo 3.8, Nile 3.6), so the drainage code was
sound and the terrain was not: stable interiors were flat, ranges rose without limit, and desert
slopes drained to the sea however little water they carried.

## What Changes

- **Basins and swells** (`ElevationParams::basin_swell_m`, `basin_swell_km`; 1500 m and
  1200 km): old stable crust carries broad sags and rises, as Africa's basins lie between its
  swells, subdued behind active margins and absent in forelands and rifts; the dynamic lift far from margins (`plateau_m`) falls from 1800 to 1200 m so the
  land's height stays Earth's.
- **Peaks held** (`ErosionParams::peak_limit_m`, 9000 m at Earth's gravity, × 1/g;
  `erosion::hold_peak`): heights above three-quarters of the limit bend toward it.
- **Dying desert rivers** (`HydrologyParams::dead_runoff_mm`, 10 mm/yr): a basin of at least
  `min_basin_km2` whose land yields less on average has no river to the sea.
- **Earth's reference corrected**: rivers over 10 000 m³/s to the sea are 16, not 8; the
  reference for over 1000 m³/s, unsourced, is dropped.

## Impact

- `geoforge-surface`: `elevation`, `erosion`, `hydrology` (new params with serde defaults:
  worlds without them in `world.toml` take the new defaults when regenerated).
- `geoforge-pipeline`: `stats`; new test `earthlike_terrain`.
- Every Earth-like world's elevation, rivers, lakes and basins change, and with them its
  climate a little (land 0.276 of the surface on average, was 0.29).
