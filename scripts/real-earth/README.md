# Real-Earth check

Runs the climate, ocean, hydrology and biome layers on Earth's actual topography and
compares each layer with the observed Earth. The world is an `--archetype earth` world. Its
elevation is replaced by ETOPO1 through `import --authored`, and everything downstream of
elevation is generated as usual.

```sh
cargo build --profile iter -p geoforge-cli
GEOFORGE=target/iter/geoforge-v2 scripts/real-earth/run.sh /tmp/realearth
```

The first run downloads about 30 MB into `./real-earth-data` (see `fetch.sh`). A 2048×1024
world takes about a minute. The script prints a summary and writes the full report to
`WORLD/earth_check.json`, the world's stats to `WORLD/stats.json`, and renders to
`WORLD/renders/`.

Python needs `numpy`, `scipy`, `xarray` with `netCDF4`, `pillow` and `pypng`.

## What is compared

| Layer | Reference |
|---|---|
| Winds, January and July | ERA-Interim 850 hPa monthly means: direction agreement over the ocean, zonal bands, the northern Indian Ocean monsoon |
| Sea-surface temperature | ERSST v5, 1971–2000: bias, RMSE, correlation of the departure from the zonal mean, warm and cold current regions |
| Land climate | Köppen–Geiger map (Rubel et al., from the `kgcpy` package), with Köppen classes derived from the model's mean, range, precipitation and warm-season share; 31 station normals |
| Ocean rainfall | Zonal means, roughly as in GPCP |
| Rivers | Discharge near 18 river mouths |
| Biomes | Shares of land; compare with `geoforge-v2 stats` Earth references |

Station normals, river discharges and the ocean rainfall table are approximate, rounded
reference values. They are good to about ±1 °C and ±15%.

## Caveats

- Cells are about 20 km across, so a station is compared with its cell. That cell can be
  higher than the station (Tehran's cell is 600 m higher, Lhasa's 600 m) or include
  mountains beside it.
- ETOPO1 is not hydrologically conditioned. At 20 km, narrow gorges are blocked, so the
  Congo ponds behind its rapids, and the upper Mississippi and middle Yangtze end in closed
  basins. River figures are indicative only.
- Below-sea land the ocean does not reach is raised to +1 m. The Black Sea and Caspian
  stay water.
- Tectonics and provinces are still generated, so they do not match Earth. Only
  hydrology's rift lakes read them.
