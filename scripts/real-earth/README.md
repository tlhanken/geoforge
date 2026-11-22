# Real-Earth check

Runs the climate, ocean, hydrology and biome layers on Earth's actual topography and
compares each layer with the observed Earth. The world is an `--archetype earth` world. Its
`rough_elevation` is replaced by ETOPO1 through `import --authored`, erosion is switched off
(`k = 0`) so `elevation` is the same map, and everything downstream is generated as usual.

```sh
cargo build --profile iter -p geoforge-cli
GEOFORGE=target/iter/geoforge scripts/real-earth/run.sh /tmp/realearth
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
| Biomes | Shares of land; compare with `geoforge stats` Earth references |

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

## Layer by layer: `ab.sh`

`ab.sh OUT_DIR [DATA_DIR] [OBS_DIR] [WIDTH]` runs each layer on Earth's real inputs and
sets it beside Earth's real result, so a miss belongs to that layer alone:

| World | Real inputs | Model layers compared |
|---|---|---|
| `plates` | PB2002 plates and Euler poles (the `known_earth_tectonics` run, 720 × 360) | geology and elevation, against ETOPO1 |
| `chain` | ETOPO1 (`run.sh`) | winds and currents, against ERA-Interim and ERSST |
| `air` | ETOPO1, ERA-Interim 850 hPa winds, ERSST departures (`inject_obs.py air`) | climate, against ERA-Interim 2 m temperature and GPCP |
| `real` | all of the above plus real temperature and rain (`inject_obs.py climate`) | the Köppen classifier, rivers, biomes |

It writes a `geoforge report` of each world and the A/B page, `OUT_DIR/page`
(`ab_page.py`), whose maps are geoforge's renders of both the real and the model field. The
extra observations come from NCAR GeoCAT's data files (raw.githubusercontent.com/NCAR/
geocat-datafiles/main/netcdf_files/): `V22_GPCP.1979-2010.nc` and
`T2M_ERAINT_rectilinear_grid_2D.nc`, in OBS_DIR.

