# Reference data for the Earth checks

The observed Earth each layer is checked against: what the A/B page (`geoforge compare`,
TODO item 1) and the known-case suites need, where it comes from, and whether this
project's cloud sessions can fetch it (their network reaches GitHub and the package
registries only; anything else is downloaded by hand into `real-earth-data/`, which git
ignores).

Python converts each source once into world layers on the Earth grid; everything after that
is Rust.

## In use

| Layer | Source | File in `real-earth-data/` |
|---|---|---|
| Elevation | ETOPO1 / ETOPO 2022, 10′ | `earth-topography-10arcmin.nc` |
| Winds, pressure | ERA-Interim 850 hPa monthly means | `eraint_uvz.nc` |
| Sea surface temperature | NOAA ERSST v5 | `ersstv5.nc` |
| Köppen–Geiger | Beck et al. 2018 (via `kgcpy`) | `kgcpy/` |
| Plates | Bird 2003, PB2002 | in the tectonics fixtures |
| Rain | GPCP v2.3 (station normals in the climate suite) | constants in the tests |

## Missing baselines

| Layer | Source | Licence | Format, size | Reachable here |
|---|---|---|---|---|
| Biomes | RESOLVE Ecoregions 2017 (Dinerstein et al. 2017): 846 ecoregions in 14 biomes; the update of Olson et al. 2001 (WWF) | CC BY 4.0 | shapefile, ~150 MB zip | no: download `Ecoregions2017.zip` from ecoregions.appspot.com, or export from Earth Engine (`RESOLVE/ECOREGIONS/2017`) |
| Ocean currents | OSCAR v2.0, 0.25°, daily, 1993–2022 (ESR, PO.DAAC); its older 1/3° 5-day product has yearly files | NASA open data, cite the DOI | netCDF; a climatology needs one year of files at least | no: NASA Earthdata login |
| Sea ice | NSIDC Sea Ice Index G02135 v4: monthly extent and the 1981–2010 median edge as shapefiles; monthly totals as spreadsheets | free, cite Fetterer et al. | shapefile, GeoTIFF, xlsx | no: `noaadata.apps.nsidc.org` / `masie_web.apps.nsidc.org` |
| Tropical cyclones | IBTrACS v04r01 (NOAA NCEI): every track since 1842, 3-hourly | public domain (NOAA) | CSV, ~65 MB compressed for all; `ibtracs.since1980.list.v04r01.csv` is smaller | no: ncei.noaa.gov |
| Earthquakes | USGS ComCat, M ≥ 5 since 1973 | public domain | CSV/GeoJSON from the FDSN event API | no: earthquake.usgs.gov |
| Volcanoes | Smithsonian GVP, Volcanoes of the World v5.4 (~1,340 Holocene volcanoes) | free, cite the DOI | Excel XML from the Holocene list page | no: volcano.si.edu |
| Rivers drawn | HydroRIVERS v1.0 (Lehner & Grill 2013), every river over 10 km² or 0.1 m³/s, with discharge | free incl. commercial, HydroSHEDS licence | shapefile/geodatabase, by continent | no: hydrosheds.org |
| Rivers drawn (light) | Natural Earth rivers and lake centrelines, 1:50 m and 1:10 m, ranked | public domain | GeoJSON in `nvkelso/natural-earth-vector` | yes (GitHub) |
| Stars | HYG v4.1 (Hipparcos, Yale, Gliese), ~120,000 stars; v4.2 and AT-HYG now live on Codeberg | CC BY-SA 4.0 | CSV in `astronexus/HYG-Database` | yes (GitHub) |
| Planets' positions | JPL DE440 ephemeris | public domain | SPICE kernel, ~115 MB | no: ssd.jpl.nasa.gov |

## How each is used

- **Biomes:** rasterise `BIOME_NUM` onto the Earth grid and map the 14 WWF biomes onto
  ours (tundra, boreal forest, temperate broadleaf and so on; mangroves and flooded
  grasslands to wetland). The A/B page's biome section then shows the real map beside the
  model fed real climate, and the climate suite's 20 biome places can read the map instead
  of hand-written expectations.
- **Currents:** the annual mean of OSCAR's u and v at 1°, for direction agreement and the
  boundary-current boxes, with the real arrows on the page.
- **Sea ice:** March and September median edges for the sea ice map; the suite already
  checks the totals.
- **Hazards:** cyclone track density per 10⁴ km² a year, earthquake rate by M ≥ 5 per cell,
  Holocene volcano points; each drawn beside the model's.
- **Rivers:** Natural Earth for the drawn network now; HydroRIVERS where discharge along the
  channel is needed.
- **Sky:** HYG for the star chart's real sky (Earth's own run only).
