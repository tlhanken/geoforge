#!/usr/bin/env bash
# Download the real-Earth reference data (≈30 MB) into DATA_DIR (default: ./real-earth-data).
#   ETOPO1 topography, 10′          Fatiando a Terra mirror of NOAA ETOPO1 (public domain)
#   ERA-Interim winds, Jan & Jul    xarray tutorial data (ECMWF ERA-Interim, monthly means)
#   ERSST v5 sea-surface temps      xarray tutorial data (NOAA, public domain)
#   Köppen–Geiger map (present)     kgcpy package (Rubel et al. 2017; BSD licence)
set -euo pipefail
DATA=${1:-real-earth-data}
mkdir -p "$DATA"
cd "$DATA"
get() { [ -s "$2" ] || curl -sSfL -o "$2" "$1"; echo "  $2"; }
get https://github.com/fatiando-data/earth-topography-10arcmin/releases/download/v1/earth-topography-10arcmin.nc earth-topography-10arcmin.nc
get https://raw.githubusercontent.com/pydata/xarray-data/master/eraint_uvz.nc eraint_uvz.nc
get https://raw.githubusercontent.com/pydata/xarray-data/master/ersstv5.nc ersstv5.nc
if [ ! -s kgcpy/kmz_int_reshape.png ]; then
  python3 -m pip download --no-deps --quiet -d wheel kgcpy==1.1.8
  python3 -m zipfile -e wheel/kgcpy-1.1.8-py3-none-any.whl wheel/x
  mkdir -p kgcpy && cp wheel/x/kgcpy/kmz_int_reshape.png wheel/x/kgcpy/kg_zoneNum.csv kgcpy/
  rm -rf wheel
fi
echo "  kgcpy/kmz_int_reshape.png"
