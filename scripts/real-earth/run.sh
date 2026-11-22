#!/usr/bin/env bash
# Build a world on real Earth topography and check every generated layer against Earth.
#   scripts/real-earth/run.sh WORLD_DIR [DATA_DIR] [WIDTH]
# Needs: geoforge on PATH (or GEOFORGE=path), python3 with numpy, scipy, xarray, netCDF4,
# pillow, pypng; the data from fetch.sh.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
WORLD=$1
DATA=${2:-real-earth-data}
WIDTH=${3:-2048}
G=${GEOFORGE:-geoforge}
[ -s "$DATA/earth-topography-10arcmin.nc" ] || "$HERE/fetch.sh" "$DATA"
rm -rf "$WORLD"
"$G" new "$WORLD" --seed 1 --width "$WIDTH" --archetype earth >/dev/null
"$G" generate "$WORLD" --until elevation >/dev/null
"$G" export "$WORLD" elevation >/dev/null
python3 -I "$HERE/inject.py" "$DATA/earth-topography-10arcmin.nc" "$WORLD"
"$G" import "$WORLD" elevation --authored >/dev/null
"$G" generate "$WORLD" --render >/dev/null
"$G" export "$WORLD" all >/dev/null
"$G" stats "$WORLD" > "$WORLD/stats.json"
python3 -I "$HERE/check.py" "$WORLD" "$DATA" "$WORLD/earth_check.json"
