#!/usr/bin/env bash
# Earth A/B, layer by layer: each layer run on real inputs, its result beside the real one.
#   scripts/real-earth/ab.sh OUT_DIR [DATA_DIR] [OBS_DIR] [WIDTH]
# Worlds made in OUT_DIR:
#   chain   real topography, every later layer the model's (run.sh)
#   air     chain with real winds and SST departures: the climate on real air
#   real    air with real temperature and rain: rivers, biomes and Köppen on a real climate
#   plates  the model's mountains on Earth's real plates (known_earth_tectonics' run, 720×360)
# then a report of each (geoforge report) and the A/B page, OUT_DIR/page (ab_page.py).
# Needs the data of fetch.sh in DATA_DIR and NCAR GeoCAT's GPCP and ERA-Interim 2 m
# temperature files in OBS_DIR (see inject_obs.py).
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
OUT=$1
DATA=${2:-real-earth-data}
OBS=${3:-obs-data}
WIDTH=${4:-1024}
G=${GEOFORGE:-geoforge}
mkdir -p "$OUT"
CITIES=(--place "London:51.5,-0.1" --place "Lima:-12.05,-77.05" --place "Phoenix:33.45,-112.07" --place "Yakutsk:62.0,129.7" --place "Mumbai:19.1,72.9")

echo "chain: real topography"
"$HERE/run.sh" "$OUT/chain" "$DATA" "$WIDTH"

echo "air: real winds and SST departures"
rm -rf "${OUT:?}/air" && cp -r "$OUT/chain" "$OUT/air"
python3 -I "$HERE/inject_obs.py" "$OUT/air" "$DATA" "$OBS" air
for p in wind_jan wind_jul sst_anomaly; do "$G" import "$OUT/air" "$p" --authored >/dev/null; done
"$G" generate "$OUT/air" >/dev/null
"$G" export "$OUT/air" all >/dev/null
python3 -I "$HERE/check.py" "$OUT/air" "$DATA" "$OUT/air/earth_check.json"

echo "real: real temperature and rain"
rm -rf "${OUT:?}/real" && cp -r "$OUT/air" "$OUT/real"
python3 -I "$HERE/inject_obs.py" "$OUT/real" "$DATA" "$OBS" climate
for p in temperature temperature_range precipitation precipitation_warm_share; do "$G" import "$OUT/real" "$p" --authored >/dev/null; done
"$G" generate "$OUT/real" >/dev/null
"$G" export "$OUT/real" all >/dev/null
python3 -I "$HERE/check.py" "$OUT/real" "$DATA" "$OUT/real/earth_check.json"

echo "plates: mountains on Earth's plates"
(cd "$ROOT" && KNOWN_CASES_DUMP=1 cargo test -q --release -p geoforge-pipeline --test known_earth_tectonics >/dev/null)
rm -rf "${OUT:?}/plates"
"$G" new "$OUT/plates" --seed 1 --width 720 --archetype earth >/dev/null
"$G" generate "$OUT/plates" --until elevation >/dev/null
"$G" export "$OUT/plates" rough_elevation >/dev/null
python3 -I "$HERE/inject_dump.py" "$ROOT/target/known-cases/tectonics_elevation.f32" "$OUT/plates" rough_elevation
"$G" import "$OUT/plates" rough_elevation --authored >/dev/null
sed -i 's/^k = .*/k = 0.0/' "$OUT/plates/world.toml"
"$G" generate "$OUT/plates" >/dev/null

for w in chain air real; do "$G" report "$OUT/$w" "${CITIES[@]}" >/dev/null; done
"$G" report "$OUT/plates" >/dev/null
python3 -I "$HERE/ab_page.py" "$OUT" "$DATA" "$ROOT"
echo "page: $OUT/page/index.html"
