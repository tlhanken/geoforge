#!/usr/bin/env bash
# Earth A/B, layer by layer: each layer run on real inputs, its result beside the real one.
#   scripts/real-earth/ab.sh OUT_DIR [DATA_DIR] [OBS_DIR] [WIDTH]
# Worlds made in OUT_DIR:
#   chain   real topography, every later layer the model's (run.sh)
#   air     chain with real winds and SST departures: the climate on real air
#   real    air with real temperature and rain: rivers, biomes and Köppen on a real climate
#   plates  the model's mountains on Earth's real plates (known_earth_tectonics' run, 720×360)
# then the A/B page, OUT_DIR/page (geoforge compare): every product the real world holds as
# observed data beside the same product in each model world.
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
SOURCES=(
  --source "rough_elevation=ETOPO1, 10′" --source "wind_jan=ERA-Interim, 850 hPa (above the friction layer)" --source "wind_jul=ERA-Interim, 850 hPa (above the friction layer)"
  --source "sst_anomaly=ERSST v5, departure from each latitude's mean" --source "temperature=ERA-Interim 2 m, 2005"
  --source "temperature_range=ERA-Interim 2 m, 2005" --source "precipitation=GPCP v2.2, 1979–2010"
  --source "koppen=read by geoforge from the observed temperature and rain"
)
# Stations, the mountain ranges and the boundary-current boxes the page lists.
PROBES=(
  --probe "London:51.5,-0.1" --probe "Moscow:55.75,37.6" --probe "Yakutsk:62.0,129.7" --probe "Mumbai:19.1,72.9"
  --probe "Lima:-12.05,-77.05" --probe "Phoenix:33.45,-112.07" --probe "Tehran:35.7,51.4" --probe "Riyadh:24.7,46.7"
  --probe "Perth:-31.95,115.85" --probe "Astana:51.15,71.45" --probe "Manaus:-3.1,-60.0" --probe "Nairobi:-1.3,36.8"
  --probe "Rockies:37..50,-115..-104" --probe "Colorado Plateau:35..39,-112..-107" --probe "Appalachians:35..45,-84..-74"
  --probe "Altiplano:-22..-14,-70..-65" --probe "Himalaya and Tibet:28..36,80..95" --probe "Alps:44..48,6..15"
  --probe "Ethiopian Highlands:7..14,36..40" --probe "East African Plateau:-5..2,30..37" --probe "Urals:50..65,57..61"
  --probe "Brazilian Highlands:-23..-15,-48..-42" --probe "Iranian Plateau:29..36,52..60" --probe "Great Dividing Range:-37..-20,145..152"
  --probe "NE Atlantic:45..65,-30..10" --probe "Kuroshio:30..40,130..150" --probe "California Current:25..40,-130..-115"
  --probe "Humboldt:-25..-5,-85..-72" --probe "Benguela:-30..-15,5..14" --probe "Canary:15..30,-25..-12" --probe "W Pacific warm pool:-10..10,120..160"
)

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
"$G" generate "$OUT/plates" --until elevation >/dev/null

"$G" compare "$OUT/real" "$OUT/chain" "$OUT/air" "$OUT/plates" --label "model air and sea" --label "real air and sea" --label "real plates" "${PROBES[@]}" "${SOURCES[@]}" --to "$OUT/page"
echo "page: $OUT/page/index.html"
