#!/usr/bin/env bash
# Generate a v2 world with preview renders and realism stats.
#
# Usage:
#   ./scripts/v2-generate.sh                 # seed 42 → ./outputs/world-42
#   ./scripts/v2-generate.sh 99              # seed 99 → ./outputs/world-99
#   ./scripts/v2-generate.sh 99 ./worlds     # seed 99 → ./worlds/world-99
#   EARTH=1 ./scripts/v2-generate.sh 42      # Earth-like reference body
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

SEED="${1:-42}"
OUT_DIR="${2:-./outputs}"
WORLD="$OUT_DIR/world-$SEED"

if ! command -v cargo >/dev/null 2>&1; then
  echo "error: cargo not found in PATH" >&2
  exit 1
fi

echo "Building geoforge (release)..."
cargo build --release -p geoforge-cli -q
BIN="$ROOT/target/release/geoforge"

mkdir -p "$OUT_DIR"
if [[ ! -f "$WORLD/world.toml" ]]; then
  if [[ "${EARTH:-}" == "1" ]]; then
    "$BIN" new "$WORLD" --seed "$SEED" --earth-like
  else
    "$BIN" new "$WORLD" --seed "$SEED"
  fi
fi
"$BIN" generate "$WORLD" --render
"$BIN" stats "$WORLD" > "$WORLD/stats.json"
echo
echo "World:   $WORLD"
echo "Renders: $WORLD/renders"
echo "Stats:   $WORLD/stats.json"
echo "Edit:    $BIN export $WORLD elevation   (then paint edit/elevation.png and: $BIN import $WORLD elevation)"
