#!/usr/bin/env bash
# Verified wrapper — avoids Cargo swallowing --export before `--`.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

SEED="${1:-42}"
OUT_DIR="${2:-./inspect}"
FMT="${3:-text}"

if ! command -v cargo >/dev/null 2>&1; then
  echo "error: cargo not found in PATH" >&2
  exit 1
fi

cargo build -p geoforge_v2 -q

# IMPORTANT: everything after `--` goes to geoforge_v2, not to cargo.
cargo run -p geoforge_v2 -- generate \
  --seed "$SEED" \
  --format "$FMT" \
  --output-directory "$OUT_DIR"

echo ""
echo "Done. Inspect:"
ls -la "$OUT_DIR"/cosmology_seed"${SEED}".* 2>/dev/null || ls -la "$OUT_DIR"/
