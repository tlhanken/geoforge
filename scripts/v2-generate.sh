#!/usr/bin/env bash
# Generate pre-tectonic pipeline output for review (JSON + text under outputs/).
#
# Usage:
#   ./scripts/v2-generate.sh              # seed 42 → ./outputs/
#   ./scripts/v2-generate.sh 99           # seed 99 → ./outputs/
#   ./scripts/v2-generate.sh 99 ./inspect both
#   ./scripts/v2-generate.sh 42 ./outputs text -v
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

SEED="${1:-42}"
OUT_DIR="${2:-./outputs}"
FMT="${3:-both}"
EXTRA=()
if [[ "${4:-}" == "-v" ]] || [[ "${4:-}" == "--verbose" ]]; then
  EXTRA=(-v)
fi

if ! command -v cargo >/dev/null 2>&1; then
  echo "error: cargo not found in PATH" >&2
  exit 1
fi

echo "Building geoforge_v2..."
cargo build -p geoforge_v2 -q

mkdir -p "$OUT_DIR"

# Flags after `--` are passed to the binary, not Cargo.
cargo run -p geoforge_v2 -q -- generate \
  --seed "$SEED" \
  --format "$FMT" \
  --output-directory "$OUT_DIR" \
  "${EXTRA[@]}"

echo ""
echo "Review files:"
ls -la "$OUT_DIR"/cosmology_seed"${SEED}".* 2>/dev/null || ls -la "$OUT_DIR"/
