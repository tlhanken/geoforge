#!/usr/bin/env bash
# Diagnose why geoforge commands fail locally.
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

echo "Branch: $(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo '?')"
if [[ ! -d "$ROOT/crates/bins/geoforge-cli" ]]; then
  echo "ERROR: crates/bins/geoforge-cli missing — check out the right branch."
  exit 1
fi
echo "OK: crates/bins/geoforge-cli exists"

if ! command -v cargo >/dev/null 2>&1; then
  echo "ERROR: cargo not in PATH (nix develop, or install rustup)."
  exit 1
fi
echo "cargo: $(cargo --version)"

echo "Building geoforge-cli..."
if ! cargo build --release -p geoforge-cli -q; then
  echo "ERROR: build failed (see above)."
  exit 1
fi
BIN="$ROOT/target/release/geoforge"
"$BIN" --version || exit 1
echo
echo "Try:"
echo "  $BIN new /tmp/demo --seed 42 && $BIN generate /tmp/demo --render"
