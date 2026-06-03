#!/usr/bin/env bash
# Diagnose why geoforge_v2 commands fail locally.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

echo "=== Geoforge v2 doctor ==="
echo "Repo: $ROOT"
echo ""

# Git branch
BRANCH="$(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo 'unknown')"
echo "Git branch: $BRANCH"
if [[ "$BRANCH" != "cursor/v2-types-foundation-b88c" ]]; then
  echo "  WARNING: cosmology export lives on branch cursor/v2-types-foundation-b88c"
  echo "  main has NO geoforge_v2 package. cleanup_and_nix is an older stub."
  echo "  Fix: git fetch origin && git checkout cursor/v2-types-foundation-b88c"
fi
echo ""

# Workspace
if [[ ! -d "$ROOT/crates/bins/geoforge_v2" ]]; then
  echo "ERROR: crates/bins/geoforge_v2 missing — you are on the wrong branch."
  echo "  Run: git fetch origin && git checkout cursor/v2-types-foundation-b88c"
  exit 1
fi
echo "OK: crates/bins/geoforge_v2 exists"
echo ""

# Rust
if ! command -v cargo >/dev/null 2>&1; then
  echo "ERROR: cargo not in PATH"
  exit 1
fi
echo "cargo: $(command -v cargo)"
cargo --version
rustc --version
echo ""

# Build
echo "Building geoforge_v2..."
if ! cargo build -p geoforge_v2 -q 2>/dev/null; then
  echo "ERROR: build failed. If you see edition2024, update Rust:"
  echo "  rustup update stable"
  exit 1
fi
echo "OK: build succeeded"
BIN="$ROOT/target/debug/geoforge_v2"
echo "Binary: $BIN"
echo ""

# CLI flags
echo "Export-related flags in --help:"
if ! "$BIN" generate --help | grep -q 'format'; then
  echo "ERROR: --format not in help — binary is too old. git pull && cargo build -p geoforge_v2"
  exit 1
fi
"$BIN" generate --help | grep -E 'format|export|no-export' || true
echo ""

# Smoke run
OUT="$ROOT/target/doctor-out"
rm -rf "$OUT"
mkdir -p "$OUT"
"$BIN" generate --seed 42 --format text -d "$OUT" | tail -3
if [[ -f "$OUT/cosmology_seed42.txt" ]]; then
  echo "OK: wrote $OUT/cosmology_seed42.txt ($(wc -c < "$OUT/cosmology_seed42.txt") bytes)"
else
  echo "ERROR: export file not created"
  exit 1
fi

echo ""
echo "All checks passed. Use:"
echo "  cargo run -p geoforge_v2 -- generate --seed 42 --format text -d ./inspect"
echo "  (note the -- before generate)"
