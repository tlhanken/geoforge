#!/usr/bin/env bash
# Bundle the parts of three.js the viewer uses (entry.js) into one minified module.
# Run from the repository root; needs npm. Pins the three.js version below.
set -euo pipefail
version=0.187.0
work=$(mktemp -d)
trap 'rm -rf "${work:?}"' EXIT
cp scripts/vendor-three/entry.js "$work/"
(cd "$work" && npm init -y >/dev/null && npm install --silent "three@$version" esbuild@0.25 >/dev/null \
  && npx esbuild entry.js --bundle --minify --format=esm --legal-comments=none --outfile=three.js)
{ echo "// three.js $version (MIT License, Copyright © 2010-2026 three.js authors), the parts entry.js names."; cat "$work/three.js"; } > crates/libs/geoforge-view/web/three.js
echo "wrote crates/libs/geoforge-view/web/three.js ($(wc -c < crates/libs/geoforge-view/web/three.js) bytes)"
