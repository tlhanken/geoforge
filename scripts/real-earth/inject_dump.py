"""Write a field dumped by the known-case tests (KNOWN_CASES_DUMP) into a world's edit view.

usage: python3 -I inject_dump.py FIELD.f32 WORLD_DIR PRODUCT
The dump is little-endian f32 in grid-cell order on an equirectangular grid of the world's size.
"""
import json
import sys
from pathlib import Path

import numpy as np
from PIL import Image


def main():
    dump, world, product = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
    side = json.loads((world / "edit" / f"{product}.edit.json").read_text())
    w, h = side["grid"]["width"], side["grid"]["height"]
    values = np.fromfile(dump, dtype="<f4").astype(np.float64)
    if values.size != w * h:
        sys.exit(f"{dump.name} holds {values.size} cells, the world {w}×{h}")
    enc = side["encoding"]
    px = np.clip(np.rint((values.reshape(h, w) - enc["offset"]) / enc["scale"]), 0, 65535).astype(np.uint16)
    Image.fromarray(px).save(world / "edit" / f"{product}.png")
    print(f"{dump.name} → {product}: {values.min():.0f}…{values.max():.0f}")


if __name__ == "__main__":
    main()
