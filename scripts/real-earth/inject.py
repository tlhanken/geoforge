"""Write real Earth topography into a world's elevation edit view.

usage: python3 -I inject.py ETOPO_NC WORLD_DIR
(then: geoforge import WORLD_DIR elevation --authored)

ETOPO1 is sampled at the world's pixel centres. Land below sea level that the ocean does not
reach is raised to +1 m (the Dead Sea, Qattara, Turpan, the Great Lakes' floors and Lake
Eyre become dry depressions that hydrology may fill as lakes), except the Black Sea and the
Caspian, which stay water.
"""
import json
import sys
from pathlib import Path

import numpy as np
import xarray as xr
from PIL import Image
from scipy import ndimage

KEEP_WATER = {"Pacific": (0.0, -150.0), "Black Sea": (43.5, 34.5), "Caspian": (42.0, 50.5)}


def resample(nc, width, height):
    lat = 90 - (np.arange(height) + 0.5) * 180 / height
    lon = -180 + (np.arange(width) + 0.5) * 360 / width
    topo = xr.load_dataset(nc).topography
    return topo.interp(latitude=xr.DataArray(lat, dims="y"), longitude=xr.DataArray(lon, dims="x")).values.astype(np.float64)


def fill_enclosed(e):
    h, w = e.shape
    lab, _ = ndimage.label(e < 0)
    for y in range(h):  # one basin across the dateline
        a, b = lab[y, 0], lab[y, -1]
        if a and b and a != b:
            lab[lab == b] = a
    keep = {lab[int((90 - la) / 180 * h), int((lo + 180) / 360 * w)] for la, lo in KEEP_WATER.values()} - {0}
    enclosed = (lab > 0) & ~np.isin(lab, list(keep))
    e[enclosed] = 1.0
    return int(enclosed.sum())


def main():
    nc, world = sys.argv[1], Path(sys.argv[2])
    product = sys.argv[3] if len(sys.argv) > 3 else "rough_elevation"
    side = json.loads((world / "edit" / f"{product}.edit.json").read_text())
    w, h = side["grid"]["width"], side["grid"]["height"]
    if side["grid"]["kind"] != "equirect":
        sys.exit("inject.py needs an equirectangular world")
    e = resample(nc, w, h)
    raised = fill_enclosed(e)
    enc = side["encoding"]
    px = np.clip(np.rint((e - enc["offset"]) / enc["scale"]), 0, 65535).astype(np.uint16)
    Image.fromarray(px).save(world / "edit" / f"{product}.png")
    print(f"ETOPO1 at {w}×{h}: {raised} enclosed below-sea cells raised to +1 m; range {e.min():.0f}…{e.max():.0f} m")


if __name__ == "__main__":
    main()
