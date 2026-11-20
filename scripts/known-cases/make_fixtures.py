"""Build the known-case test fixtures (crates/libs/geoforge-pipeline/tests/fixtures).

usage: python3 -I make_fixtures.py DATA_DIR OUT_DIR

DATA_DIR holds earth-topography-10arcmin.nc (scripts/real-earth/fetch.sh) and PB2002's
GeoJSON plates and steps plus PB2002_poles.dat.txt (github.com/fraxen/tectonicplates).
Writes, on equirectangular grids with pixel-centred cells:
  earth_etopo1_360x180.png     ETOPO1 block means, 16-bit, metres + 12000
  pb2002_plates_720x360.png    plate index per half-degree cell, 8-bit (row order of
                               pb2002_plates.csv)
  pb2002_plates.csv            code, name, Euler pole relative to the Pacific (PB2002)
  pb2002_steps.csv             boundary steps: midpoint, class, plate pair
"""
import json
import sys
from pathlib import Path

import numpy as np
import xarray as xr
from matplotlib.path import Path as MPath
from PIL import Image
from scipy import ndimage

W, H = 360, 180
LAT = 90 - (np.arange(H) + 0.5) * 180 / H
LON = -180 + (np.arange(W) + 0.5) * 360 / W


def topography(nc):
    """Block means of ETOPO1 over one-degree cells (from a 6×6 sample of each)."""
    k = 6
    lat = 90 - (np.arange(H * k) + 0.5) * 180 / (H * k)
    lon = -180 + (np.arange(W * k) + 0.5) * 360 / (W * k)
    t = xr.load_dataset(nc).topography.interp(latitude=xr.DataArray(lat, dims="y"), longitude=xr.DataArray(lon, dims="x")).values
    e = t.reshape(H, k, W, k).mean(axis=(1, 3))
    e = burn_rivers(t, e, k)
    # Below-sea land the ocean does not reach is raised to +1 m, unless it is a sea in its
    # own right (≥ 100 000 km²: the Caspian, and any sea the coarse grid cut off).
    lab, _ = ndimage.label(e < 0)
    for y in range(H):
        a, b = lab[y, 0], lab[y, -1]
        if a and b and a != b:
            lab[lab == b] = a
    area = np.repeat((np.cos(np.radians(LAT)) * (np.pi / 180 * 6371) ** 2)[:, None], W, axis=1)
    sizes = ndimage.sum(area, lab, index=np.arange(lab.max() + 1))
    ocean = np.argmax(sizes[1:]) + 1
    small = (lab > 0) & (lab != ocean) & (sizes[lab] < 1e5)
    e[small] = 1.0
    return e


def burn_rivers(t, e, k, min_basin_km2=500_000, depression_m=5.0):
    """Lower each block a big river crosses to the river's level.

    A block mean dams a valley narrower than the block: the Nile's floor below Aswan is
    176 m, but its one-degree blocks average 230–260 m. The river is traced on the fine
    grid (priority-flood from the ocean: each cell drains to the cell that flooded it), and
    a block a river crosses takes the lowest river level in it where that is below its
    mean. Only river cells outside filled depressions count, so real closed basins (the
    Great Basin, the Tarim) keep their rims.
    """
    import heapq

    h, w = t.shape
    lab, _ = ndimage.label(t < 0)
    lat = 90 - (np.arange(h) + 0.5) * 180 / h
    cell_km2 = np.repeat((np.cos(np.radians(lat)) * (np.pi / 180 * 6371 / (h / 180)) ** 2)[:, None], w, axis=1)
    sizes = ndimage.sum(cell_km2, lab, index=np.arange(lab.max() + 1))
    sea = (lab > 0) & (sizes[lab] >= 1e5)
    filled = np.where(sea, t, np.inf)
    to = np.full(t.size, -1, np.int64)
    heap = [(t[y, x], y * w + x) for y, x in zip(*np.nonzero(sea))]
    heapq.heapify(heap)
    order = []
    tf, ff = t.ravel(), filled.ravel()
    while heap:
        z, i = heapq.heappop(heap)
        order.append(i)
        y, x = divmod(i, w)
        for dy in (-1, 0, 1):
            ny = y + dy
            if not 0 <= ny < h:
                continue
            for dx in (-1, 0, 1):
                j = ny * w + (x + dx) % w
                if ff[j] == np.inf:
                    ff[j] = max(tf[j], z)
                    to[j] = i
                    heapq.heappush(heap, (ff[j], j))
    acc = cell_km2.ravel().copy()
    for i in reversed(order):
        if to[i] >= 0:
            acc[to[i]] += acc[i]
    river = (acc >= min_basin_km2) & (ff - tf < depression_m) & (tf >= 0)
    level = np.where(river, ff, np.inf).reshape(h, w)
    level = level.reshape(h // k, k, w // k, k).min(axis=(1, 3))
    burned = level < e
    print(f"rivers: {burned.sum()} blocks lowered, by {np.median((e - level)[burned]):.0f} m (median)")
    return np.where(burned, level, e)


def plates(geojson, poles_txt, w=720, h=360):
    feats = json.loads(Path(geojson).read_text())["features"]
    codes = []
    for f in feats:
        c = f["properties"]["Code"]
        if c not in codes:
            codes.append(c)
    names = {f["properties"]["Code"]: f["properties"]["PlateName"] for f in feats}
    poles = {}
    for line in Path(poles_txt).read_text().splitlines():
        parts = line.split()
        if len(parts) >= 4 and len(parts[0]) == 2 and parts[0].isupper():
            poles[parts[0]] = (float(parts[1]), float(parts[2]), float(parts[3]))
    lat_c = 90 - (np.arange(h) + 0.5) * 180 / h
    lon_c = -180 + (np.arange(w) + 0.5) * 360 / w
    lon, lat = np.meshgrid(lon_c, lat_c)
    pts = np.column_stack([lon.ravel(), lat.ravel()])
    idx = np.full(w * h, 255, np.uint8)
    for f in feats:
        g = f["geometry"]
        polys = [g["coordinates"]] if g["type"] == "Polygon" else g["coordinates"]
        i = codes.index(f["properties"]["Code"])
        for poly in polys:
            ring = np.array(poly[0])
            inside = MPath(ring).contains_points(pts, radius=1e-9)
            idx[inside & (idx == 255)] = i
    idx = idx.reshape(h, w)
    gaps = idx == 255
    if gaps.any():  # cells on a polygon edge: take the nearest plate
        _, (iy, ix) = ndimage.distance_transform_edt(gaps, return_indices=True)
        idx = idx[iy, ix]
    rows = [(c, names[c], *poles[c]) for c in codes]
    return idx, rows, int(gaps.sum())


def steps(geojson):
    out = []
    for f in json.loads(Path(geojson).read_text())["features"]:
        p = f["properties"]
        lat = (p["STARTLAT"] + p["FINALLAT"]) / 2
        lon = (p["STARTLONG"] + p["FINALLONG"]) / 2
        if abs(p["STARTLONG"] - p["FINALLONG"]) > 180:
            lon = ((lon + 360) % 360) - 180
        out.append((round(lat, 3), round(lon, 3), p["STEPCLASS"], p["PLATEBOUND"], p["STEPLENGTH"]))
    return out


def main():
    data, out = Path(sys.argv[1]), Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    e = topography(data / "earth-topography-10arcmin.nc")
    Image.fromarray(np.clip(np.rint(e + 12000), 0, 65535).astype(np.uint16)).save(out / "earth_etopo1_360x180.png")
    idx, rows, gaps = plates(data / "PB2002_plates.json", data / "PB2002_poles.dat.txt")
    Image.fromarray(idx).save(out / "pb2002_plates_720x360.png")
    with open(out / "pb2002_plates.csv", "w") as f:
        f.write("code,name,pole_lat,pole_lon,deg_per_myr\n")
        for r in rows:
            f.write(f"{r[0]},{r[1]},{r[2]},{r[3]},{r[4]}\n")
    st = steps(data / "PB2002_steps.json")
    with open(out / "pb2002_steps.csv", "w") as f:
        f.write("lat,lon,class,plates,length_km\n")
        for s in st:
            f.write(f"{s[0]},{s[1]},{s[2]},{s[3]},{s[4]}\n")
    land = (e >= 0).mean()
    print(f"topography: land {land:.3f} of cells; plates: {len(rows)}, {gaps} edge cells filled; steps: {len(st)}")


if __name__ == "__main__":
    main()
