"""Read a world's edit views (`geoforge export DIR all`) as numpy arrays."""
import json
from pathlib import Path

import numpy as np
import png
from PIL import Image


class World:
    def __init__(self, world_dir):
        self.edit = Path(world_dir) / "edit"
        grid = self.side("elevation")["grid"]
        self.W, self.H = grid["width"], grid["height"]
        self.lat = 90 - (np.arange(self.H) + 0.5) * 180 / self.H
        self.lon = -180 + (np.arange(self.W) + 0.5) * 360 / self.W
        # Relative cell area (cos latitude).
        self.area = np.repeat(np.cos(np.radians(self.lat))[:, None], self.W, axis=1)
        self._cache = {}

    def side(self, product):
        return json.loads((self.edit / f"{product}.edit.json").read_text())

    def field(self, product):
        """A continuous field, in its unit."""
        if product not in self._cache:
            enc = self.side(product)["encoding"]
            px = np.array(Image.open(self.edit / f"{product}.png")).astype(np.float64)
            self._cache[product] = enc["offset"] + enc["scale"] * px
        return self._cache[product]

    def vector(self, product):
        """(eastward, northward) components of a 16-bit RGB vector view."""
        enc = self.side(product)["encoding"]
        w, h, rows, info = png.Reader(filename=str(self.edit / f"{product}.png")).asDirect()
        px = np.vstack([np.asarray(r, dtype=np.float64) for r in rows]).reshape(h, w, info["planes"])
        return enc["offset"] + enc["scale"] * px[..., 0], enc["offset"] + enc["scale"] * px[..., 1]

    def partition(self, product):
        """(key per cell, {key: row}); 0 where no entity owns the cell."""
        if product not in self._cache:
            table = json.loads((self.edit / f"{product}.table.json").read_text())
            rgb = np.array(Image.open(self.edit / f"{product}.png").convert("RGB")).astype(np.int64)
            packed = (rgb[..., 0] << 16) | (rgb[..., 1] << 8) | rgb[..., 2]
            keys = np.zeros(packed.shape, np.int64)
            for t in table:
                keys[packed == int(t["color"][1:], 16)] = t["key"]
            self._cache[product] = (keys, {t["key"]: t["row"] for t in table})
        return self._cache[product]

    def rivers(self):
        return json.loads((self.edit / "rivers.json").read_text())

    def cell(self, lat, lon):
        y = int(np.clip((90 - lat) / 180 * self.H, 0, self.H - 1))
        x = int(((lon + 180) / 360 * self.W) % self.W)
        return y, x

    def nearest_land(self, lat, lon, radius=4):
        """The cell at (lat, lon), or the nearest land cell within `radius` cells."""
        y, x = self.cell(lat, lon)
        e = self.field("elevation")
        if e[y, x] >= 0:
            return y, x
        best = None
        for dy in range(-radius, radius + 1):
            for dx in range(-radius, radius + 1):
                yy, xx = y + dy, (x + dx) % self.W
                if 0 <= yy < self.H and e[yy, xx] >= 0 and (best is None or dy * dy + dx * dx < best[0]):
                    best = (dy * dy + dx * dx, yy, xx)
        return (best[1], best[2]) if best else (y, x)
