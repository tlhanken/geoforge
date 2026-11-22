"""Write observed Earth fields into a world's edit views, so a layer can run on real inputs.

usage: python3 -I inject_obs.py WORLD_DIR DATA_DIR OBS_DIR air|climate
(after: geoforge export WORLD_DIR <products>; then: geoforge import WORLD_DIR <product> --authored)

air      wind_jan, wind_jul    ERA-Interim 850 hPa winds, January and July (DATA_DIR/eraint_uvz.nc)
         sst_anomaly           ERSST v5 1971–2000 mean, departure from each latitude's ocean mean
climate  temperature           ERA-Interim 2 m temperature, 2005 annual mean (OBS_DIR/T2M_ERAINT_rectilinear_grid_2D.nc)
         temperature_range     its warmest minus coldest month
         precipitation         GPCP v2.2, 1979–2010 annual mean (OBS_DIR/V22_GPCP.1979-2010.nc)
         precipitation_warm_share  GPCP's share in the summer half-year (April–September north,
                               October–March south)

OBS_DIR files come from NCAR's GeoCAT data files (github.com/NCAR/geocat-datafiles).
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import png
import xarray as xr
from PIL import Image

sys.path.insert(0, str(Path(__file__).parent))
from gfview import World  # noqa: E402


def to_grid(da, lat_name, lon_name, w, fill=True):
    """Interpolate a (lat, lon) array to the world's pixel centres, longitudes wrapped; with
    `fill`, gaps take the nearest value in their column."""
    da = da.assign_coords({lon_name: ((da[lon_name] + 180) % 360) - 180}).sortby(lon_name).sortby(lat_name)
    # Pad one column each side so interpolation wraps the dateline.
    lon = da[lon_name].values
    da = xr.concat([da.isel({lon_name: [-1]}).assign_coords({lon_name: [lon[-1] - 360]}), da, da.isel({lon_name: [0]}).assign_coords({lon_name: [lon[0] + 360]})], dim=lon_name)
    out = da.interp({lat_name: xr.DataArray(w.lat, dims="y"), lon_name: xr.DataArray(w.lon, dims="x")}).values
    # Rows beyond the data's last latitude (near the poles) and missing values.
    return pd.DataFrame(out).ffill().bfill().values if fill else out


def write_gray(w, product, values):
    enc = w.side(product)["encoding"]
    px = np.clip(np.rint((np.nan_to_num(values) - enc["offset"]) / enc["scale"]), 0, 65535).astype(np.uint16)
    Image.fromarray(px).save(w.edit / f"{product}.png")


def write_vector(w, product, u, v):
    enc = w.side(product)["encoding"]
    q = lambda x: np.clip(np.rint((np.nan_to_num(x) - enc["offset"]) / enc["scale"]), 0, 65535).astype(np.uint16)
    rgb = np.stack([q(u), q(v), np.zeros_like(q(u))], axis=-1)
    with open(w.edit / f"{product}.png", "wb") as f:
        png.Writer(w.W, w.H, bitdepth=16, greyscale=False).write(f, rgb.reshape(w.H, w.W * 3).tolist())


def air(w, data):
    ds = xr.open_dataset(data / "eraint_uvz.nc").sel(level=850)
    for month, prod in ((1, "wind_jan"), (7, "wind_jul")):
        u = to_grid(ds.u.sel(month=month), "latitude", "longitude", w)
        v = to_grid(ds.v.sel(month=month), "latitude", "longitude", w)
        write_vector(w, prod, u, v)
    sst = xr.open_dataset(data / "ersstv5.nc").sst.sel(time=slice("1971-01-01", "2000-12-31")).mean("time")
    obs = to_grid(sst, "lat", "lon", w, fill=False)
    ocean = (w.field("elevation") < 0) & np.isfinite(obs)
    anomaly = np.zeros_like(obs)
    for y in range(w.H):
        row = ocean[y]
        if row.sum() > 10:
            anomaly[y, row] = obs[y, row] - obs[y, row].mean()
    write_gray(w, "sst_anomaly", anomaly)
    print(f"air: 850 hPa winds, SST departures {np.nanmin(anomaly):.1f}…{np.nanmax(anomaly):.1f} K")


def climate(w, obs):
    t2m = xr.open_dataset(obs / "T2M_ERAINT_rectilinear_grid_2D.nc", decode_times=False).T2M - 273.15
    months = np.stack([to_grid(t2m.isel(time=m), "lat", "lon", w) for m in range(12)])
    write_gray(w, "temperature", months.mean(0))
    write_gray(w, "temperature_range", months.max(0) - months.min(0))
    gpcp = xr.open_dataset(obs / "V22_GPCP.1979-2010.nc", decode_times=False).PREC
    # Month of year of each record (the series starts in January 1979).
    clim = gpcp.groupby(xr.DataArray(np.arange(gpcp.sizes["time"]) % 12, dims="time", name="m")).mean("time")
    rain = np.stack([np.maximum(to_grid(clim.isel(m=m), "lat", "lon", w), 0.0) for m in range(12)])
    days = np.array([31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])[:, None, None]
    monthly_mm = rain * days
    annual = monthly_mm.sum(0)
    north = monthly_mm[3:9].sum(0)
    summer = np.where(w.lat[:, None] >= 0, north, annual - north)
    write_gray(w, "precipitation", annual)
    write_gray(w, "precipitation_warm_share", np.where(annual > 0, summer / np.maximum(annual, 1e-9), 0.5))
    print(f"climate: mean {np.average(months.mean(0), weights=w.area):.1f} °C, rain {np.average(annual, weights=w.area):.0f} mm a year")


def main():
    world, data, obs, which = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), sys.argv[4]
    w = World(world)
    if which == "air":
        air(w, data)
    elif which == "climate":
        climate(w, obs)
    else:
        sys.exit("air or climate")


if __name__ == "__main__":
    main()
