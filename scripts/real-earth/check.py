"""Check a world built on real Earth topography (see run.sh) against the observed Earth.

usage: python3 -I check.py WORLD_DIR DATA_DIR [REPORT.json]

Layer by layer: winds (ERA-Interim 850 hPa, January and July), sea-surface temperature and its
departure from the zonal mean (ERSST v5, 1971–2000), land climate (Köppen–Geiger map, Rubel
et al.; station normals), ocean rainfall by latitude, rivers, biomes. The station normals,
river discharges and zonal ocean rainfall below are approximate, rounded reference values.
"""
import json
import sys
from pathlib import Path

import numpy as np
import xarray as xr
from PIL import Image

sys.path.insert(0, str(Path(__file__).parent))
from gfview import World  # noqa: E402

Image.MAX_IMAGE_PIXELS = None

# Station normals (approximate): lat, lon, annual mean °C, warmest − coldest month °C, mm/yr,
# Köppen class.
STATIONS = {
    "London": (51.5, -0.1, 11.3, 13.5, 620, "Cfb"),
    "Moscow": (55.75, 37.6, 5.8, 25.7, 700, "Dfb"),
    "Yakutsk": (62.0, 129.7, -8.8, 58.0, 240, "Dfd"),
    "Chicago": (41.9, -87.6, 10.0, 27.9, 950, "Dfa"),
    "Beijing": (39.9, 116.4, 12.9, 30.0, 570, "Dwa"),
    "Vancouver": (49.25, -123.1, 10.5, 14.5, 1450, "Cfb"),
    "Reykjavik": (64.1, -21.9, 4.7, 11.0, 850, "Cfc"),
    "Los Angeles": (34.05, -118.25, 18.6, 9.0, 380, "Csb"),
    "Perth": (-31.95, 115.85, 18.7, 11.0, 730, "Csa"),
    "Cape Town": (-33.9, 18.45, 17.0, 8.5, 515, "Csb"),
    "Buenos Aires": (-34.6, -58.4, 17.9, 13.5, 1200, "Cfa"),
    "Cairo": (30.05, 31.25, 22.0, 14.5, 25, "BWh"),
    "Riyadh": (24.7, 46.7, 26.5, 22.0, 100, "BWh"),
    "In Salah": (27.2, 2.5, 25.5, 25.0, 15, "BWh"),
    "Timbuktu": (16.8, -3.0, 29.0, 12.0, 200, "BWh"),
    "Phoenix": (33.45, -112.07, 24.0, 21.0, 200, "BWh"),
    "Alice Springs": (-23.7, 133.9, 21.0, 16.0, 285, "BWh"),
    "Lima": (-12.05, -77.05, 19.4, 7.0, 15, "BWh"),
    "Tehran": (35.7, 51.4, 18.5, 27.0, 230, "BSk"),
    "Dodge City": (37.75, -100.0, 13.0, 26.0, 560, "BSk"),
    "Ulaanbaatar": (47.9, 106.9, -0.4, 40.0, 270, "Dwc"),
    "Astana": (51.15, 71.45, 3.5, 37.0, 330, "Dfb"),
    "Manaus": (-3.1, -60.0, 27.0, 1.8, 2300, "Af"),
    "Singapore": (1.35, 103.8, 27.5, 1.8, 2500, "Af"),
    "Kinshasa": (-4.3, 15.3, 25.3, 4.0, 1400, "Aw"),
    "Mumbai": (19.1, 72.9, 27.2, 6.0, 2200, "Am"),
    "Nairobi": (-1.3, 36.8, 17.8, 4.0, 870, "Cwb"),
    "Lhasa": (29.65, 91.1, 8.0, 17.0, 430, "Cwb"),
    "Utqiagvik": (71.3, -156.8, -11.0, 31.0, 140, "ET"),
    "Summit Greenland": (72.6, -38.5, -30.0, 30.0, 250, "EF"),
    "Vostok": (-78.45, 106.8, -55.0, 35.0, 20, "EF"),
}

# River mouths (lat, lon), mean discharge m³/s and basin area km² (approximate).
RIVERS = {
    "Amazon": (-0.5, -50.0, 209000, 7.0e6),
    "Congo": (-6.0, 12.3, 41000, 3.7e6),
    "Ganges–Brahmaputra": (22.0, 90.0, 38000, 1.7e6),
    "Orinoco": (8.6, -61.0, 33000, 0.95e6),
    "Yangtze": (31.4, 121.8, 30000, 1.8e6),
    "Paraná (La Plata)": (-34.5, -58.0, 25000, 2.6e6),
    "Yenisei": (71.5, 82.5, 19600, 2.6e6),
    "Lena": (72.5, 127.0, 17000, 2.5e6),
    "Mississippi": (29.2, -89.3, 16800, 3.2e6),
    "Mekong": (9.8, 106.6, 16000, 0.8e6),
    "Ob": (66.6, 71.5, 12500, 3.0e6),
    "Amur": (52.9, 141.1, 11400, 1.9e6),
    "Mackenzie": (68.9, -135.5, 10300, 1.8e6),
    "St. Lawrence": (48.6, -69.0, 10100, 1.3e6),
    "Volga (to Caspian)": (45.8, 47.9, 8000, 1.4e6),
    "Niger": (4.8, 6.0, 6000, 2.1e6),
    "Danube": (45.2, 29.6, 6500, 0.8e6),
    "Nile": (31.4, 31.2, 2800, 3.3e6),
}

# Regions for SST departures from the zonal mean: (lat0, lat1, lon0, lon1), sign expected.
SST_REGIONS = {
    "NE Atlantic (Gulf Stream / N Atlantic Drift)": ((45, 65), (-30, 10), +1),
    "Kuroshio": ((30, 40), (130, 150), +1),
    "California Current": ((25, 40), (-130, -115), -1),
    "Humboldt": ((-25, -5), (-85, -72), -1),
    "Benguela": ((-30, -15), (5, 14), -1),
    "Canary": ((15, 30), (-25, -12), -1),
    "W Pacific warm pool": ((-10, 10), (120, 160), +1),
}

KOPPEN_GROUPS = ["A", "BW", "BS", "C", "D", "E"]


def koppen(t_mean, t_range, precip, warm_share):
    """Köppen (Kottek et al. 2006 thresholds) from annual mean, annual range and precipitation."""
    warm = t_mean + t_range / 2
    cold = t_mean - t_range / 2
    pth = 10 * (2 * t_mean + np.where(warm_share >= 2 / 3, 28, np.where(warm_share <= 1 / 3, 0, 14)))
    out = np.full(t_mean.shape, "", dtype=object)
    e = warm < 10
    b = ~e & (precip < pth)
    a = ~e & ~b & (cold >= 18)
    c = ~e & ~b & ~a & (cold > -3)
    d = ~e & ~b & ~a & ~c
    out[e] = np.where(warm[e] < 0, "EF", "ET")
    out[b] = np.where(precip[b] < pth[b] / 2, "BW", "BS")
    out[b] = out[b] + np.where(t_mean[b] >= 18, "h", "k")
    out[a], out[c], out[d] = "A", "C", "D"
    return out


def group(k):
    k = str(k)
    if k.startswith("B"):
        return k[:2]
    return k[:1]


def observed_koppen(w, ref):
    img = np.array(Image.open(ref / "kgcpy" / "kmz_int_reshape.png"))
    names = {}
    for line in (ref / "kgcpy" / "kg_zoneNum.csv").read_text().splitlines()[1:]:
        _, k, n = line.split(",")
        names[int(n)] = k
    H, W = img.shape
    ys = np.clip(np.round((90 - w.lat) * H / 180 - 0.5).astype(int), 0, H - 1)
    xs = np.clip(np.round((w.lon + 180) * W / 360 - 0.5).astype(int), 0, W - 1)
    codes = img[np.ix_(ys, xs)]
    return np.vectorize(lambda c: names.get(int(c), "Ocean"))(codes)


def stations(w):
    t, tr, p, ws = (w.field(x) for x in ("temperature", "temperature_range", "precipitation", "precipitation_warm_share"))
    elev = w.field("elevation")
    bk, brow = w.partition("biome")
    rows = []
    for name, (la, lo, t_obs, r_obs, p_obs, k_obs) in STATIONS.items():
        y, x = w.nearest_land(la, lo)
        k = koppen(np.array([t[y, x]]), np.array([tr[y, x]]), np.array([p[y, x]]), np.array([ws[y, x]]))[0]
        rows.append(
            {
                "station": name,
                "elev_m": round(float(elev[y, x])),
                "t": [round(float(t[y, x]), 1), t_obs],
                "range": [round(float(tr[y, x]), 1), r_obs],
                "precip": [round(float(p[y, x])), p_obs],
                "koppen": [k, k_obs],
                "biome": brow.get(int(bk[y, x]), {}).get("name", "?"),
            }
        )
    return rows


def koppen_confusion(w, ref):
    t, tr, p, ws = (w.field(x) for x in ("temperature", "temperature_range", "precipitation", "precipitation_warm_share"))
    land = w.field("elevation") >= 0
    obs = observed_koppen(w, ref)
    model = koppen(t, tr, p, ws)
    both = land & (obs != "Ocean")
    g_obs = np.vectorize(group)(obs)
    g_mod = np.vectorize(group)(model)
    a = w.area
    tot = a[both].sum()
    shares = {g: [round(float(a[both & (g_mod == g)].sum() / tot), 3), round(float(a[both & (g_obs == g)].sum() / tot), 3)] for g in KOPPEN_GROUPS}
    agree = round(float(a[both & (g_mod == g_obs)].sum() / tot), 3)
    conf = {go: {gm: round(float(a[both & (g_obs == go) & (g_mod == gm)].sum() / tot), 3) for gm in KOPPEN_GROUPS} for go in KOPPEN_GROUPS}
    bh = {k: [round(float(a[both & (model == k)].sum() / tot), 3), round(float(a[both & (obs == k)].sum() / tot), 3)] for k in ("BWh", "BWk", "BSh", "BSk")}
    return {"shares_model_obs": shares, "b_subtypes_model_obs": bh, "agreement": agree, "confusion_obs_to_model": conf}, g_obs, g_mod


def sst(w, ref):
    ds = xr.open_dataset(ref / "ersstv5.nc")
    clim = ds.sst.sel(time=slice("1971-01-01", "2000-12-31")).mean("time")
    clim = clim.assign_coords(lon=((clim.lon + 180) % 360) - 180).sortby("lon")
    obs = clim.interp(lat=xr.DataArray(w.lat, dims="y"), lon=xr.DataArray(w.lon, dims="x"), method="nearest").values
    model = w.field("temperature")
    ocean = (w.field("elevation") < 0) & np.isfinite(obs) & (np.abs(w.lat)[:, None] < 60)
    a = w.area
    diff = model - obs
    bias = float((diff * a)[ocean].sum() / a[ocean].sum())
    rmse = float(np.sqrt(((diff**2) * a)[ocean].sum() / a[ocean].sum()))
    zonal = {}
    for lo, hi in [(-60, -40), (-40, -20), (-20, 0), (0, 20), (20, 40), (40, 60)]:
        m = ocean & (w.lat[:, None] >= lo) & (w.lat[:, None] < hi)
        zonal[f"{lo}..{hi}"] = [round(float((model * a)[m].sum() / a[m].sum()), 1), round(float((obs * a)[m].sum() / a[m].sum()), 1)]

    def departure(f):
        out = np.full(f.shape, np.nan)
        for y in range(f.shape[0]):
            row = ocean[y]
            if row.sum() > 10:
                out[y, row] = f[y, row] - f[y, row].mean()
        return out

    dm, do = departure(model), departure(obs)
    an = w.field("sst_anomaly")
    ok = ocean & np.isfinite(dm) & np.isfinite(do)
    corr_t = float(np.corrcoef(dm[ok], do[ok])[0, 1])
    corr_a = float(np.corrcoef(an[ok], do[ok])[0, 1])
    regions = {}
    for name, ((la0, la1), (lo0, lo1), sign) in SST_REGIONS.items():
        m = ok & (w.lat[:, None] >= la0) & (w.lat[:, None] <= la1) & (w.lon[None, :] >= lo0) & (w.lon[None, :] <= lo1)
        regions[name] = {"model_departure": round(float(np.nanmean(dm[m])), 1), "observed_departure": round(float(np.nanmean(do[m])), 1), "model_anomaly_product": round(float(np.nanmean(an[m])), 1), "expected_sign": sign}
    return {"bias_c": round(bias, 2), "rmse_c": round(rmse, 2), "zonal_model_obs": zonal, "departure_corr_temperature": round(corr_t, 2), "departure_corr_sst_anomaly": round(corr_a, 2), "regions": regions}


def winds(w, ref):
    ds = xr.open_dataset(ref / "eraint_uvz.nc").sel(level=850)
    out = {}
    for month, prod in ((1, "wind_jan"), (7, "wind_jul")):
        u, v = w.vector(prod)
        ou = ds.u.sel(month=month).interp(latitude=xr.DataArray(w.lat, dims="y"), longitude=xr.DataArray(w.lon, dims="x")).values
        ov = ds.v.sel(month=month).interp(latitude=xr.DataArray(w.lat, dims="y"), longitude=xr.DataArray(w.lon, dims="x")).values
        ocean = w.field("elevation") < 0
        bands = {}
        for lo, hi in [(-70, -50), (-50, -30), (-30, -10), (-10, 10), (10, 30), (30, 50), (50, 70)]:
            m = ocean & (w.lat[:, None] >= lo) & (w.lat[:, None] < hi) & np.isfinite(ou)
            bands[f"{lo}..{hi}"] = {"u_model": round(float(u[m].mean()), 1), "u_era850": round(float(ou[m].mean()), 1), "v_model": round(float(v[m].mean()), 1), "v_era850": round(float(ov[m].mean()), 1)}
        ok = np.isfinite(ou)
        dir_agree = float(((u * ou + v * ov) > 0)[ok & ocean].mean())
        # Monsoon: Arabian Sea and Bay of Bengal.
        m = ocean & (w.lat[:, None] >= 8) & (w.lat[:, None] <= 18) & (w.lon[None, :] >= 55) & (w.lon[None, :] <= 95)
        out[prod] = {"bands": bands, "direction_agreement_ocean": round(dir_agree, 2), "north_indian_ocean_uv_model": [round(float(u[m].mean()), 1), round(float(v[m].mean()), 1)], "north_indian_ocean_uv_era850": [round(float(ou[m].mean()), 1), round(float(ov[m].mean()), 1)]}
    return out


def rivers(w):
    r = w.rivers()
    elev = w.field("elevation").ravel()
    wk, _ = w.partition("watersheds")
    wkf = wk.ravel()
    a = (w.area * (2 * np.pi * 6371 / w.W) ** 2).ravel()
    basin_area = np.bincount(wkf, weights=a)
    mouths = [(e["from"], e["to"], e["attrs"]["discharge_m3_s"]) for e in r["edges"]]
    q_into = {}
    for f, t, q in mouths:
        q_into[f] = max(q_into.get(f, 0.0), q)
    cells = np.array(list(q_into.keys()))
    qs = np.array([q_into[c] for c in cells])
    cy, cx = cells // w.W, cells % w.W
    out = []
    for name, (la, lo, q_obs, a_obs) in RIVERS.items():
        y, x = w.cell(la, lo)
        dx = np.minimum(np.abs(cx - x), w.W - np.abs(cx - x))
        near = (np.abs(cy - y) <= 12) & (dx <= 12)
        if not near.any():
            out.append({"river": name, "q": [None, q_obs]})
            continue
        i = np.argmax(np.where(near, qs, -1))
        c = int(cells[i])
        key = int(wkf[c])
        out.append({"river": name, "q": [round(float(qs[i])), q_obs], "basin_km2": [round(float(basin_area[key])) if key else None, a_obs]})
    top = np.argsort(-qs)
    seen, biggest = [], []
    for i in top:
        k = int(wkf[cells[i]])
        if k in seen:
            continue
        seen.append(k)
        biggest.append({"lat": round(float(w.lat[cy[i]]), 1), "lon": round(float(w.lon[cx[i]]), 1), "q": round(float(qs[i])), "basin_km2": round(float(basin_area[k]))})
        if len(biggest) == 12:
            break
    return {"named": out, "largest_by_discharge": biggest}


# Zonal-mean ocean precipitation, mm/yr, approximately as in the GPCP climatology.
OCEAN_ZONAL_PRECIP = {(-60, -50): 900, (-50, -40): 1100, (-40, -30): 900, (-30, -20): 700, (-20, -10): 1000, (-10, 0): 1500, (0, 10): 2600, (10, 20): 1300, (20, 30): 700, (30, 40): 1000, (40, 50): 1200, (50, 60): 1100}


def ocean_precip(w):
    p = w.field("precipitation")
    ocean = w.field("elevation") < 0
    a = w.area
    out = {}
    for (lo, hi), ref in OCEAN_ZONAL_PRECIP.items():
        m = ocean & (w.lat[:, None] >= lo) & (w.lat[:, None] < hi)
        out[f"{lo}..{hi}"] = [round(float((p * a)[m].sum() / a[m].sum())), ref]
    return out


def biome_shares(w):
    bk, brow = w.partition("biome")
    land = w.field("elevation") >= 0
    a = w.area
    tot = a[land].sum()
    return {brow[k]["slug"]: round(float(a[land & (bk == k)].sum() / tot), 3) for k in sorted(brow) if a[land & (bk == k)].sum() > 0}


def summarize(report):
    lines = ["Stations (model / observed):"]
    for s in report["stations"]:
        lines.append(f"  {s['station']:17s} {s['elev_m']:5d} m  T {s['t'][0]:6.1f}/{s['t'][1]:6.1f}  range {s['range'][0]:5.1f}/{s['range'][1]:5.1f}  P {s['precip'][0]:5d}/{s['precip'][1]:5d}  Köppen {s['koppen'][0]:4s}/{s['koppen'][1]:4s}  {s['biome']}")
    st = [s for s in report["stations"] if s["station"] != "Vostok"]

    def rmse(key, log=False):
        d = [(np.log1p(s[key][0]) - np.log1p(s[key][1])) if log else (s[key][0] - s[key][1]) for s in st]
        return round(float(np.sqrt(np.mean(np.square(d)))), 2)

    lines.append(f"  RMSE: temperature {rmse('t')} °C, range {rmse('range')} °C, log precipitation {rmse('precip', True)}")
    k = report["koppen"]
    lines.append(f"Köppen major groups, area share of land model/observed: {k['shares_model_obs']}")
    lines.append(f"  desert and steppe subtypes: {k['b_subtypes_model_obs']}; cell agreement {k['agreement']}")
    for m in ("wind_jan", "wind_jul"):
        wd = report["winds"][m]
        lines.append(f"{m}: direction agrees with ERA-Interim over {wd['direction_agreement_ocean']:.0%} of the ocean; northern Indian Ocean (u, v) model {wd['north_indian_ocean_uv_model']} vs {wd['north_indian_ocean_uv_era850']} m/s")
    ss = report["sst"]
    lines.append(f"SST: bias {ss['bias_c']} °C, RMSE {ss['rmse_c']} °C; departure-from-zonal-mean correlation {ss['departure_corr_sst_anomaly']} (anomaly product), {ss['departure_corr_temperature']} (temperature)")
    for name, r in ss["regions"].items():
        lines.append(f"  {name}: model {r['model_anomaly_product']:+.1f} °C, observed {r['observed_departure']:+.1f} °C")
    lines.append(f"Ocean rainfall by latitude, model/reference mm: {report['ocean_precip']}")
    lines.append("Rivers at the mouth, model/observed m³/s:")
    for r in report["rivers"]["named"]:
        lines.append(f"  {r['river']:20s} {r['q'][0]}/{r['q'][1]}")
    lines.append(f"Biomes, share of land: {report['biomes']}")
    return "\n".join(lines)


def main():
    w = World(sys.argv[1])
    ref = Path(sys.argv[2])
    report = {"stations": stations(w)}
    report["koppen"], _, _ = koppen_confusion(w, ref)
    report["sst"] = sst(w, ref)
    report["winds"] = winds(w, ref)
    report["ocean_precip"] = ocean_precip(w)
    report["rivers"] = rivers(w)
    report["biomes"] = biome_shares(w)
    if len(sys.argv) > 3:
        Path(sys.argv[3]).write_text(json.dumps(report, indent=1))
    print(summarize(report))


if __name__ == "__main__":
    main()
