"""The Earth A/B page: each layer run on real inputs, its result beside the real one.

usage: python3 -I ab_page.py OUT_DIR DATA_DIR REPO_DIR
Reads the worlds ab.sh makes in OUT_DIR (chain, air, real, plates), their reports and checks,
and writes OUT_DIR/page (index.html and img/). Maps come from geoforge's own renders, so real
and model fields share palettes and keys.
"""
import html
import json
import re
import shutil
import sys
from pathlib import Path

import numpy as np
from PIL import Image

sys.path.insert(0, str(Path(__file__).parent))
from gfview import World  # noqa: E402

e = html.escape
OUT, DATA, REPO = (Path(a) for a in sys.argv[1:4])
PAGE = OUT / "page"

# Boxes for the mountain check: (lat0, lat1, lon0, lon1).
RANGES = {
    "Rockies": (37, 50, -115, -104),
    "Colorado Plateau": (35, 39, -112, -107),
    "Appalachians": (35, 45, -84, -74),
    "Andes, Altiplano": (-22, -14, -70, -65),
    "Himalaya and Tibet": (28, 36, 80, 95),
    "Alps": (44, 48, 6, 15),
    "Ethiopian Highlands": (7, 14, 36, 40),
    "East African Plateau": (-5, 2, 30, 37),
    "Urals": (50, 65, 57, 61),
    "Brazilian Highlands": (-23, -15, -48, -42),
    "Iranian Plateau": (29, 36, 52, 60),
    "Great Dividing Range": (-37, -20, 145, 152),
}


def report_parts(world):
    """CSS, the row-labelling script, and each figure's key, from a world's report page."""
    page = (OUT / world / "report" / "index.html").read_text()
    css = re.search(r"<style>(.*?)</style>", page, re.S).group(1)
    script = re.search(r"<script>(.*?)</script>", page, re.S).group(1)
    keys = {m.group(1): m.group(2) for m in re.finditer(r'src="img/([a-z_]+)\.jpg" alt="[^"]*"><figcaption><b>[^<]*</b>(.*?)</figcaption>', page, re.S)}
    return css, script, keys


CSS, SCRIPT, KEYS = report_parts("real")
_, _, CHAIN_KEYS = report_parts("chain")


def img(world, name):
    """Copy a world's report image into the page; its path there."""
    src = OUT / world / "report" / "img" / f"{name}.jpg"
    dst = PAGE / "img" / f"{world}_{name}.jpg"
    shutil.copy(src, dst)
    return f"img/{dst.name}"


def fig(src, title, note=""):
    return f'<figure><img loading="lazy" src="{src}" alt="{e(title)}"><figcaption><b>{e(title)}</b>{(" " + e(note)) if note else ""}</figcaption></figure>'


def pair(name, left, right, key_from=None):
    """Two figures and their shared key: left and right are (world or src, title, note)."""
    figs = "".join(fig(src if src.startswith("img/") else img(src, name), t, n) for src, t, n in (left, right))
    key = (key_from or KEYS).get(name, "")
    return f'<div class="pair">{figs}</div>{key and "<div class=card style=margin-top:8px>" + key + "</div>"}'


def num(x, d=1):
    return "—" if x is None else f"{x:,.{d}f}".replace("-", "−")


# ------------------------------------------------------------------ mountains on real plates
def mountains():
    chain = World(OUT / "chain")
    real = chain.field("elevation")
    H, W = 360, 720
    model = np.fromfile(REPO / "target" / "known-cases" / "tectonics_elevation.f32", dtype="<f4").astype(float).reshape(H, W)
    ys = ((np.arange(H) + 0.5) * chain.H / H).astype(int)
    xs = ((np.arange(W) + 0.5) * chain.W / W).astype(int)
    real = real[np.ix_(ys, xs)]
    lat = 90 - (np.arange(H) + 0.5) * 180 / H
    lon = -180 + (np.arange(W) + 0.5) * 360 / W
    area = np.repeat(np.cos(np.radians(lat))[:, None], W, axis=1)

    def stats(f):
        land = f > 0
        a = area[land]
        return {
            "land": a.sum() / area.sum(),
            "mean": (f[land] * a).sum() / a.sum(),
            "above1": a[f[land] > 1000].sum() / a.sum(),
            "above2": a[f[land] > 2000].sum() / a.sum(),
        }

    r, m = stats(real), stats(model)
    both = (real > 0) & (model > 0)
    corr = float(np.corrcoef(real[both], model[both])[0, 1])
    rows = "".join(
        f'<tr><td>{e(k)}</td><td class="n">{num(r[k_] * s, 0)}{u}</td><td class="n">{num(m[k_] * s, 0)}{u}</td></tr>'
        for k, k_, s, u in [("Land, share of the surface", "land", 100, "%"), ("Mean height of land", "mean", 1, " m"), ("Land above 1 km", "above1", 100, "%"), ("Land above 2 km", "above2", 100, "%")]
    )
    ranges = ""
    for name, (la0, la1, lo0, lo1) in RANGES.items():
        box = (lat[:, None] >= la0) & (lat[:, None] <= la1) & (lon[None, :] >= lo0) & (lon[None, :] <= lo1)
        rb, mb = box & (real > 0), box & (model > 0)
        rm = float(np.mean(real[rb])) if rb.any() else None
        mm = float(np.mean(model[mb])) if mb.any() else None
        r95 = float(np.percentile(real[rb], 95)) if rb.any() else None
        m95 = float(np.percentile(model[mb], 95)) if mb.any() else None
        ok = rm and mm and 0.6 <= mm / rm <= 1.6
        ranges += f'<tr><td>{e(name)}</td><td class="n">{num(rm, 0)} / {num(r95, 0)}</td><td class="n {"ok" if ok else "miss"}">{num(mm, 0)} / {num(m95, 0)}</td></tr>'

    # The model's provinces on the real plates, in the provinces map's colours.
    src = (REPO / "crates/libs/geoforge-core/src/surface.rs").read_text()
    block = src[src.index("category!(GeologicProvince"):]
    colours = {int(n): (int(a), int(b), int(c)) for n, a, b, c in re.findall(r"\w+ = (\d+), \"[^\"]+\", \[(\d+), (\d+), (\d+)\]", block[: block.index("});")])}
    names = {int(n): t for n, t in re.findall(r"\w+ = (\d+), \"([^\"]+)\"", block[: block.index("});")])}
    prov = np.fromfile(REPO / "target" / "known-cases" / "tectonics_provinces.f32", dtype="<f4").astype(int).reshape(H, W)
    rgb = np.zeros((H, W, 3), np.uint8)
    for k, c in colours.items():
        rgb[prov == k] = c
    Image.fromarray(rgb).resize((1200, 600), Image.NEAREST).save(PAGE / "img" / "plates_provinces.jpg", quality=88)
    present = sorted({int(k) for k in np.unique(prov) if k in names})
    key = '<ul class="key">' + "".join(f'<li><i style="background:#{colours[k][0]:02x}{colours[k][1]:02x}{colours[k][2]:02x}"></i>{e(names[k])}</li>' for k in present) + "</ul>"
    body = pair("elevation", ("chain", "Real: ETOPO1", "Earth's measured relief."), ("plates", "Model: elevation on Earth's real plates", "PB2002 plates and Euler poles in, mountains out (erosion on)."))
    body += f'<div class="pair">{fig("img/plates_provinces.jpg", "Model: geologic provinces on the real plates", "What the elevation layer builds from.")}<div class="card"><h3>Province key</h3>{key}</div></div>'
    body += f"""<div class="cols"><div class="card"><h3>Hypsometry</h3><table><thead><tr><th></th><th>Real</th><th>Model</th></tr></thead><tbody>{rows}</tbody></table>
<p class="note">Height correlation where both are land: {corr:.2f}.</p></div>
<div class="card"><h3>Ranges: mean / 95th percentile of land, m</h3><table><thead><tr><th>Range</th><th>Real</th><th>Model</th></tr></thead><tbody>{ranges}</tbody></table>
<p class="note">Green: the model's mean within 0.6–1.6 of the real one. Ranges far from plate edges (Rockies, Appalachians, Ethiopia, Brazil) have nothing to build them yet: roadmap 8e.</p></div></div>"""
    return body, corr, r, m


# ------------------------------------------------------------------ air and sea
def air():
    c = json.loads((OUT / "chain" / "earth_check.json").read_text())
    w = c["winds"]
    body = ""
    for prod, month in (("wind_jan", "January"), ("wind_jul", "July")):
        x = w[prod]
        note = f"Direction agrees over {100 * x['direction_agreement_ocean']:.0f}% of the ocean; northern Indian Ocean (u, v) {x['north_indian_ocean_uv_model']} m/s against {x['north_indian_ocean_uv_era850']}."
        body += f"<h3>Surface winds, {month}</h3>" + pair(prod, ("real", "Real: ERA-Interim, 850 hPa", "Above the friction layer, so stronger and less turned than surface winds."), ("chain", "Model: on real topography", note))
    s = c["sst"]
    rows = "".join(f'<tr><td>{e(k)}</td><td class="n">{v["observed_departure"]:+.1f}</td><td class="n {"ok" if np.sign(v["model_anomaly_product"]) == np.sign(v["observed_departure"]) and abs(v["model_anomaly_product"]) >= 0.5 * abs(v["observed_departure"]) else "miss"}">{v["model_anomaly_product"]:+.1f}</td></tr>' for k, v in s["regions"].items())
    body += "<h3>Sea-surface temperature, departure from each latitude's mean</h3>" + pair("sst_anomaly", ("real", "Real: ERSST v5, 1971–2000", ""), ("chain", "Model: currents on real topography", f"Pattern correlation {s['departure_corr_sst_anomaly']:.2f}."))
    body += f'<div class="card"><h3>Boundary currents, K from the latitude\'s mean</h3><table><thead><tr><th>Region</th><th>Real</th><th>Model</th></tr></thead><tbody>{rows}</tbody></table><p class="note">Green: right sign and at least half the size. The cold eastern currents run at about half their real strength, as the gyre question suspected.</p></div>'
    return body, w, s


# ------------------------------------------------------------------ climate on real air
def land_metrics(world, real):
    m = World(OUT / world)
    land = real.field("elevation") >= 0
    a = real.area[land]
    out = {}
    for prod in ("temperature", "temperature_range"):
        d = (m.field(prod) - real.field(prod))[land]
        out[prod] = ((d * a).sum() / a.sum(), np.sqrt((d * d * a).sum() / a.sum()))
    pm, pr = m.field("precipitation")[land], real.field("precipitation")[land]
    lm, lr = np.log(pm + 50), np.log(pr + 50)
    out["precip_corr"] = float(np.corrcoef(lm, lr)[0, 1])
    out["precip_ratio"] = float(np.exp(np.median(lm - lr)))
    dry = pr < 250
    out["too_wet_deserts"] = float(a[dry & (pm > 500)].sum() / max(a[dry].sum(), 1e-9))
    return out


def climate():
    real = World(OUT / "real")
    ch, ai = land_metrics("chain", real), land_metrics("air", real)
    rows = f"""<tr><td>Temperature bias / RMSE, °C</td><td class="n">{ch['temperature'][0]:+.1f} / {ch['temperature'][1]:.1f}</td><td class="n">{ai['temperature'][0]:+.1f} / {ai['temperature'][1]:.1f}</td></tr>
<tr><td>Seasonal range bias / RMSE, K</td><td class="n">{ch['temperature_range'][0]:+.1f} / {ch['temperature_range'][1]:.1f}</td><td class="n">{ai['temperature_range'][0]:+.1f} / {ai['temperature_range'][1]:.1f}</td></tr>
<tr><td>Rain: pattern correlation (log)</td><td class="n">{ch['precip_corr']:.2f}</td><td class="n">{ai['precip_corr']:.2f}</td></tr>
<tr><td>Rain: typical model ÷ real</td><td class="n">{ch['precip_ratio']:.2f}</td><td class="n">{ai['precip_ratio']:.2f}</td></tr>
<tr><td>Real deserts (under 250 mm) given over 500 mm</td><td class="n">{100 * ch['too_wet_deserts']:.0f}%</td><td class="n">{100 * ai['too_wet_deserts']:.0f}%</td></tr>"""
    sc = {s["station"]: s for s in json.loads((OUT / "chain" / "earth_check.json").read_text())["stations"]}
    sa = {s["station"]: s for s in json.loads((OUT / "air" / "earth_check.json").read_text())["stations"]}
    st = ""
    for name in ["Lima", "Phoenix", "Tehran", "Riyadh", "Perth", "Yakutsk", "Astana", "Moscow", "London", "Mumbai", "Manaus", "Nairobi"]:
        c, a = sc[name], sa[name]
        st += f'<tr><td>{e(name)}</td><td class="n">{c["t"][1]:.1f} / {c["precip"][1]:,} / {e(c["koppen"][1])}</td><td class="n">{c["t"][0]:.1f} / {c["precip"][0]:,} / {e(c["koppen"][0])}</td><td class="n">{a["t"][0]:.1f} / {a["precip"][0]:,} / {e(a["koppen"][0])}</td></tr>'
    body = ""
    for prod, title, note in [
        ("temperature", "Annual mean temperature", "Real: ERA-Interim 2 m, 2005."),
        ("temperature_range", "Seasonal range", "Real: ERA-Interim 2 m, 2005."),
        ("precipitation", "Annual precipitation", "Real: GPCP v2.2, 1979–2010."),
    ]:
        body += f"<h3>{title}</h3>" + pair(prod, ("real", f"Real: {note.split(': ')[1]}", ""), ("air", "Model: climate on real topography, winds and sea temperatures", ""))
    body += f"""<div class="cols"><div class="card"><h3>Climate layer over land</h3><table><thead><tr><th></th><th>Model air</th><th>Real air</th></tr></thead><tbody>{rows}</tbody></table>
<p class="note">Left, the climate on the model's own winds and currents; right, on ERA-Interim's winds and ERSST's sea temperatures. What stays wrong on the right is the climate layer's own.</p></div>
<div class="card"><h3>Stations: °C / mm / Köppen</h3><div class="scroll"><table><thead><tr><th>Station</th><th>Real</th><th>Model air</th><th>Real air</th></tr></thead><tbody>{st}</tbody></table></div></div></div>"""
    return body, ch, ai


# ------------------------------------------------------------------ derived on a real climate
def koppen_page():
    """The observed Köppen map in the model's colours."""
    src = (REPO / "crates/libs/geoforge-analysis/src/koppen.rs").read_text()
    colour = {code: (int(r), int(g), int(b)) for code, _, r, g, b in re.findall(r'\("(\w+)", "([^"]+)", \[(\d+), (\d+), (\d+)\]\)', src)}
    names = {}
    for line in (DATA / "kgcpy" / "kg_zoneNum.csv").read_text().splitlines()[1:]:
        _, k, n = line.split(",")
        names[int(n)] = k
    zones = np.array(Image.open(DATA / "kgcpy" / "kmz_int_reshape.png"))
    H, W = 600, 1200
    ys = ((np.arange(H) + 0.5) * zones.shape[0] / H).astype(int)
    xs = ((np.arange(W) + 0.5) * zones.shape[1] / W).astype(int)
    z = zones[np.ix_(ys, xs)]
    rgb = np.full((H, W, 3), (225, 232, 240), np.uint8)
    for n, k in names.items():
        if k in colour:
            rgb[z == n] = colour[k]
    Image.fromarray(rgb).save(PAGE / "img" / "observed_koppen.jpg", quality=88)
    checks = {w: json.loads((OUT / w / "earth_check.json").read_text())["koppen"] for w in ("chain", "air", "real")}
    shares = "".join(
        f'<tr><td>{g}</td><td class="n">{100 * checks["real"]["shares_model_obs"][g][1]:.0f}%</td>' + "".join(f'<td class="n">{100 * checks[w]["shares_model_obs"][g][0]:.0f}%</td>' for w in ("real", "air", "chain")) + "</tr>"
        for g in ["A", "BW", "BS", "C", "D", "E"]
    )
    agree = "".join(f'<td class="n">{100 * checks[w]["agreement"]:.0f}%</td>' for w in ("real", "air", "chain"))
    body = f"""<div class="pair">{fig("img/observed_koppen.jpg", "Real: Köppen–Geiger map (Rubel et al.)", "")}{fig(img("real", "koppen"), "Model classifier on the real climate", "ERA-Interim temperature and GPCP rain, classified by geoforge.")}</div>
<div class="pair">{fig(img("air", "koppen"), "Model climate on real air", "")}{fig(img("chain", "koppen"), "Model climate on model air", "")}</div>
<div class="card" style="margin-top:8px">{KEYS.get("koppen", "")}</div>
<div class="card" style="margin-top:14px"><h3>Köppen groups, share of land</h3><div class="scroll"><table><thead><tr><th>Group</th><th>Real map</th><th>Classifier on real climate</th><th>Model on real air</th><th>Model on model air</th></tr></thead><tbody>{shares}<tr><td>Cells agreeing with the map</td><td></td>{agree}</tr></tbody></table></div>
<p class="note">The classifier on the real climate is the ceiling: what's left there is the classifier and the data's own resolution (2005 temperature, 2.5° rain).</p></div>"""
    return body, checks


def rivers_page():
    rc = {r["river"]: r for r in json.loads((OUT / "real" / "earth_check.json").read_text())["rivers"]["named"]}
    cc = {r["river"]: r for r in json.loads((OUT / "chain" / "earth_check.json").read_text())["rivers"]["named"]}
    rows = ""
    for name, r in rc.items():
        obs = r["q"][1]
        q_real, q_model = r["q"][0], cc[name]["q"][0]
        cls = lambda q: "ok" if q and 0.5 <= q / obs <= 2 else "miss"
        rows += f'<tr><td>{e(name)}</td><td class="n">{obs:,}</td><td class="n {cls(q_real)}">{num(q_real, 0)}</td><td class="n {cls(q_model)}">{num(q_model, 0)}</td></tr>'
    body = pair("rivers", ("real", "Model rivers on real rain", "GPCP rain on ETOPO relief."), ("chain", "Model rivers on model rain", ""))
    body += f'<div class="card"><h3>Discharge at the mouth, m³/s</h3><div class="scroll"><table><thead><tr><th>River</th><th>Real</th><th>On real rain</th><th>On model rain</th></tr></thead><tbody>{rows}</tbody></table></div><p class="note">Green within a factor of two. Where both columns miss, it is the routing (gorges the grid cannot resolve, filled basins); where only the right misses, it is the rain.</p></div>'
    biomes = pair("biome", ("real", "Model biomes on the real climate", ""), ("chain", "Model biomes on the model climate", ""))
    return body, biomes, rc, cc


def sky_page():
    """The chain world's sky: the real Solar System on Earth's real calendar."""
    page = (OUT / "chain" / "report" / "index.html").read_text()
    section = re.search(r'<section class="band" id="sky"><div class="wrap"><h2>[^<]*</h2>(.*?)</div></section>', page, re.S).group(1)
    return re.sub(r'src="img/([a-z_]+)\.jpg"', lambda m: f'src="{img("chain", m.group(1))}"', section)


def main():
    if PAGE.exists():
        shutil.rmtree(PAGE)
    (PAGE / "img").mkdir(parents=True)
    m_html, corr, rh, mh = mountains()
    a_html, winds, sst = air()
    c_html, ch, ai = climate()
    k_html, kc = koppen_page()
    r_html, b_html, rc, cc = rivers_page()
    s_html = sky_page()
    ok = lambda q, o: q and 0.5 <= q / o <= 2
    rivers_real = sum(1 for r in rc.values() if ok(r["q"][0], r["q"][1]))
    rivers_chain = sum(1 for r in cc.values() if ok(r["q"][0], r["q"][1]))
    summary = f"""<tr><td>Elevation from plates</td><td>PB2002 plates</td><td>ETOPO1</td><td>Height correlation {corr:.2f}; land above 1 km {100 * mh['above1']:.0f}% vs {100 * rh['above1']:.0f}%</td></tr>
<tr><td>Winds</td><td>ETOPO1</td><td>ERA-Interim 850 hPa</td><td>Direction right over {100 * winds['wind_jan']['direction_agreement_ocean']:.0f}% / {100 * winds['wind_jul']['direction_agreement_ocean']:.0f}% of the ocean (Jan / Jul)</td></tr>
<tr><td>Currents → SST</td><td>ETOPO1</td><td>ERSST v5</td><td>Departure pattern correlation {sst['departure_corr_sst_anomaly']:.2f}; cold currents half strength</td></tr>
<tr><td>Climate</td><td>ETOPO1 + real winds and SST</td><td>ERA-Interim 2 m, GPCP</td><td>Temperature RMSE {ai['temperature'][1]:.1f} °C (on model air {ch['temperature'][1]:.1f}); rain correlation {ai['precip_corr']:.2f} ({ch['precip_corr']:.2f})</td></tr>
<tr><td>Köppen classifier</td><td>Real temperature and rain</td><td>Köppen–Geiger map</td><td>{100 * kc['real']['agreement']:.0f}% of land agrees (model climate on real air {100 * kc['air']['agreement']:.0f}%, on model air {100 * kc['chain']['agreement']:.0f}%)</td></tr>
<tr><td>Rivers</td><td>ETOPO1 + GPCP rain</td><td>Discharge of 18 rivers</td><td>{rivers_real} of 18 within a factor of two (on model rain {rivers_chain})</td></tr>"""
    page = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Earth A/B, layer by layer</title>
<link rel="preconnect" href="https://fonts.googleapis.com"><link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;600&family=Spectral:wght@500&display=swap" rel="stylesheet">
<style>{CSS}
.pair{{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,420px),1fr));gap:14px;margin-top:12px}}</style></head><body>
<header class="band cover"><div class="wrap"><h1>Earth A/B, layer by layer</h1>
<p class="lede">Each layer gets Earth's real inputs and is set beside Earth's real result, so a miss belongs to that layer alone. Left is always the real Earth; right is geoforge. Maps are drawn by geoforge's own renderer for both, with the same key. Grid 1024 × 512 (39 km cells); the plates run 720 × 360.</p>
<img src="{img("real", "koppen")}" alt="Köppen climates classified from the real climate"></div></header>
<nav class="nav"><a href="#summary">Summary</a><a href="#mountains">Mountains</a><a href="#air">Air and sea</a><a href="#climate">Climate</a><a href="#koppen">Köppen</a><a href="#rivers">Rivers</a><a href="#biomes">Biomes</a><a href="#sky">Sky</a></nav>
<section class="band" id="summary"><div class="wrap"><h2>Each layer on real inputs</h2><div class="scroll"><table class="stack"><thead><tr><th>Layer</th><th>Fed</th><th>Compared with</th><th>Result</th></tr></thead><tbody>{summary}</tbody></table></div>
<p class="note">Not here yet: ocean currents (no surface-current data in hand), hazards (no storm or quake catalogue), biomes (no observed map; shown on both climates instead).</p></div></section>
<section class="band" id="mountains"><div class="wrap"><h2>Mountains from Earth's real plates</h2>{m_html}</div></section>
<section class="band" id="air"><div class="wrap"><h2>Air and sea on real topography</h2>{a_html}</div></section>
<section class="band" id="climate"><div class="wrap"><h2>Climate on real topography, winds and sea temperatures</h2>{c_html}</div></section>
<section class="band" id="koppen"><div class="wrap"><h2>Köppen climates</h2>{k_html}</div></section>
<section class="band" id="rivers"><div class="wrap"><h2>Rivers on real rain</h2>{r_html}</div></section>
<section class="band" id="biomes"><div class="wrap"><h2>Biomes on the real climate</h2>{b_html}</div></section>
<section class="band" id="sky"><div class="wrap"><h2>Earth's sky in the real Solar System, day 30</h2>
<p class="note">Planets, the Moon and the eclipses come from published orbital elements set at the 2026 March equinox; the stars are still drawn from the galaxy model, not a catalogue.</p>{s_html}</div></section>
<script>{SCRIPT}</script></body></html>"""
    (PAGE / "index.html").write_text(page)
    print(f"written {PAGE / 'index.html'}")


if __name__ == "__main__":
    main()
