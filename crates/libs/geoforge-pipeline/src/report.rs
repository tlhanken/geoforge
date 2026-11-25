//! A report on a world, to take its pulse at a glance: the world whole in true colour, its
//! main maps by layer in pipeline order (each with a key) and its numbers against Earth's.
//! The full report ([`Options::full`]) adds every other map (hazards, seasons, soils,
//! resources) and the system's bodies. Places, their weather and the sky have commands of their own (`place`,
//! `almanac`, `weather`, `sky`).
//!
//! [`build`] makes the page and its images; `geoforge_store::report::write` puts them in a
//! folder (`index.html`, `img/`).

use crate::{stats, views, PipelineError, Result, View};
use geoforge_almanac::ephemeris::{real_date, Calendar};
use geoforge_core::grid::{AnyGrid, Grid};
use geoforge_core::product::{ProductId, Products};
use geoforge_core::BodyContext;
use geoforge_store::render::{self, Legend};
use geoforge_store::World;
use image::RgbImage;
use std::fmt::Write;

/// What goes into a report.
#[derive(Debug, Clone, PartialEq)]
pub struct Options {
    /// Width of the maps, px.
    pub width: u32,
    /// A description of the system and world (the `show` command's), shown as is in the full
    /// report.
    pub system_text: Option<String>,
    /// Every map and the system, not only the main maps ([`LITE`]).
    pub full: bool,
}

impl Default for Options {
    fn default() -> Self {
        Self { width: 1200, system_text: None, full: false }
    }
}

/// A report page and its images (`img/<name>.jpg` in the page).
#[derive(Debug, Clone)]
pub struct Report {
    pub title: String,
    pub html: String,
    pub images: Vec<(String, RgbImage)>,
}

/// The maps of the lite report: what the world looks like and how it works, in one each.
pub const LITE: [&str; 13] = ["appearance", "plates", "provinces", "lithology", "elevation", "wind_jan", "wind_jul", "ocean_currents", "temperature", "precipitation", "koppen", "rivers", "biome"];

/// The maps, by layer in the order the pipeline makes them, each derived view after the layer
/// it is read from: anchor, heading, and product or view names.
pub(crate) const GROUPS: [(&str, &str, &[&str]); 11] = [
    ("appearance", "From space", &["appearance"]),
    ("tectonics", "Tectonics", &["plates", "earthquakes", "volcanoes", "tsunamis"]),
    ("geology", "Geology", &["provinces", "lithology", "rock_age"]),
    ("elevation", "Elevation", &["rough_elevation"]),
    ("circulation", "Circulation", &["wind_jan", "wind_jul", "pressure_jan", "pressure_jul", "ocean_currents", "sst_anomaly"]),
    (
        "climate",
        "Climate",
        &["temperature", "temperature_range", "precipitation", "koppen", "sea_ice", "reefs", "productivity", "cyclones", "tornadoes", "blizzards", "dust_storms", "fog", "lightning"],
    ),
    ("erosion", "Erosion", &["elevation"]),
    ("hydrology", "Hydrology", &["rivers", "lakes", "floods", "dunes"]),
    ("biomes", "Biomes and soils", &["biome", "soil", "fertility", "timber", "wildfire"]),
    ("resources", "Resources", &["deposits", "iron", "copper", "gold", "coal", "salt"]),
    ("space", "Space weather", &["auroras"]),
];

pub(crate) fn caption(name: &str) -> String {
    let c = match name {
        "appearance" => "True colour from space, without clouds: seas by depth, sea ice, land by what covers it, snow",
        "rough_elevation" => "Elevation from the plates, before erosion",
        "elevation" => "Elevation after erosion",
        "plates" => "Tectonic plates",
        "provinces" => "Geologic provinces",
        "lithology" => "Surface rock",
        "rock_age" => "Age of surface rock",
        "temperature" => "Annual mean temperature",
        "temperature_range" => "Seasonal range, warmest to coldest month",
        "precipitation" => "Annual precipitation",
        "koppen" => "Köppen–Geiger climates",
        "wind_jan" => "Surface winds, January",
        "wind_jul" => "Surface winds, July",
        "ocean_currents" => "Ocean currents",
        "sst_anomaly" => "Sea-surface temperature anomaly",
        "pressure_jan" => "Sea-level pressure, January",
        "pressure_jul" => "Sea-level pressure, July",
        "sea_ice" => "Sea ice, seasonal and perennial",
        "reefs" => "Coral reefs and atolls",
        "productivity" => "Ocean productivity",
        "biome" => "Biomes",
        "rivers" => "Rivers by size",
        "lakes" => "Lakes",
        "dunes" => "Sand seas by their dunes' form",
        "soil" => "Soil groups",
        "fertility" => "Soil fertility",
        "timber" => "Timber",
        "deposits" => "Giant ore deposits",
        "iron" | "copper" | "gold" | "coal" | "salt" => return format!("{}: chance of a deposit per 10,000 km²", upper_first(name)),
        "earthquakes" => "Earthquakes M6+ in 10 years",
        "volcanoes" => "Eruptions in 10 years",
        "tsunamis" => "Damaging tsunamis in 100 years",
        "cyclones" => "Tropical cyclones in a year",
        "tornadoes" => "Tornadoes in a year",
        "floods" => "River floods in 10 years",
        "wildfire" => "Large wildfires in 10 years",
        "blizzards" => "Blizzard days",
        "dust_storms" => "Dust storm days",
        "fog" => "Days of thick fog",
        "lightning" => "Lightning",
        "auroras" => "Nights with auroras overhead",
        other => other,
    };
    c.to_string()
}

/// Build the report of `world` (its layers generated).
///
/// # Errors
/// If the world's layers are not generated yet, or on store errors.
pub fn build(world: &World, opts: &Options) -> Result<Report> {
    let products = world.load_products()?;
    let body = world.load_body()?;
    let grid = world.manifest.grid.build()?;
    let system = world.load_system()?;
    let (radius, sea) = (body.exterior.radius_km, body.interior.has_surface_liquid_water());
    let calendar = Calendar::of_body(&body);
    let summary = stats::summarize(&products, &grid, radius, sea);
    let stat = |sec: &str, k: &str| summary.get(sec).and_then(|s| s.get(k)).map(|s| s.value);

    // Maps, each with its key.
    let atm = geoforge_almanac::weather::atmosphere(world, &products, &grid, &body).ok();
    let mut source = views::Source::new(&products, &grid, &body);
    if let Some(a) = &atm {
        source = source.with_storms(a);
    }
    let sea_cells = source.sea();
    let mut images: Vec<(String, RgbImage)> = Vec::new();
    let mut groups_html = String::new();
    let mut nav = Vec::new();
    // The lite report's few maps run in one section in `LITE`'s order; the full report's by
    // layer.
    let mut lite_figs: Vec<(usize, String)> = Vec::new();
    for (anchor, heading, names) in GROUPS {
        let mut figs = String::new();
        for &name in names.iter().filter(|n| opts.full || LITE.contains(n)) {
            let Some((img, key)) = map(name, &products, &grid, &body, &source, &sea_cells)? else { continue };
            let fig = figure(name, &caption(name), key.as_ref());
            if let Some(k) = LITE.iter().position(|n| *n == name).filter(|_| !opts.full) {
                lite_figs.push((k, fig));
            } else {
                figs += &fig;
            }
            images.push((name.to_string(), fit(&img, opts.width)));
        }
        if !figs.is_empty() {
            nav.push((anchor, heading));
            let _ = write!(groups_html, r#"<section class="band" id="{anchor}"><div class="wrap"><h2>{}</h2><div class="plates">{figs}</div></div></section>"#, esc(heading));
        }
    }

    if !lite_figs.is_empty() {
        lite_figs.sort_by_key(|(k, _)| *k);
        let figs: String = lite_figs.into_iter().map(|(_, f)| f).collect();
        groups_html = format!(r#"<section class="band" id="maps"><div class="wrap"><h2>Maps</h2><div class="plates">{figs}</div></div></section>"#);
    }

    // The world whole, as the viewer's globe shows it: centred on the longitude with the
    // most land.
    let globe = map("appearance", &products, &grid, &body, &source, &sea_cells)?.map(|(flat, _)| {
        let mut land = [0usize; 36];
        for (c, sea) in sea_cells.iter().enumerate() {
            if !sea {
                let lon = grid.center(c as geoforge_core::grid::Cell).to_lat_lon().lon_deg;
                land[((lon + 180.0) / 10.0).clamp(0.0, 35.0) as usize] += 1;
            }
        }
        let best = (0..36).max_by_key(|&k| (0..9).map(|d| land[(k + d) % 36]).sum::<usize>()).unwrap_or(0);
        render::render_globe(&flat, opts.width.min(900), 18.0, best as f64 * 10.0 + 45.0 - 180.0, [15, 43, 51])
    });
    if let Some(g) = globe {
        images.push(("globe".into(), g));
    }

    // Cover and numbers.
    let star = system.as_ref().and_then(|s| geoforge_almanac::ephemeris::light_source(s, world.manifest.focus).and_then(|k| s.node(k))).and_then(|n| match &n.node {
        geoforge_core::body::NodeKind::Star(st) => Some(st.spectral_class.clone()),
        _ => None,
    });
    let title = format!("{}: {}", body.name, geoforge_cosmos::describe(&body));
    let mut figs = vec![
        (format!("{:.0} km", radius), format!("radius ({:.2} Earth)", radius / 6371.0)),
        (format!("{:.2} g", body.gravity_m_s2() / 9.81), "gravity".to_string()),
        (format!("{:.2} bar", body.interior.surface_pressure_bar()), "air at the surface".to_string()),
        (calendar.solar_day_h.map_or_else(|| "none".into(), |d| if d < 72.0 { format!("{d:.1} h") } else { format!("{:.1} days", d / 24.0) }), "solar day".to_string()),
        (format!("{:.0} days", calendar.year_days), star.map_or_else(|| "year".to_string(), |c| format!("year around a {c} star"))),
    ];
    if let Some(t) = stat("temperature", "global_mean_c") {
        figs.push((format!("{t:.1} °C"), "global mean (Earth 14.5)".into()));
    }
    if let Some(a) = system.as_ref().and_then(|s| geoforge_cosmos::asteroids::of(s, world.manifest.focus)) {
        figs.push((format!("{}-type", a.taxon.code()), a.leading(4)));
    }
    if let Some(l) = stat("elevation", "land_fraction") {
        figs.push((format!("{:.0}%", 100.0 * l), "land (Earth 29%)".into()));
    }
    let figs_html: String = figs.iter().map(|(a, b)| format!(r#"<div class="fig"><b>{}</b><span>{}</span></div>"#, esc(a), esc(b))).collect();
    let present = system.as_ref().and_then(|s| real_date(s, 0.0)).map_or(String::new(), |d| format!(" Its present is {d}."));
    let (w, h) = match world.manifest.grid {
        geoforge_core::grid::GridSpec::Equirect { width, height } => (width, height),
        _ => (0, 0),
    };
    let at = world.manifest.address.map_or_else(|| format!("Seed {}", world.manifest.root_seed.0), |a| a.to_string());
    let lede = format!("{at} · {w} × {h} grid, {:.0} km cells · geoforge {}.{present}", std::f64::consts::TAU * radius / f64::from(w.max(1)), world.manifest.code_version);
    let cover_img = if images.iter().any(|(n, _)| n == "globe") { format!(r#"<img class="globe" src="img/globe.jpg" alt="{} in true colour">"#, esc(&body.name)) } else { String::new() };
    let rows = stats_rows(&summary);

    // The full report adds its links and the system's bodies; the lite one is short enough
    // to scroll.
    let (nav_html, system_html) = if opts.full {
        let mut items: Vec<(&str, &str)> = vec![("system", "System")];
        items.extend(nav.iter().copied());
        items.push(("numbers", "Numbers"));
        let links: String = items.iter().map(|(a, t)| format!(r##"<a href="#{a}">{}</a>"##, esc(t))).collect();
        let text = opts.system_text.as_deref().map_or(String::new(), |t| format!("<pre>{}</pre>", esc(t)));
        (format!(r#"<nav class="nav">{links}</nav>"#), format!(r#"<section class="band" id="system"><div class="wrap"><h2>The system</h2>{text}</div></section>"#))
    } else {
        (String::new(), String::new())
    };
    let open = if opts.full { " open" } else { "" };
    let html = format!(
        r#"<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{t}</title>
<link rel="preconnect" href="https://fonts.googleapis.com"><link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;600&family=Spectral:wght@500&display=swap" rel="stylesheet">
<style>{CSS}</style></head><body>
<header class="band cover"><div class="wrap"><div><h1>{t}</h1><p class="lede">{lede}</p><div class="figs">{figs_html}</div></div>{cover_img}</div></header>
{nav_html}
{system_html}
{groups_html}
<section class="band" id="numbers"><div class="wrap"><details{open}><summary><h2>{name} against Earth's numbers</h2></summary>
<p class="note">Green is within about a third of Earth's figure (temperatures within 2.5 K), coral further off: where the world differs, not necessarily where it is wrong.</p>
<div class="scroll"><table class="stack"><thead><tr><th>Layer</th><th>Measure</th><th>{name}</th><th>Earth</th></tr></thead><tbody>{rows}</tbody></table></div></details></div></section>
<script>{SCRIPT}</script></body></html>"#,
        t = esc(&title),
        name = esc(&body.name),
    );
    Ok(Report { title, html, images })
}

/// A map by name (a product or a view) and its key; `None` if the world lacks it.
pub(crate) fn map(name: &str, products: &Products, grid: &AnyGrid, body: &BodyContext, source: &views::Source, sea_cells: &[bool]) -> Result<Option<(RgbImage, Option<Legend>)>> {
    let (radius, sea) = (body.exterior.radius_km, body.interior.has_surface_liquid_water());
    if let Some(p) = ProductId::from_name(name) {
        return Ok(render::render(p, products, grid, radius, sea, body.exterior.class)?.map(|img| (img, render::legend(p, products, sea))));
    }
    let Some(v) = View::from_name(name) else { return Err(PipelineError::Invalid(format!("no product or view `{name}`"))) };
    let relief = Some((products, radius));
    Ok(match views::compute(v, source) {
        Some(views::Map::Classes(c)) => {
            let key = Legend::Classes(c.legend.iter().filter(|k| k.area_share > 0.0).map(|k| (format!("{} {}", k.code, k.name).trim().to_string(), k.color)).collect());
            Some((render::render_classes(&c.colors(), grid, relief)?, Some(key)))
        }
        Some(views::Map::Chances(c)) => Some((render::render_chances(&c.cells, sea_cells, grid, relief)?, Some(render::chance_legend()))),
        Some(views::Map::Quantity(q)) => {
            let span = 6.0 * q.contour;
            Some((render::render_quantity(&q.cells, sea_cells, q.contour, span, grid, relief)?, Some(render::quantity_legend(&q.cells, span, q.units))))
        }
        Some(views::Map::Colours(c)) => Some((render::render_classes(&c.into_iter().map(Some).collect::<Vec<_>>(), grid, relief)?, None)),
        None => None,
    })
}

/// `img` scaled to `width` (half as tall).
pub(crate) fn fit(img: &RgbImage, width: u32) -> RgbImage {
    if img.width() == width {
        return img.clone();
    }
    image::imageops::resize(img, width, (width * img.height() / img.width().max(1)).max(1), image::imageops::FilterType::Triangle)
}

fn figure(name: &str, caption: &str, key: Option<&Legend>) -> String {
    format!(r#"<figure><img loading="lazy" src="img/{name}.jpg" alt="{c}"><figcaption><b>{c}</b>{}</figcaption></figure>"#, key.map_or(String::new(), key_html), c = esc(caption))
}

fn hex(c: [u8; 3]) -> String {
    format!("#{:02x}{:02x}{:02x}", c[0], c[1], c[2])
}

/// A key: named swatches, or a ramp's colours under their values.
pub(crate) fn key_html(key: &Legend) -> String {
    match key {
        Legend::Classes(items) => {
            let li: String = items.iter().map(|(n, c)| format!(r#"<li><i style="background:{}"></i>{}</li>"#, hex(*c), esc(n))).collect();
            if items.len() > 30 {
                format!(r#"<details><summary>Key: {} classes</summary><ul class="key">{li}</ul></details>"#, items.len())
            } else {
                format!(r#"<ul class="key">{li}</ul>"#)
            }
        }
        Legend::Ramp { stops, units } => {
            let chance = *units == "chance";
            // Every other value when they would crowd a phone's width.
            let thin = stops.len() > 8;
            let label = |i: usize, v: f64| match () {
                () if thin && i % 2 == 1 => String::new(),
                () if chance => format!("{:.0}%", 100.0 * v),
                () => short(v),
            };
            let cells: String = stops.iter().enumerate().map(|(i, (v, c))| format!(r#"<span><i style="background:{}"></i>{}</span>"#, hex(*c), label(i, *v))).collect();
            let units = if chance { "chance" } else { units };
            format!(r#"<div class="ramp">{cells}</div><div class="units">{}</div>"#, esc(units))
        }
    }
}

pub(crate) fn num(x: f64) -> String {
    if x.abs() >= 1000.0 {
        let s = format!("{:.0}", x.abs());
        let mut out = String::new();
        for (i, ch) in s.chars().enumerate() {
            if i > 0 && (s.len() - i) % 3 == 0 {
                out.push(',');
            }
            out.push(ch);
        }
        format!("{}{out}", if x < 0.0 { "−" } else { "" })
    } else if x.abs() >= 10.0 || x == x.round() {
        format!("{x:.0}").replace('-', "−")
    } else if x.abs() >= 1.0 {
        format!("{x:.1}").replace('-', "−")
    } else {
        format!("{x:.2}").replace('-', "−")
    }
}

/// A ramp label: thousands as "4.5k".
fn short(x: f64) -> String {
    if x.abs() >= 1000.0 {
        let k = x / 1000.0;
        let s = if (k - k.round()).abs() < 0.05 { format!("{k:.0}k") } else { format!("{k:.1}k") };
        s.replace('-', "−")
    } else {
        num(x)
    }
}

fn stats_rows(summary: &stats::Summary) -> String {
    let order = ["temperature", "precipitation", "elevation", "hydrology", "biomes", "plates", "crust"];
    let mut secs: Vec<&String> = summary.keys().collect();
    secs.sort_by_key(|s| order.iter().position(|o| o == s).unwrap_or(order.len()));
    let mut rows = String::new();
    for sec in secs {
        for (k, v) in &summary[sec] {
            let Some(earth) = v.earth else { continue };
            let close = (earth != 0.0 && (0.75..=1.33).contains(&(v.value / earth))) || (sec == "temperature" && (v.value - earth).abs() < 2.5);
            let label = k.replace('_', " ").replace(" of land", " (share of land)");
            let _ = write!(rows, r#"<tr><td>{}</td><td>{}</td><td class="n {}">{}</td><td class="n">{}</td></tr>"#, esc(sec), esc(&label), if close { "ok" } else { "miss" }, num(v.value), num(earth));
        }
    }
    rows
}

pub(crate) fn lat_lon(lat: f64, lon: f64) -> String {
    format!("{:.0}°{} {:.0}°{}", lat.abs(), if lat >= 0.0 { "N" } else { "S" }, lon.abs(), if lon >= 0.0 { "E" } else { "W" })
}

fn upper_first(s: &str) -> String {
    let mut c = s.chars();
    c.next().map_or_else(String::new, |f| f.to_uppercase().collect::<String>() + c.as_str())
}

pub(crate) fn esc(s: &str) -> String {
    s.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;").replace('"', "&quot;")
}

/// The page's style: a field atlas, deep sea and chart paper with one coral signal.
pub(crate) const CSS: &str = r#":root{color-scheme:light dark;--bg:light-dark(#f4f6f3,#0c1d22);--ink:light-dark(#10232a,#e4ece9);--soft:light-dark(#5a6b70,#9fb2b4);--line:light-dark(#d5ddda,#21383f);--card:light-dark(#ffffff,#112830);--sea:#0f2b33;--hot:#e0533a;--ok:light-dark(#2f7d5b,#5fc297)}
*{box-sizing:border-box;min-width:0}
body{margin:0;background:var(--bg);color:var(--ink);font:16px/1.55 Inter,system-ui,sans-serif;overflow-x:clip}
h1,h2,h3{font-family:Spectral,Georgia,serif;font-weight:500;line-height:1.08;text-wrap:balance;margin:0 0 .4em}
h1{font-size:clamp(30px,5vw,48px)}h2{font-size:clamp(24px,3vw,32px);margin-top:.2em}h3{font-size:19px;margin-top:1em}
.band{padding:clamp(28px,4vw,56px) max(5vw,16px)}.wrap{max-width:1240px;margin:0 auto}
.cover{background:var(--sea);color:#e8f1ef}.cover .lede{color:#b9cdcc;max-width:68ch;font-size:13px;overflow-wrap:anywhere}
.cover .wrap{display:grid;grid-template-columns:minmax(0,1.1fr) minmax(0,1fr);gap:clamp(20px,4vw,56px);align-items:center}
.cover img.globe{width:100%;max-width:560px;display:block;justify-self:center;aspect-ratio:1;border-radius:50%}
@media(max-width:760px){.cover .wrap{grid-template-columns:minmax(0,1fr)}.cover img.globe{max-width:420px}}
.figs{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,150px),1fr));gap:10px;margin-top:22px}
.fig{background:#ffffff10;border:1px solid #ffffff22;border-radius:6px;padding:10px 12px}
.fig b{display:block;font:500 clamp(20px,2.4vw,26px)/1.1 Spectral,serif;color:#fff}.fig span{font-size:13px;color:#b9cdcc}
.nav{display:flex;gap:6px;overflow-x:auto;padding:10px max(5vw,16px);border-bottom:1px solid var(--line);position:sticky;top:0;background:var(--bg);z-index:2}
.nav a{flex:none;color:var(--ink);text-decoration:none;font-size:14px;padding:8px 12px;border-radius:999px;border:1px solid var(--line)}
.nav a:hover{border-color:var(--hot)}
.plates{display:grid;grid-template-columns:repeat(auto-fill,minmax(min(100%,380px),1fr));gap:18px;margin-top:14px}
figure{margin:0;background:var(--card);border:1px solid var(--line);border-radius:6px;overflow:hidden}
figure img{width:100%;display:block;aspect-ratio:2/1;object-fit:cover;background:#0f2b33}
figcaption{padding:8px 12px 12px;font-size:14px;color:var(--soft)}figcaption b{color:var(--ink);font-weight:600}
ul.key{list-style:none;padding:0;margin:8px 0 0;display:flex;flex-wrap:wrap;gap:4px 12px;font-size:13px}
ul.key li{display:flex;align-items:center;gap:5px}ul.key i,.ramp i{display:inline-block;width:14px;height:14px;border-radius:3px;border:1px solid #0002;flex:none}
.ramp{display:flex;margin-top:8px}.ramp span{flex:1;display:flex;flex-direction:column;align-items:center;font-size:12px;font-variant-numeric:tabular-nums}.ramp span i{width:100%;height:12px;border-radius:0;border:0}
.units{font-size:12px;text-align:center;margin-top:2px}
details summary{cursor:pointer;font-size:13px;min-height:32px;display:flex;align-items:center}
details summary h2{margin:0;font-size:inherit}#numbers summary{font-size:clamp(24px,3vw,32px);font-family:Spectral,Georgia,serif;list-style-position:inside}
p.note{color:var(--soft);max-width:75ch}
table{border-collapse:collapse;width:100%;font-size:14.5px;margin-top:10px}
th,td{padding:7px 10px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}
th{font-weight:600;font-size:13px;color:var(--soft)}td.n{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}
.ok{color:var(--ok)}.miss{color:var(--hot)}
.cols{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,340px),1fr));gap:22px;margin-top:18px}
.card{background:var(--card);border:1px solid var(--line);border-radius:6px;padding:16px 18px}
pre{background:var(--card);border:1px solid var(--line);border-radius:6px;padding:12px;font-size:12.5px;line-height:1.45;overflow-x:auto}
.scroll{overflow-x:auto}
@media(max-width:600px){table{font-size:13px}.stack thead{display:none}.stack tr{display:grid;padding-block:10px}.stack td{display:flow-root;border:0;padding:3px 0;padding-inline-start:8em;text-align:start}.stack td::before{content:attr(data-l);float:inline-start;width:7.5em;margin-inline-start:-8em;white-space:normal;color:var(--soft)}.stack td:empty{display:none}}"#;

/// Labels each cell of a stacked table with its column's heading, for narrow screens.
pub(crate) const SCRIPT: &str =
    "for(const c of document.querySelectorAll('.stack td')){let i=0,p=c;while(p=p.previousElementSibling)i+=p.colSpan;c.dataset.l=c.closest('table').rows[0].cells[i]?.textContent||''}";
