//! A report on a world, to take its pulse at a glance: its maps by layer in pipeline order, each with
//! a key; its numbers against Earth's; three places with their weather; the weather map and
//! the night sky with keys.
//!
//! [`build`] makes the page and its images; `geoforge_store::report::write` puts them in a
//! folder (`index.html`, `img/`).

use crate::{place, stats, views, PipelineError, Result, View};
use geoforge_almanac::ephemeris::{real_date, Calendar, EclipseKind};
use geoforge_almanac::sky::Kind;
use geoforge_core::grid::{AnyGrid, Cell, Grid};
use geoforge_core::product::{ProductId, Products};
use geoforge_core::BodyContext;
use geoforge_store::render::{self, Legend, SkyPoint};
use geoforge_store::World;
use image::RgbImage;
use std::fmt::Write;

/// A place to report on.
#[derive(Debug, Clone, PartialEq)]
pub struct Place {
    pub name: String,
    pub lat: f64,
    pub lon: f64,
}

/// What goes into a report.
#[derive(Debug, Clone, PartialEq)]
pub struct Options {
    /// Width of the maps and charts, px.
    pub width: u32,
    /// Places to report on; empty picks a temperate, a tropical and a dry one.
    pub places: Vec<Place>,
    /// The day (Earth days from the present) of the weather, weather map and sky.
    pub day: f64,
    /// Days of weather at each place.
    pub days: usize,
    /// A description of the system and world (the `show` command's), shown as is.
    pub system_text: Option<String>,
}

impl Default for Options {
    fn default() -> Self {
        Self { width: 1200, places: Vec::new(), day: 30.0, days: 14, system_text: None }
    }
}

/// A report page and its images (`img/<name>.jpg` in the page).
#[derive(Debug, Clone)]
pub struct Report {
    pub title: String,
    pub html: String,
    pub images: Vec<(String, RgbImage)>,
}

/// The maps, by layer in the order the pipeline makes them, each derived view after the layer
/// it is read from: anchor, heading, and product or view names.
pub(crate) const GROUPS: [(&str, &str, &[&str]); 10] = [
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
    ("hydrology", "Hydrology", &["rivers", "lakes", "floods"]),
    ("biomes", "Biomes and soils", &["biome", "soil", "fertility", "timber", "wildfire"]),
    ("resources", "Resources", &["deposits", "iron", "copper", "gold", "coal", "salt"]),
    ("space", "Space weather", &["auroras"]),
];

pub(crate) fn caption(name: &str) -> String {
    let c = match name {
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
    for (anchor, heading, names) in GROUPS {
        let mut figs = String::new();
        for &name in names {
            let Some((img, key)) = map(name, &products, &grid, &body, &source, &sea_cells)? else { continue };
            figs += &figure(name, &caption(name), key.as_ref());
            images.push((name.to_string(), fit(&img, opts.width)));
        }
        if !figs.is_empty() {
            nav.push((anchor, heading));
            let _ = write!(groups_html, r#"<section class="band" id="{anchor}"><div class="wrap"><h2>{}</h2><div class="plates">{figs}</div></div></section>"#, esc(heading));
        }
    }

    // Places and their weather.
    let places = if opts.places.is_empty() { pick_places(&products, &grid, &source) } else { opts.places.clone() };
    let mut cards = String::new();
    let mut weather = String::new();
    let dates = geoforge_almanac::calendar::of_world(world, &geoforge_almanac::calendar::Choice::World).ok();
    for p in &places {
        if let Ok(r) = place::report(world, p.lat, p.lon, 100.0, 0.0, calendar.year_days) {
            cards += &place_card(&p.name, &r);
        }
        if let Ok(f) = geoforge_almanac::weather::forecast(world, p.lat, p.lon, opts.day, opts.days) {
            weather += &weather_table(&format!("{} from day {}", p.name, opts.day), &f, system.as_ref(), dates.as_deref());
        }
    }
    let wmap = geoforge_almanac::weather::weather_map_image(world, opts.day, opts.width).ok().map(|(img, storms)| {
        images.push(("weather_map".into(), img));
        storms.len()
    });

    // The sky.
    let sky_html = match &system {
        Some(sys) => match geoforge_almanac::sky::starfield(world)
            .map_err(|e| e.to_string())
            .and_then(|field| geoforge_almanac::sky::sky(sys, world.manifest.focus, opts.day, calendar.year_days, &field).map(|s| (s, field)))
        {
            Ok((s, field)) => {
                let (gw, gh) = (360, 180);
                let lines: Vec<[(f64, f64); 2]> = s.figures(&world.load_constellations().unwrap_or_default(), &field).iter().flat_map(|f| f.lines.iter().copied()).collect();
                let points: Vec<SkyPoint> = s.stars.iter().chain(s.bodies.iter().filter(|b| b.v <= geoforge_cosmos::sky::NAKED_EYE)).map(sky_point).collect();
                images.push(("sky".into(), render::render_sky(opts.width, &s.glow(gw, gh), gw, &points, &s.ecliptic, &lines)));
                let key = |points: &[SkyPoint]| render::render_sky_key(56, KEY_SCALE, points);
                images.push(("sky_key_magnitude".into(), key(&MAGNITUDES.map(|v| key_point(v, 5800.0, false)))));
                images.push(("sky_key_colour".into(), key(&SPECTRAL.map(|(_, t)| key_point(1.0, t, false)))));
                images.push(("sky_key_planet".into(), key(&[key_point(0.0, 5800.0, true)])));
                sky_section(&body, sys, &s, opts.day)
            }
            Err(e) => format!(r#"<p class="note">No sky: {}</p>"#, esc(&e)),
        },
        None => r#"<p class="note">This world has no solar system (a hand-made body), so it has no sky to draw.</p>"#.into(),
    };

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
    if let Some(l) = stat("elevation", "land_fraction") {
        figs.push((format!("{:.0}%", 100.0 * l), "land (Earth 29%)".into()));
    }
    let figs_html: String = figs.iter().map(|(a, b)| format!(r#"<div class="fig"><b>{}</b><span>{}</span></div>"#, esc(a), esc(b))).collect();
    let present = system.as_ref().and_then(|s| real_date(s, 0.0)).map_or(String::new(), |d| format!(" Its present is {d}."));
    let (w, h) = match world.manifest.grid {
        geoforge_core::grid::GridSpec::Equirect { width, height } => (width, height),
        _ => (0, 0),
    };
    let lede = format!(
        "Seed {}, a {w} × {h} grid ({:.0} km cells at the equator), every layer generated by geoforge {}.{present} The maps below come straight from the world's layers and derived views; each has its key.",
        world.manifest.root_seed.0,
        std::f64::consts::TAU * radius / f64::from(w.max(1)),
        world.manifest.code_version
    );
    let cover_img = ["biome", "elevation"].into_iter().find(|n| images.iter().any(|(m, _)| m == n)).map_or(String::new(), |n| format!(r#"<img src="img/{n}.jpg" alt="{}">"#, esc(&caption(n))));
    let rows = stats_rows(&summary);

    let mut nav_items: Vec<(&str, &str)> = vec![("system", "System")];
    nav_items.extend(nav.iter().copied());
    nav_items.extend([("numbers", "Numbers"), ("places", "Places"), ("weather", "Weather"), ("sky", "Sky")]);
    let nav_html: String = nav_items.iter().map(|(a, t)| format!(r##"<a href="#{a}">{}</a>"##, esc(t))).collect();
    let system_text = opts.system_text.as_deref().map_or(String::new(), |t| format!("<pre>{}</pre>", esc(t)));
    let wmap_html = wmap.map_or(String::new(), |n| {
        format!(
            r#"<figure class="sq"><img src="img/weather_map.jpg" alt="Weather map"><figcaption><b>Weather map, day {}{}.</b> Isobars every 4 hPa, rain in green, warm fronts red, cold fronts blue, low centres as dots; {n} storms.</figcaption></figure>"#,
            opts.day,
            system.as_ref().and_then(|s| real_date(s, opts.day)).map_or(String::new(), |d| format!(" ({d})"))
        )
    });
    let html = format!(
        r#"<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{t}</title>
<link rel="preconnect" href="https://fonts.googleapis.com"><link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;600&family=Spectral:wght@500&display=swap" rel="stylesheet">
<style>{CSS}</style></head><body>
<header class="band cover"><div class="wrap"><h1>{t}</h1><p class="lede">{lede}</p><div class="figs">{figs_html}</div>{cover_img}</div></header>
<nav class="nav">{nav_html}</nav>
<section class="band" id="system"><div class="wrap"><h2>The system</h2>{system_text}</div></section>
{groups_html}
<section class="band" id="numbers"><div class="wrap"><h2>{name} against Earth's numbers</h2>
<p class="note">Green is within about a third of Earth's figure (temperatures within 2.5 K), coral further off: where the world differs, not necessarily where it is wrong.</p>
<div class="scroll"><table class="stack"><thead><tr><th>Layer</th><th>Measure</th><th>{name}</th><th>Earth</th></tr></thead><tbody>{rows}</tbody></table></div></div></section>
<section class="band" id="places"><div class="wrap"><h2>Places</h2><div class="cols">{cards}</div></div></section>
<section class="band" id="weather"><div class="wrap"><h2>Weather from day {day}</h2>{wmap_html}{weather}</div></section>
<section class="band" id="sky"><div class="wrap"><h2>The night sky, day {day}</h2>{sky_html}</div></section>
<script>{SCRIPT}</script></body></html>"#,
        t = esc(&title),
        name = esc(&body.name),
        day = opts.day,
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

fn place_card(name: &str, r: &place::PlaceReport) -> String {
    let mut deposits: Vec<_> = r.deposits.iter().collect();
    deposits.sort_by(|a, b| b.chance.total_cmp(&a.chance));
    let dep = deposits.iter().take(5).map(|d| format!("{} {:.0}%", d.kind.replace('_', " "), 100.0 * d.chance)).collect::<Vec<_>>().join(", ");
    let mut hazards: Vec<_> = r.hazards.iter().collect();
    hazards.sort_by(|a, b| b.chance.total_cmp(&a.chance));
    let hz: String = hazards
        .iter()
        .filter(|h| h.chance >= 0.005)
        .take(5)
        .map(|h| format!(r#"<tr><td>{}</td><td class="n">{:.3} a year</td><td class="n">{:.0}%</td></tr>"#, esc(h.event), h.per_year, 100.0 * h.chance))
        .collect();
    format!(
        r#"<div class="card"><h3>{}</h3><p class="note">{:.1}°, {:.1}° · {:.0} m · Köppen <b>{}</b> · {} (fertility {:.1}) · timber {:.1}</p><p><b>Deposits within 100 km</b>, chance of at least one: {}</p>{}</div>"#,
        esc(name),
        r.lat,
        r.lon,
        r.elevation_m,
        r.koppen.unwrap_or("—"),
        r.soil.unwrap_or("—"),
        r.fertility.unwrap_or(0.0),
        r.timber.unwrap_or(0.0),
        if dep.is_empty() { "none likely".into() } else { esc(&dep) },
        if hz.is_empty() {
            r#"<p class="note">No hazard with even a 1 in 200 chance a year within 100 km.</p>"#.to_string()
        } else {
            format!("<table><thead><tr><th>Most likely hazards within 100 km</th><th>Rate</th><th>In a year</th></tr></thead><tbody>{hz}</tbody></table>")
        }
    )
}

fn weather_table(title: &str, f: &geoforge_almanac::weather::Forecast, system: Option<&geoforge_core::body::SolarSystem>, dates: Option<&dyn geoforge_core::time::Calendar>) -> String {
    let rows: String = f
        .days
        .iter()
        .map(|d| {
            let date = system.and_then(|s| real_date(s, d.day)).map_or_else(
                || dates.map_or_else(|| format!("{:.0}", d.local_day.floor()), |c| c.format_date(geoforge_core::time::Instant::from_earth_days(d.day + 0.5 * f.step_days))),
                |s| s[..10].to_string(),
            );
            let ground = match (f.at_sea, d.snowpack_mm >= 1.0) {
                (true, _) => String::new(),
                (false, true) => format!("{:.0}%, {:.0} mm snow", 100.0 * d.soil, d.snowpack_mm),
                (false, false) => format!("{:.0}%", 100.0 * d.soil),
            };
            let precip = if d.precipitation_mm > 0.0 { format!("{:.1} mm{}", d.precipitation_mm, if d.snow { " snow" } else { "" }) } else { "—".into() };
            let mut sky = d.storm.as_ref().map_or_else(|| d.sky.to_string(), |name| format!("{}, {name} passing", d.sky));
            for e in d.events.iter().filter(|e| e.kind != "thunderstorm" || d.sky != "thunderstorm") {
                sky += &format!("; {}", e.what);
            }
            format!(
                r#"<tr><td class="n">{date}</td><td class="n">{:.0}° / {:.0}°</td><td class="n">{precip}</td><td class="n">{:.0} m/s</td><td class="n">{ground}</td><td>{}</td></tr>"#,
                d.low_c,
                d.high_c,
                d.wind_m_s,
                esc(&sky)
            )
        })
        .collect();
    let day = if system.and_then(|s| s.root.present_jd_tt).is_some() || dates.is_some() { "Date" } else { "Local day" };
    format!(
        r#"<h3>{}</h3><div class="scroll"><table class="wx"><thead><tr><th>{day}</th><th>Low / high</th><th>Rain</th><th>Wind</th><th>Soil water</th><th>Sky</th></tr></thead><tbody>{rows}</tbody></table></div>"#,
        esc(title)
    )
}

/// The sky's key draws its stars this much larger than the chart does, to be read on a phone.
const KEY_SCALE: f64 = 2.5;
/// Magnitudes and spectral classes shown in the sky's key.
const MAGNITUDES: [f64; 8] = [-1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
const SPECTRAL: [(&str, f64); 7] = [("O", 35_000.0), ("B", 15_000.0), ("A", 9_000.0), ("F", 6_800.0), ("G", 5_800.0), ("K", 4_500.0), ("M", 3_200.0)];

fn key_point(v: f64, temperature_k: f64, ring: bool) -> SkyPoint {
    SkyPoint { ra_deg: 0.0, dec_deg: 0.0, v, temperature_k, ring }
}

fn sky_point(o: &geoforge_almanac::sky::Object) -> SkyPoint {
    SkyPoint { ra_deg: o.ra_deg, dec_deg: o.dec_deg, v: o.v, temperature_k: o.temperature_k, ring: matches!(o.kind, Kind::Planet | Kind::Moon) }
}

fn sky_section(body: &BodyContext, system: &geoforge_core::body::SolarSystem, s: &geoforge_almanac::sky::Sky, day: f64) -> String {
    let labels = |items: Vec<String>| -> String { items.into_iter().map(|l| format!("<span>{}</span>", esc(&l))).collect() };
    let mags = labels(MAGNITUDES.iter().map(|v| num(*v)).collect());
    let classes = labels(SPECTRAL.iter().map(|(c, t)| format!("{c} {}", short(*t))).collect());
    let star_rows: String = s
        .stars
        .iter()
        .take(12)
        .map(|o| {
            format!(
                r#"<tr><td>{}</td><td class="n">{:.2}</td><td>{}</td><td class="n">{:.0} ly</td><td class="n">{:.1}° / {:+.1}°</td></tr>"#,
                esc(o.name.as_deref().unwrap_or("")),
                o.v,
                esc(o.spectral_class.as_deref().unwrap_or("")),
                o.distance,
                o.ra_deg,
                o.dec_deg
            )
        })
        .collect();
    let body_rows: String = s
        .bodies
        .iter()
        .map(|o| {
            let seen = match o.kind {
                Kind::Sun => "the star".to_string(),
                _ if o.v > geoforge_cosmos::sky::NAKED_EYE => "too faint to see".into(),
                _ if !o.visible() => "lost in the glare".into(),
                _ => "in view at night".into(),
            };
            let lit = o.lit.map_or(String::new(), |l| format!("{:.0}%", 100.0 * l));
            let elong = o.elongation_deg.map_or(String::new(), |e| format!("{e:.0}°"));
            format!(r#"<tr><td>{}</td><td class="n">{:.1}</td><td class="n">{elong}</td><td class="n">{lit}</td><td>{seen}</td></tr>"#, esc(o.name.as_deref().unwrap_or("")), o.v)
        })
        .collect();
    let eclipses: String = s
        .eclipses
        .iter()
        .map(|e| {
            let when = real_date(system, e.day).unwrap_or_else(|| format!("day {:.1}", e.day));
            let what = if e.kind == EclipseKind::Solar {
                format!("{:?} solar eclipse of {} by {}", e.coverage, e.star, e.body)
            } else {
                format!("{:?} lunar eclipse of {} by {}", e.coverage, e.body, body.name)
            };
            format!(r#"<li>{}: {} (magnitude {:.2})</li>"#, esc(&when), esc(&what), e.magnitude)
        })
        .collect();
    let eclipses = if eclipses.is_empty() { r#"<p class="note">No eclipses in the year.</p>"#.to_string() } else { format!("<h3>Eclipses in the year</h3><ul>{eclipses}</ul>") };
    let date = real_date(system, day).map_or(String::new(), |d| format!(" ({d})"));
    format!(
        r#"<figure class="sq"><img src="img/sky.jpg" alt="Star chart"><figcaption><b>Star chart of {name}'s sky, day {day}{date}:</b> right ascension and declination on the world's equator, RA 360° at the left to 0° at the right, as the sky looks from below. {n} naked-eye stars, drawn from the {galaxy} model {r:.1} kpc from its centre (not a catalogue).</figcaption></figure>
<div class="cols skykey"><div class="card"><h3>Key</h3>
<p class="note">Brightness by visual magnitude (a G star; drawn 2.5 times the chart's size):</p><img class="keystrip" src="img/sky_key_magnitude.jpg" alt="Star sizes by magnitude"><div class="keylabels">{mags}</div>
<p class="note">Colour by spectral class and surface temperature in K (magnitude 1):</p><img class="keystrip" src="img/sky_key_colour.jpg" alt="Star colours by class"><div class="keylabels">{classes}</div>
<p class="note"><img class="keydot" src="img/sky_key_planet.jpg" alt="Planet">A green ring marks a planet or moon of the system. The dotted gold line is the star's yearly path (the ecliptic); grid lines every 30°, the equator brighter; the blue-grey haze is the galaxy's band, its dust lane dark.</p></div>
<div class="card"><h3>Brightest stars</h3><div class="scroll"><table><thead><tr><th>Star</th><th>Mag</th><th>Class</th><th>Distance</th><th>RA / Dec</th></tr></thead><tbody>{star_rows}</tbody></table></div></div></div>
<h3>The system tonight</h3><div class="scroll"><table class="stack"><thead><tr><th>Body</th><th>Mag</th><th>From the star</th><th>Lit</th><th>Seen</th></tr></thead><tbody>{body_rows}</tbody></table></div>
{eclipses}"#,
        name = esc(&body.name),
        n = s.stars.len(),
        galaxy = esc(&s.galaxy),
        r = s.radius_kpc,
    )
}

/// A temperate, a tropical and a dry place, by Köppen class: of the land in the world's
/// commonest class of each kind, the cell with the median rainfall.
fn pick_places(products: &Products, grid: &AnyGrid, source: &views::Source) -> Vec<Place> {
    let (Some(views::Map::Classes(k)), Some(rain)) = (views::compute(View::Koppen, source), products.precipitation()) else { return Vec::new() };
    let mut out = Vec::new();
    for (label, prefixes) in [("Temperate", ["Cf", "C", "D"]), ("Tropical", ["Af", "A", "A"]), ("Dry", ["BWh", "BW", "B"])] {
        // The first prefix any land has, then that prefix's commonest class.
        let Some(class) = prefixes.iter().find_map(|pre| {
            let mut counts = std::collections::BTreeMap::new();
            for c in k.cells.iter().flatten().filter(|&&c| k.legend[c].code.starts_with(pre)) {
                *counts.entry(*c).or_insert(0usize) += 1;
            }
            counts.into_iter().max_by_key(|&(c, n)| (n, std::cmp::Reverse(c))).map(|(c, _)| c)
        }) else {
            continue;
        };
        let mut cells: Vec<Cell> = grid.cells().filter(|&c| k.cells[c as usize] == Some(class)).collect();
        cells.sort_by(|a, b| rain.values[*a as usize].total_cmp(&rain.values[*b as usize]).then(a.cmp(b)));
        let Some(&c) = cells.get(cells.len() / 2) else { continue };
        let ll = grid.center(c).to_lat_lon();
        let code = &k.legend[class].code;
        out.push(Place { name: format!("{label} ({code}), {}", lat_lon(ll.lat_deg, ll.lon_deg)), lat: ll.lat_deg, lon: ll.lon_deg });
    }
    out
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
.cover{background:var(--sea);color:#e8f1ef}.cover .lede{color:#b9cdcc;max-width:68ch;font-size:17px}
.cover img{width:100%;display:block;border-radius:6px;margin-top:22px;box-shadow:0 10px 40px #0007}
.figs{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,150px),1fr));gap:10px;margin-top:22px}
.fig{background:#ffffff10;border:1px solid #ffffff22;border-radius:6px;padding:10px 12px}
.fig b{display:block;font:500 clamp(20px,2.4vw,26px)/1.1 Spectral,serif;color:#fff}.fig span{font-size:13px;color:#b9cdcc}
.nav{display:flex;gap:6px;overflow-x:auto;padding:10px max(5vw,16px);border-bottom:1px solid var(--line);position:sticky;top:0;background:var(--bg);z-index:2}
.nav a{flex:none;color:var(--ink);text-decoration:none;font-size:14px;padding:8px 12px;border-radius:999px;border:1px solid var(--line)}
.nav a:hover{border-color:var(--hot)}
.plates{display:grid;grid-template-columns:repeat(auto-fill,minmax(min(100%,380px),1fr));gap:18px;margin-top:14px}
figure{margin:0;background:var(--card);border:1px solid var(--line);border-radius:6px;overflow:hidden}
figure img{width:100%;display:block;aspect-ratio:2/1;object-fit:cover;background:#0f2b33}
figure.sq img{aspect-ratio:auto}
figcaption{padding:8px 12px 12px;font-size:14px;color:var(--soft)}figcaption b{color:var(--ink);font-weight:600}
ul.key{list-style:none;padding:0;margin:8px 0 0;display:flex;flex-wrap:wrap;gap:4px 12px;font-size:13px}
ul.key li{display:flex;align-items:center;gap:5px}ul.key i,.ramp i{display:inline-block;width:14px;height:14px;border-radius:3px;border:1px solid #0002;flex:none}
.ramp{display:flex;margin-top:8px}.ramp span{flex:1;display:flex;flex-direction:column;align-items:center;font-size:12px;font-variant-numeric:tabular-nums}.ramp span i{width:100%;height:12px;border-radius:0;border:0}
.units{font-size:12px;text-align:center;margin-top:2px}
details summary{cursor:pointer;font-size:13px;min-height:32px;display:flex;align-items:center}
p.note{color:var(--soft);max-width:75ch}
table{border-collapse:collapse;width:100%;font-size:14.5px;margin-top:10px}
th,td{padding:7px 10px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}
th{font-weight:600;font-size:13px;color:var(--soft)}td.n{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}
.ok{color:var(--ok)}.miss{color:var(--hot)}
.cols{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,340px),1fr));gap:22px;margin-top:18px}
.card{background:var(--card);border:1px solid var(--line);border-radius:6px;padding:16px 18px}
img.keystrip{width:100%;max-width:448px;display:block;border-radius:4px}
.keylabels{display:flex;max-width:448px;font-size:12px}.keylabels span{flex:1;text-align:center}
table.wx td{padding:5px 8px}table.wx{font-size:13.5px}
img.keydot{width:40px;height:40px;vertical-align:middle;margin-right:8px;border-radius:4px}
pre{background:var(--card);border:1px solid var(--line);border-radius:6px;padding:12px;font-size:12.5px;line-height:1.45;overflow-x:auto}
.scroll{overflow-x:auto}
@media(max-width:600px){table{font-size:13px}.stack thead{display:none}.stack tr{display:grid;padding-block:10px}.stack td{display:flow-root;border:0;padding:3px 0;padding-inline-start:8em;text-align:start}.stack td::before{content:attr(data-l);float:inline-start;width:7.5em;margin-inline-start:-8em;white-space:normal;color:var(--soft)}.stack td:empty{display:none}}"#;

/// Labels each cell of a stacked table with its column's heading, for narrow screens.
pub(crate) const SCRIPT: &str =
    "for(const c of document.querySelectorAll('.stack td')){let i=0,p=c;while(p=p.previousElementSibling)i+=p.colSpan;c.dataset.l=c.closest('table').rows[0].cells[i]?.textContent||''}";
