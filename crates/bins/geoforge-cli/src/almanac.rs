//! A place and a day: `place`, `weather`, `almanac`, `weather-map`, `sky`.

use crate::system::calendar_of;
use crate::R;
use geoforge_pipeline as pipe;
use geoforge_store::World;
use std::path::{Path, PathBuf};

pub fn place(dir: &Path, lat: f64, lon: f64, radius_km: f64, from_day: f64, days: Option<f64>, json: bool) -> R<()> {
    let w = World::open(dir)?;
    let days = days.unwrap_or_else(|| w.load_body().map_or(365.25, |b| b.irradiance.year_days));
    let r = pipe::place::report(&w, lat, lon, radius_km, from_day, days)?;
    if json {
        println!("{}", serde_json::to_string_pretty(&r)?);
    } else {
        print_place(&r);
    }
    Ok(())
}

pub fn weather(dir: &Path, lat: f64, lon: f64, day: f64, days: usize, json: bool) -> R<()> {
    let f = geoforge_almanac::weather::forecast(&World::open(dir)?, lat, lon, day, days)?;
    if json {
        println!("{}", serde_json::to_string_pretty(&f)?);
    } else {
        print_weather(&f);
    }
    Ok(())
}

pub fn almanac(dir: &Path, lat: f64, lon: f64, year: i64, json: bool) -> R<()> {
    let a = geoforge_almanac::year::year(&World::open(dir)?, lat, lon, year)?;
    if json {
        println!("{}", serde_json::to_string_pretty(&a)?);
        return Ok(());
    }
    let hm = |h: f64| format!("{:02}:{:02}", h.floor() as i64, ((h.fract()) * 60.0).floor() as i64);
    println!("Almanac of {:.2}°, {:.2}°, year {}", a.lat, a.lon, a.year);
    let moons: Vec<&str> = a.days.first().map_or_else(Vec::new, |d| d.moons.iter().map(|m| m.0.as_str()).collect());
    println!("  {:>6} {:>5} {:>5} {:>6} {:>6} {:>8}  {}sky", "day", "rise", "set", "low", "high", "precip", moons.iter().map(|m| format!("{:>7} ", trim(m, 7))).collect::<String>());
    for d in &a.days {
        let w = &d.weather;
        let (rise, set) = match (d.sun, a.day_h) {
            (Some((r, s)), Some(h)) if s - r >= h - 1e-6 => ("sun".into(), "up".into()),
            (Some((r, s)), _) => (hm(r), hm(s)),
            (None, Some(_)) => ("—".into(), "—".into()),
            (None, None) => (String::new(), String::new()),
        };
        let precip = if w.precipitation_mm > 0.0 { format!("{:.1}", w.precipitation_mm) } else { "—".into() };
        let lit: String = d.moons.iter().map(|m| format!("{:>6.0}% ", 100.0 * m.1)).collect();
        let mut sky = w.storm.as_ref().map_or_else(|| w.sky.to_string(), |name| format!("{} ({name})", w.sky));
        for e in w.events.iter().filter(|e| e.kind != "thunderstorm" || w.sky != "thunderstorm") {
            sky += &format!("; {}", e.what);
        }
        println!("  {:>6.0} {rise:>5} {set:>5} {:>6.1} {:>6.1} {precip:>8}  {lit}{sky}", (w.local_day + 1e-6).floor(), w.low_c, w.high_c);
    }
    println!("the year's months: rain, drought index (SPI-3, against the {} years before), soil water", a.droughts.years);
    for m in &a.months {
        let spi = m.spi3.map_or("  too dry".to_string(), |s| format!("{s:>+9.2}"));
        println!(
            "  day {:>7.0} {:>7.0} mm {spi} {:>5.0}%{}",
            m.day,
            m.rain_mm,
            100.0 * m.soil,
            m.spi3.filter(|&s| s <= geoforge_almanac::drought::MODERATE_SPI).map_or("", |s| if s <= geoforge_almanac::drought::SEVERE_SPI { "  severe drought" } else { "  drought" })
        );
    }
    print_droughts(&a.droughts.years, a.droughts.moderate_share, a.droughts.severe_share, a.droughts.longest_months, a.droughts.dry_soil_share);
    Ok(())
}

fn trim(s: &str, n: usize) -> String {
    s.chars().take(n).collect()
}

fn print_droughts(years: &usize, moderate: f64, severe: f64, longest: usize, dry_soil: Option<f64>) {
    println!(
        "droughts over {years} years: {:.0}% of months moderate (SPI-3 ≤ −1), {:.0}% severe (≤ −1.5), longest {longest} months{}",
        100.0 * moderate,
        100.0 * severe,
        dry_soil.map_or(String::new(), |s| format!("; soil under a quarter full in {:.0}% of growing-season months", 100.0 * s))
    );
}

pub fn weather_map(dir: &Path, day: f64, to: Option<PathBuf>, width: u32) -> R<()> {
    let w = World::open(dir)?;
    let to = to.unwrap_or_else(|| w.root.join("renders").join(format!("weather_day{day}.png")));
    if let Some(parent) = to.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let storms = geoforge_almanac::weather::weather_map(&w, day, width, &to)?;
    println!("{} storms on day {day}:", storms.len());
    for s in &storms {
        let kind = if s.kind == geoforge_almanac::storms::StormKind::Low { "low" } else { "tropical cyclone" };
        println!("  {:<12} {kind:<16} {:>6.1}° {:>7.1}°  {:>4.0} hPa deep, {:.1} days old", s.name, s.lat, s.lon, s.depth_hpa, s.age_days);
    }
    println!("written {}", to.display());
    Ok(())
}

fn print_place(r: &pipe::place::PlaceReport) {
    println!("{:.2}°, {:.2}° within {:.0} km ({:.0} km²): {} {:.0} m", r.lat, r.lon, r.radius_km, r.area_km2, if r.at_sea { "sea," } else { "land," }, r.elevation_m);
    if let Some(k) = r.koppen {
        let soil = r.soil.unwrap_or("—");
        println!("climate {k}; soil {soil} (fertility {:.1}); timber {:.1}", r.fertility.unwrap_or(0.0), r.timber.unwrap_or(0.0));
    }
    if !r.deposits.is_empty() {
        println!("deposits: expected here, chance of any, how many a local map holds; settings at the centre");
        for d in &r.deposits {
            let settings = d.settings.iter().map(|(n, l)| format!("{n} {l:.2}")).collect::<Vec<_>>().join(", ");
            println!("  {:<18} {:>6.2} {:>5.0}% {:>3}  {settings}", d.kind, d.expected, 100.0 * d.chance, d.drawn);
        }
    }
    for (d, km) in &r.giants {
        println!("  giant {} deposit ({}), {:.0} Mt of ore, {:.0} m down, host rock {:.0} Ma, {km:.0} km away", d.kind, d.setting, d.ore_mt, d.depth_m, d.age_ma);
    }
    println!("hazards over {:.0} days from day {:.0}: a year, in the window, chance of any", r.days, r.from_day);
    for h in &r.hazards {
        println!("  {:<12} {:>9.3} {:>9.3} {:>5.1}%  ({})", h.hazard, h.per_year, h.expected, 100.0 * h.chance, h.event);
    }
    println!("  magnitude 7+: {:.2}%, 8+: {:.3}%; largest likely {:.1}", 100.0 * r.quake_7, 100.0 * r.quake_8, r.largest_quake);
    if let Some(d) = &r.droughts {
        print_droughts(&d.years, d.moderate_share, d.severe_share, d.longest_months, d.dry_soil_share);
        if let Some(spi) = d.spi_now {
            println!("  drought index (SPI-3) as the window starts: {spi:+.2}");
        }
    }
}

fn print_weather(f: &geoforge_almanac::weather::Forecast) {
    let step = if (f.step_days - 1.0).abs() < 1e-9 { "Earth day".to_string() } else { format!("local day of {:.1} h", 24.0 * f.step_days) };
    println!("Weather at {:.2}°, {:.2}°{}, one row per {step}", f.lat, f.lon, if f.at_sea { " (at sea)" } else { "" });
    println!("  {:>9} {:>12} {:>6} {:>6} {:>6} {:>8} {:>10}  sky", "day", "date", "low", "mean", "high", "precip", "wind");
    for d in &f.days {
        let precip = if d.precipitation_mm > 0.0 { format!("{:.1} mm", d.precipitation_mm) } else { "—".into() };
        println!(
            "  {:>9.2} {:>12} {:>6.1} {:>6.1} {:>6.1} {precip:>8} {:>4.0} m/s {:>3}  {}",
            d.day,
            format!("y{} d{:.0}", d.year, (d.local_day + 1e-6).floor()),
            d.low_c,
            d.mean_c,
            d.high_c,
            d.wind_m_s,
            compass(d.wind_from_deg),
            {
                let mut sky = d.storm.as_ref().map_or_else(|| d.sky.to_string(), |name| format!("{} ({name}, {:.0} hPa)", d.sky, d.pressure_hpa));
                for e in d.events.iter().filter(|e| e.kind != "thunderstorm" || d.sky != "thunderstorm") {
                    sky += &format!("; {}", e.what);
                }
                sky
            }
        );
    }
}

/// "NNE" for 22.5°.
fn compass(deg: f64) -> &'static str {
    const POINTS: [&str; 16] = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE", "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"];
    POINTS[((deg.rem_euclid(360.0) / 22.5).round() as usize) % 16]
}

/// The sky on `day`, printed, and with `to` charted and catalogued.
pub fn sky(w: &World, day: f64, to: Option<&Path>, width: u32) -> R<()> {
    use geoforge_almanac::sky::Kind;
    let sys = w.load_system()?.ok_or("this world has no solar system (a hand-made body)")?;
    let body = w.load_body()?;
    let calendar = calendar_of(&body);
    let s = geoforge_almanac::sky::sky(&sys, w.manifest.focus, day, calendar.year_days)?;
    let (year, local) = calendar.date(day);
    let date = |t: f64| geoforge_almanac::ephemeris::real_date(&sys, t).map_or(String::new(), |d| format!(", {d}"));
    println!("The sky of {} on day {day} (year {year}, day {local:.1}{}): {} naked-eye stars; {} {:.1} kpc from its centre", body.name, date(day), s.stars.len(), s.galaxy, s.radius_kpc);
    let radec = |o: &geoforge_almanac::sky::Object| format!("RA {:>5.1}° Dec {:>+5.1}°", o.ra_deg, o.dec_deg);
    println!("Brightest stars:");
    for o in s.stars.iter().take(12) {
        let name = o.name.clone().unwrap_or_default();
        println!("  {name:<12} {:>5.2}  {:<7} {:>6.0} ly  {}", o.v, o.spectral_class.as_deref().unwrap_or(""), o.distance, radec(o));
    }
    println!("The system:");
    for o in &s.bodies {
        let name = o.name.clone().unwrap_or_default();
        let kind = match o.kind {
            Kind::Sun => "star",
            Kind::Planet => "planet",
            Kind::Moon => "moon",
            Kind::Star => "",
        };
        let seen = match (o.kind, o.elongation_deg) {
            (Kind::Sun, None) => "the sun".to_string(),
            _ if o.v > geoforge_cosmos::sky::NAKED_EYE => "too faint".to_string(),
            _ if !o.visible() => "lost in the glare".to_string(),
            _ => "in view at night".to_string(),
        };
        let elong = o.elongation_deg.map_or(String::new(), |e| format!(", {e:.0}° from the sun"));
        let lit = o.lit.map_or(String::new(), |l| format!(", {:.0}% lit", 100.0 * l));
        println!("  {name:<12} {kind:<6} {:>7.2}  {}{elong}{lit}: {seen}", o.v, radec(o));
    }
    if s.eclipses.is_empty() {
        println!("No eclipses in the year from day {day}.");
    } else {
        println!("Eclipses in the year from day {day}:");
        for e in &s.eclipses {
            let (y, d) = calendar.date(e.day);
            println!(
                "  day {:>8.2} (year {y}, day {d:.1}{}): {:?} {:?} eclipse of {} by {}, magnitude {:.2}",
                e.day,
                date(e.day),
                e.coverage,
                e.kind,
                if e.kind == geoforge_almanac::ephemeris::EclipseKind::Solar { &e.star } else { &e.body },
                if e.kind == geoforge_almanac::ephemeris::EclipseKind::Solar { &e.body } else { &body.name },
                e.magnitude
            );
        }
    }
    if let Some(dir) = to {
        std::fs::create_dir_all(dir)?;
        let (gw, gh) = (360, 180);
        let glow = s.glow(gw, gh);
        let points: Vec<geoforge_store::render::SkyPoint> = s
            .stars
            .iter()
            .chain(s.bodies.iter().filter(|b| b.v <= geoforge_cosmos::sky::NAKED_EYE))
            .map(|o| geoforge_store::render::SkyPoint { ra_deg: o.ra_deg, dec_deg: o.dec_deg, v: o.v, temperature_k: o.temperature_k, ring: matches!(o.kind, Kind::Planet | Kind::Moon) })
            .collect();
        let img = geoforge_store::render::render_sky(width, &glow, gw, &points, &s.ecliptic);
        let chart = dir.join("sky.png");
        img.save(&chart)?;
        let json = dir.join("sky.json");
        std::fs::write(&json, serde_json::to_string_pretty(&s)?)?;
        println!("written {} and {}", chart.display(), json.display());
    }
    Ok(())
}
