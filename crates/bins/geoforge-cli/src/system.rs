//! The solar system and the focus: `show`, `positions`, `census`.

use crate::R;
use geoforge_almanac::ephemeris::Calendar;
use geoforge_core::body::{NodeKind, Star, StellarStage};
use geoforge_core::grid::GridSpec;
use geoforge_core::magnetism::MagneticField;
use geoforge_core::BodyContext;
use geoforge_store::World;
use std::fmt::Write;
use std::path::Path;

pub fn show(w: &World) -> R<()> {
    print!("{}", show_text(w)?);
    Ok(())
}

/// The system and the focus, as `show` prints them.
pub fn show_text(w: &World) -> R<String> {
    let mut out = String::new();
    let body = w.load_body()?;
    if let Some(sys) = w.load_system()? {
        writeln!(out, "system {} (seed {}), {:.2} Gyr old", sys.root.name, w.manifest.root_seed.0, sys.root.epoch.system_age_gyr)?;
        if let Some(jd) = sys.root.present_jd_tt {
            writeln!(out, "  the present: {} (its real calendar)", geoforge_almanac::ephemeris::utc_date(jd))?;
        }
        if let Some(site) = &sys.root.site {
            writeln!(out, "  drawn from the {}, {:.1} kpc from its centre", site.galaxy, site.radius_kpc)?;
        }
        for n in &sys.nodes {
            let a = n.elements.map_or(String::new(), |e| format!("  a={:.4} AU e={:.3}", e.a_au, e.e));
            let kind = match &n.node {
                NodeKind::Star(s) => describe_star(s),
                NodeKind::Body(b) => {
                    // A real body without an archetype is not described from a drawn interior.
                    let context = if sys.root.present_jd_tt.is_some() { geoforge_cosmos::sol::known_context(&sys, n.key) } else { geoforge_cosmos::body_context(&sys, n.key) };
                    let what = context.map_or_else(|| "measured orbit and size".into(), |c| geoforge_cosmos::describe(&c));
                    format!("{:?} {:.3} M⊕ {:.0} km — {what}", b.class, b.mass_kg / geoforge_core::body::EARTH_MASS_KG, b.radius_km)
                }
                NodeKind::Barycenter => "barycentre of the pair".into(),
                NodeKind::Belt(p) => format!("belt {:.2}–{:.2} AU", p.inner_au, p.outer_au),
                NodeKind::Ring(p) => format!("rings {:.0}–{:.0} km", p.inner_au * geoforge_core::body::AU_KM, p.outer_au * geoforge_core::body::AU_KM),
            };
            let mark = if n.key == w.manifest.focus { " ← focus" } else { "" };
            writeln!(out, "  {:<14} {kind} {a}{mark}", n.name)?;
            for note in geoforge_cosmos::remnants::notes(&sys, n.key) {
                writeln!(out, "  {:<14}   {note}", "")?;
            }
        }
    }
    writeln!(
        out,
        "
focus {}: radius {:.0} km, g {:.2} m/s², day {:.1} h, obliquity {:.1}°, flux {:.0} W/m² ({:.2}× Earth), ocean {:.2} km GEL, surface pressure {:.2} bar",
        body.name,
        body.exterior.radius_km,
        body.gravity_m_s2(),
        body.exterior.rotation_h,
        body.exterior.obliquity_deg,
        body.irradiance.flux_w_m2,
        body.irradiance.flux_w_m2 / 1361.0,
        body.interior.surface_ocean_gel_km(),
        body.interior.surface_pressure_bar()
    )?;
    if let Some(c) = geoforge_cosmos::worlds::climate_estimate(&body) {
        let h = match geoforge_cosmos::worlds::habitability(&body) {
            geoforge_cosmos::worlds::Habitability::Temperate => "temperate",
            geoforge_cosmos::worlds::Habitability::Habitable => "habitable",
            geoforge_cosmos::worlds::Habitability::None => "not habitable",
        };
        writeln!(out, "climate estimate: mean {:.1} °C, {:.0}% of the surface above freezing: {h}", c.mean_c, 100.0 * c.thawed_share)?;
    }
    if body.interior.has_atmosphere() {
        let rot = body.exterior.rotation_h;
        if let Some(s) = geoforge_surface::atmosphere::Substellar::for_body(&body) {
            writeln!(out, "circulation: tidally locked; air rises under the star (0°N {:.0}°E) and returns along the ground from the night side; equatorial jet {:.1}", s.lon, s.superrotation)?;
        } else if body.interior.has_solid_surface() {
            let cells = geoforge_surface::atmosphere::Cells::for_body(rot, body.exterior.radius_km, 0.0, 0.0);
            writeln!(out, "circulation: {} cells per hemisphere; eddy-driven jets: {}; Hadley cell to {:.0}°", cells.count(), cells.jets, cells.hadley)?;
        } else {
            let p = geoforge_surface::bands::Profile::for_body(&body, &geoforge_surface::bands::BandsParams::default());
            writeln!(out, "bands: jets per hemisphere: {}; equatorial wind {:+.0} m/s", p.jets, p.wind(0.0))?;
        }
    }
    writeln!(out, "magnetic field: {}", describe_field(&body))?;
    if body.interior.has_atmosphere() && body.interior.has_solid_surface() {
        writeln!(out, "air: {}", describe_air(&body))?;
    }
    if !body.tides.is_empty() {
        writeln!(out, "tides: {}", describe_tides(&body))?;
    }
    writeln!(out, "calendar: {}", describe_calendar(&calendar_of(&body)))?;
    if let GridSpec::Equirect { width, height } = w.manifest.grid {
        writeln!(out, "grid {width}×{height} equirectangular ({:.0} km cells at the equator)", 2.0 * std::f64::consts::PI * body.exterior.radius_km / f64::from(width))?;
    }
    Ok(out)
}

/// "30 µT dipole tilted 7°; magnetopause 9.5 radii against a 2.0 nPa wind; auroras at 71°
/// magnetic latitude", or why there is none.
fn describe_field(body: &BodyContext) -> String {
    let f = MagneticField::of(body);
    let npa = |p: f64| if p < 0.1 { format!("{p:.1e} nPa") } else { format!("{p:.1} nPa") };
    let onto = if body.interior.has_atmosphere() { "air" } else { "surface" };
    if !f.dynamo {
        return match f.lifetime_gyr {
            Some(life) => format!("none: the core dynamo ran {life:.1} Gyr and stopped {:.1} Gyr ago; the star's wind ({}) reaches the {onto}", body.epoch.system_age_gyr - life, npa(f.wind_npa)),
            None => format!("none (no convecting core); the star's wind ({}) reaches the {onto}", npa(f.wind_npa)),
        };
    }
    let shape = if f.dipolar { "dipole" } else { "multipolar" };
    let mut out = format!("{:.0} µT {shape} tilted {:.0}°", f.surface_ut, f.tilt_deg);
    if let Some(r) = f.magnetopause_radii {
        out += &format!("; magnetopause {r:.1} radii against a {} wind", npa(f.wind_npa));
    }
    out += &match f.aurora_lat_deg {
        Some(lat) => format!("; auroras at {lat:.0}° magnetic latitude"),
        None => format!("; the wind presses the field down to the {onto}: no auroral ovals"),
    };
    out
}

/// Whether people could breathe it: O₂ 0.16–0.5 bar (thin enough air starves, more is
/// toxic over time), CO₂ under 0.005 bar (5000 ppm at sea level), N₂ under 3 bar (narcosis).
fn describe_air(body: &BodyContext) -> String {
    use geoforge_core::body::Species;
    let pp = |s| body.interior.partial_pressure_bar(s);
    let (o2, co2, n2) = (pp(Species::Oxygen), pp(Species::CarbonDioxide), pp(Species::Nitrogen));
    let mut problems = Vec::new();
    if o2 < 0.16 {
        problems.push(format!("O₂ {o2:.2} bar, under 0.16"));
    } else if o2 > 0.5 {
        problems.push(format!("O₂ {o2:.2} bar, over 0.5"));
    }
    if co2 > 0.005 {
        problems.push(format!("CO₂ {co2:.4} bar, over 0.005"));
    }
    if n2 > 3.0 {
        problems.push(format!("N₂ {n2:.1} bar, over 3 (narcosis)"));
    }
    let g = body.gravity_m_s2() / 9.81;
    if problems.is_empty() {
        format!("breathable (O₂ {o2:.2} bar, CO₂ {co2:.4} bar, N₂ {n2:.2} bar); gravity {g:.2} g")
    } else {
        format!("not breathable: {}; gravity {g:.2} g", problems.join(", "))
    }
}

/// Equilibrium tide heights (∝ M/d³ · R⁴/M_planet) in units of Earth's lunar tide, and the
/// spring and neap tides of those that rise and fall (a locked world's star tide stands still).
fn describe_tides(body: &BodyContext) -> String {
    let height = |m: f64, d_km: f64| m / d_km.powi(3) * body.exterior.radius_km.powi(4) / body.exterior.mass_kg;
    let earth_moon = 7.342e22 / 384_400f64.powi(3) * 6371f64.powi(4) / 5.972e24;
    let mut each: Vec<(&str, f64)> = body.tides.iter().map(|t| (t.name.as_str(), height(t.mass_kg, t.distance_km) / earth_moon)).collect();
    // The star comes last in the list; a locked world keeps one face to it.
    let locked = body.is_tidally_locked();
    let moving: Vec<f64> = each.iter().take(each.len() - usize::from(locked)).map(|t| t.1).collect();
    each.sort_by(|a, b| b.1.total_cmp(&a.1));
    let mut out = each.iter().map(|(n, h)| format!("{n} {h:.2}")).collect::<Vec<_>>().join(", ") + " (Earth's Moon = 1)";
    if moving.len() >= 2 {
        let spring: f64 = moving.iter().sum();
        let neap = (2.0 * moving.iter().copied().fold(0.0, f64::max) - spring).abs();
        out += &format!("; spring {spring:.2}, neap {neap:.2} (Earth 1.46, 0.54)");
    }
    if locked {
        out += "; the star's is a fixed bulge (the world is locked to it)";
    }
    out
}

/// The focus's calendar; its moons are the tidal partners lighter than it.
pub fn calendar_of(body: &BodyContext) -> Calendar {
    Calendar::of_body(body)
}

/// "solar day 24.0 h; year 365.3 days, 365.2 local days; month (Moon) 29.5 days, 12.4 a year".
fn describe_calendar(c: &Calendar) -> String {
    let mut out = match (c.solar_day_h, c.days_per_year()) {
        (Some(day), Some(n)) if day < 72.0 => format!("solar day {day:.1} h; year {:.1} days, {n:.1} local days", c.year_days),
        (Some(day), Some(n)) => format!("solar day {:.1} days; year {:.1} days, {n:.2} local days", day / 24.0, c.year_days),
        _ => format!("no days (one face to the star); year {:.1} days", c.year_days),
    };
    for (name, m) in &c.months {
        out += &format!("; month ({name}) {m:.1} days, {:.1} a year", c.year_days / m);
    }
    out
}

/// Positions on `day`, printed and optionally written as JSON.
pub fn positions(w: &World, day: f64, to: Option<&Path>) -> R<()> {
    let sys = w.load_system()?.ok_or("this world has no solar system (a hand-made body)")?;
    let body = w.load_body()?;
    let calendar = calendar_of(&body);
    let (year, local) = calendar.date(day);
    let unit = if calendar.solar_day_h.is_some() { "local day" } else { "day" };
    println!("day {day} (Earth days from the present): {}, year {year}, {unit} {local:.2}", body.name);
    let placed = geoforge_almanac::ephemeris::positions(&sys, day);
    let at = |k: geoforge_core::Key| placed.iter().find(|p| p.key == k).map(|p| p.position_au);
    println!("  {:<16} {:>11} {:>11} {:>11} {:>13} {:>11}", "", "x AU", "y AU", "z AU", "from parent", "period d");
    for p in &placed {
        let (from, period) = match p.period_days {
            0.0 => ("—".to_string(), "—".to_string()),
            d if p.distance_au < 0.01 => (format!("{:.0} km", p.distance_au * geoforge_core::body::AU_KM), format!("{d:.2}")),
            d => (format!("{:.4} AU", p.distance_au), format!("{d:.2}")),
        };
        let [x, y, z] = p.position_au;
        println!("  {:<16} {x:>11.5} {y:>11.5} {z:>11.5} {from:>13} {period:>11}", p.name);
    }
    // The star (or the pair's barycentre) the focus's planet orbits sits at the root.
    let star = [0.0; 3];
    let focus = at(w.manifest.focus).ok_or("the focus is not in the system")?;
    let mut lit = Vec::new();
    for m in sys.children(w.manifest.focus) {
        if let Some(pos) = at(m.key) {
            let f = geoforge_almanac::ephemeris::lit_fraction(pos, focus, star);
            println!("  {} is {:.0}% lit", m.name, 100.0 * f);
            lit.push(serde_json::json!({ "name": m.name, "lit": f }));
        }
    }
    if let Some(path) = to {
        let bodies: Vec<serde_json::Value> = placed
            .iter()
            .map(|p| {
                serde_json::json!({
                    "key": p.key.0, "name": p.name, "parent": p.parent.map(|k| k.0),
                    "position_au": p.position_au, "distance_from_parent_au": p.distance_au, "period_days": p.period_days,
                })
            })
            .collect();
        let doc = serde_json::json!({ "day": day, "focus": body.name, "date": { "year": year, "day": local, "unit": unit }, "moons": lit, "bodies": bodies });
        std::fs::write(path, serde_json::to_string_pretty(&doc)?)?;
        println!("written to {}", path.display());
    }
    Ok(())
}

/// "white dwarf DA8, 0.59 M☉, 6300 K, 2.1e-4 L☉, [Fe/H] −0.12 (born 1.80 M☉)".
fn describe_star(s: &Star) -> String {
    let lum = if s.luminosity_solar >= 0.01 { format!("{:.3} L☉", s.luminosity_solar) } else { format!("{:.1e} L☉", s.luminosity_solar) };
    let born = match s.initial_mass_solar {
        Some(m0) if (m0 - s.mass_solar).abs() > 0.01 => format!(" (born {m0:.2} M☉)"),
        _ => String::new(),
    };
    let body = if s.stage == StellarStage::BlackHole { String::new() } else { format!(", {:.0} K, {lum}", s.temperature_k) };
    format!("{} {}, {:.2} M☉{body}, [Fe/H] {:+.2}{born}", s.stage.label(), s.spectral_class, s.mass_solar, s.metallicity_dex)
}

pub fn print_census(c: &geoforge_cosmos::census::Census) {
    let n = c.systems().max(1) as f64;
    let pct = |x: u64, of: u64| if of == 0 { "-".to_string() } else { format!("{:.0}%", 100.0 * x as f64 / of as f64) };
    let per = |x: u64, of: u64| if of == 0 { "-".to_string() } else { format!("{:.2}", x as f64 / of as f64) };
    let seeds = |v: &[u64]| v.iter().map(u64::to_string).collect::<Vec<_>>().join(" ");
    println!("{} systems drawn from the {} (seeds {}..{})\n", c.systems(), c.galaxy, c.seeds.start, c.seeds.end);
    println!("{:<13} {:>7} {:>6} {:>6} {:>6} {:>6} {:>7} {:>5} {:>7}  {:<20} wet examples", "primary", "share", "multi", "rocky", "giant", "ice", "barren", "wet", "locked", "examples");
    for k in c.kinds.iter().filter(|k| k.systems > 0) {
        let solid = k.systems - k.barren;
        println!(
            "{:<13} {:>7} {:>6} {:>6} {:>6} {:>6} {:>7} {:>5} {:>7}  {:<20} {}",
            k.name,
            format!("{:.2}%", 100.0 * k.systems as f64 / n),
            pct(k.multiples, k.systems),
            per(k.rocky, solid),
            per(k.giants, solid),
            per(k.ice_giants, solid),
            pct(k.barren, k.systems),
            pct(k.wet, k.systems),
            pct(k.locked, solid),
            seeds(&k.examples),
            seeds(&k.wet_examples)
        );
    }
    println!(
        "\nbarren {} (no solid body: `new` fails for these seeds), wet focus {}, planets per system with any: rocky {}, giant {}, ice {}",
        pct(c.total(|k| k.barren), c.systems()),
        pct(c.total(|k| k.wet), c.systems()),
        per(c.total(|k| k.rocky), c.systems() - c.total(|k| k.barren)),
        per(c.total(|k| k.giants), c.systems() - c.total(|k| k.barren)),
        per(c.total(|k| k.ice_giants), c.systems() - c.total(|k| k.barren)),
    );
    println!("columns: multi = more than one star; rocky/giant/ice = planets per system with a solid body; wet = focus has surface water; locked = focus tidally locked to its star");
}
