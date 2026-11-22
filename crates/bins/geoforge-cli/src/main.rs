//! `geoforge-v2`: make a world, generate it, look at it, edit any layer and regenerate.

use clap::{Parser, Subcommand};
use geoforge_core::grid::GridSpec;
use geoforge_core::product::{LayerId, ProductId};
use geoforge_core::body::{NodeKind, Star, StellarStage};
use geoforge_core::ephemeris::Calendar;
use geoforge_core::magnetism::MagneticField;
use geoforge_core::{BodyContext, Seed};
use geoforge_cosmos::{GalaxyPreset, Recipe};
use geoforge_pipeline::{self as pipe, GenerateOptions, LayerAction, NewWorld, State, View};
use geoforge_store::World;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

#[derive(Parser)]
#[command(name = "geoforge-v2", version, about = "Procedural worlds you can edit at any layer")]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Create a world: a solar system drawn from a galaxy and its most Earth-like body (or
    /// the habitable recipe, a real body, or your own body file).
    New {
        dir: PathBuf,
        #[arg(long, default_value_t = 42)]
        seed: u64,
        /// The galaxy to draw the system from: milky-way, elliptical, dwarf, starburst.
        #[arg(long, default_value = "milky-way", value_parser = parse_galaxy)]
        galaxy: GalaxyPreset,
        /// Force a habitable system: a Sun-like star (F8–K5) with an Earth-like world in its
        /// habitable zone, whatever the galaxy would draw.
        #[arg(long, conflicts_with_all = ["body", "archetype"])]
        habitable: bool,
        /// Grid width in cells (height is half).
        #[arg(long, default_value_t = 2048)]
        width: u32,
        /// A body file (TOML or JSON BodyContext) instead of a generated system.
        #[arg(long)]
        body: Option<PathBuf>,
        /// Start from a real body instead of a generated system: earth, mars, venus, mercury,
        /// moon, titan, europa, gas_giant, ice_giant.
        #[arg(long, conflicts_with = "body")]
        archetype: Option<geoforge_core::Archetype>,
        /// Focus on another body of the generated system: its name (see `show`), or
        /// "<belt name> #<n>" for a belt's n-th largest member.
        #[arg(long, conflicts_with_all = ["body", "archetype"])]
        focus: Option<String>,
    },
    /// Generate (or update) every product that is missing or stale.
    Generate {
        dir: PathBuf,
        /// Rebuild everything.
        #[arg(long)]
        force: bool,
        /// Stop after this layer.
        #[arg(long)]
        until: Option<String>,
        /// Rebuild this layer and everything after it.
        #[arg(long)]
        from: Option<String>,
        /// Also write preview renders to <dir>/renders.
        #[arg(long)]
        render: bool,
    },
    /// Show what generate would do and why.
    Plan { dir: PathBuf },
    /// Sources and states of every product.
    Status { dir: PathBuf },
    /// Write edit views (PNG + world file + sidecar) of products for editing elsewhere, and
    /// derived views (koppen) as PNG + legend + world file.
    Export {
        dir: PathBuf,
        /// Product names, layer names, view names, or "all".
        what: Vec<String>,
        #[arg(long)]
        to: Option<PathBuf>,
        /// Write rasters this many times the grid's size (each cell repeated), for detail.
        #[arg(long, default_value_t = 1)]
        scale: u32,
        /// Also write polar cap views (<product>.north.png, .south.png) for painting the
        /// poles without the map's stretching.
        #[arg(long)]
        caps: bool,
    },
    /// Bring an edited view back: as a patch on the generated product (default) or as authored data.
    Import {
        dir: PathBuf,
        product: String,
        #[arg(long)]
        from: Option<PathBuf>,
        /// Replace the product with your data instead of patching it.
        #[arg(long)]
        authored: bool,
    },
    /// Discard a product's edits; it is generated again.
    Reset { dir: PathBuf, product: String },
    /// Write preview renders of stored products and derived views (koppen).
    Render {
        dir: PathBuf,
        /// Product names, layer names, view names, or "all" (the default).
        what: Vec<String>,
        #[arg(long)]
        to: Option<PathBuf>,
    },
    /// Describe the system and the focus body.
    Show { dir: PathBuf },
    /// What is at a place: climate, soil, timber, the deposits likely within reach (and how
    /// many a local map holds), and the chance of each hazard over a stretch of days.
    Place {
        dir: PathBuf,
        #[arg(long, allow_hyphen_values = true)]
        lat: f64,
        #[arg(long, allow_hyphen_values = true)]
        lon: f64,
        #[arg(long, default_value_t = 100.0)]
        radius_km: f64,
        /// The stretch of days: from this day (Earth days from the present), for this many
        /// (default: one of the world's years).
        #[arg(long, default_value_t = 0.0)]
        from_day: f64,
        #[arg(long)]
        days: Option<f64>,
        /// Print JSON instead.
        #[arg(long)]
        json: bool,
    },
    /// Where every star, planet and moon is on a day (Earth days from the present), the
    /// focus's date on its own calendar, and how lit its moons are.
    Positions {
        dir: PathBuf,
        #[arg(long, default_value_t = 0.0)]
        day: f64,
        /// Also write them as JSON.
        #[arg(long)]
        to: Option<PathBuf>,
    },
    /// The focus's night sky on a day: its brightest stars, the planets and moons in view,
    /// and the eclipses of the year to come; with `--to`, a star chart (sky.png) and a
    /// catalogue of every naked-eye star (sky.json).
    Sky {
        dir: PathBuf,
        #[arg(long, default_value_t = 0.0)]
        day: f64,
        #[arg(long)]
        to: Option<PathBuf>,
        /// Chart width in pixels (half as tall).
        #[arg(long, default_value_t = 2048)]
        width: u32,
    },
    /// Day-by-day weather at a place from a day (Earth days from the present): drawn from
    /// its climate and the season of the date, the same every time for the same date.
    Weather {
        dir: PathBuf,
        #[arg(long, allow_hyphen_values = true)]
        lat: f64,
        #[arg(long, allow_hyphen_values = true)]
        lon: f64,
        #[arg(long, default_value_t = 0.0)]
        day: f64,
        /// How many days (local days on a world with days).
        #[arg(long, default_value_t = 14)]
        days: usize,
        #[arg(long)]
        json: bool,
    },
    /// Print the seed a product (or a product's entity) is generated from.
    Seed {
        dir: PathBuf,
        product: String,
        /// Owner key (defaults to the focus body).
        #[arg(long)]
        owner: Option<u64>,
    },
    /// Summary statistics of stored products (for checking realism).
    Stats { dir: PathBuf },
    /// Draw many systems from a galaxy and count them: star kinds, planets, wet worlds, and
    /// example seeds of each kind for `new`.
    Census {
        #[arg(long, default_value = "milky-way", value_parser = parse_galaxy)]
        galaxy: GalaxyPreset,
        /// How many seeds, and the first.
        #[arg(long, default_value_t = 10_000)]
        count: u64,
        #[arg(long, default_value_t = 0)]
        from: u64,
    },
}

fn parse_galaxy(name: &str) -> Result<GalaxyPreset, String> {
    GalaxyPreset::from_name(name).ok_or_else(|| format!("unknown galaxy `{name}` (one of: {})", GalaxyPreset::ALL.map(GalaxyPreset::name).join(", ")))
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

type R<T> = Result<T, Box<dyn std::error::Error>>;

/// Products and derived views by name (or a layer's products); everything for "all" or none.
fn parse_targets(what: &[String]) -> R<(Vec<ProductId>, Vec<View>)> {
    if what.is_empty() || what.iter().any(|w| w == "all") {
        return Ok((ProductId::ALL.to_vec(), View::all()));
    }
    let (mut ps, mut vs) = (Vec::new(), Vec::new());
    for w in what {
        if let Some(p) = ProductId::from_name(w) {
            ps.push(p);
        } else if let Some(l) = LayerId::from_name(w) {
            ps.extend(l.products());
        } else if let Some(v) = View::from_name(w) {
            vs.push(v);
        } else if let Some(group) = View::group(w) {
            vs.extend(group);
        } else {
            let views = View::all().into_iter().map(View::name).collect::<Vec<_>>().join(", ");
            return Err(format!("unknown product, layer or view `{w}` (views: {views}, or the groups climate_maps, oceans, resources and hazards; products: {})", ProductId::ALL.map(ProductId::name).join(", ")).into());
        }
    }
    ps.sort();
    ps.dedup();
    vs.sort();
    vs.dedup();
    Ok((ps, vs))
}

fn product(name: &str) -> R<ProductId> {
    ProductId::from_name(name).ok_or_else(|| format!("unknown product `{name}`").into())
}

fn read_body(path: &Path) -> R<BodyContext> {
    let text = std::fs::read_to_string(path)?;
    Ok(if path.extension().is_some_and(|e| e == "json") { serde_json::from_str(&text)? } else { toml::from_str(&text)? })
}

fn run(cli: Cli) -> R<()> {
    match cli.cmd {
        Cmd::New { dir, seed, galaxy, habitable, width, body, archetype, focus } => {
            if seed > i64::MAX as u64 {
                return Err("seed must be at most 9223372036854775807".into());
            }
            let mut body = body.as_deref().map(read_body).transpose()?;
            if let Some(a) = archetype {
                body = Some(a.context(geoforge_core::Key(1), a.name()));
            }
            let recipe = if habitable { Recipe::Habitable } else { Recipe::Galaxy(galaxy.galaxy()) };
            let w = pipe::create_world(&dir, &NewWorld { seed: Seed(seed), recipe, grid: GridSpec::Equirect { width, height: width / 2 }, body, focus })?;
            show(&w)?;
            println!("\ncreated {} — next: geoforge-v2 generate {}", dir.display(), dir.display());
        }
        Cmd::Generate { dir, force, until, from, render } => {
            let mut w = World::open(&dir)?;
            let until = until.map(|u| LayerId::from_name(&u).ok_or(format!("unknown layer `{u}`"))).transpose()?;
            let from = from.map(|u| LayerId::from_name(&u).ok_or(format!("unknown layer `{u}`"))).transpose()?;
            let (products, _) = pipe::generate(&mut w, &GenerateOptions { force, until, from }, |r| {
                let what = match &r.action {
                    LayerAction::Ran { seconds } => format!("generated in {seconds:.1}s"),
                    LayerAction::Reused => "up to date".into(),
                    LayerAction::NotApplicable { reason } => format!("does not apply: {reason}"),
                };
                println!("{:<10} {what}", r.layer.name());
                for n in &r.notes {
                    println!("           note: {n}");
                }
            })?;
            if render {
                let files = pipe::render_all(&w, &products, &w.root.join("renders"), &ProductId::ALL, &View::all())?;
                println!("renders: {} files in {}", files.len(), w.root.join("renders").display());
            }
        }
        Cmd::Plan { dir } => {
            let w = World::open(&dir)?;
            let plan = pipe::plan(&w)?;
            for lp in &plan.layers {
                let by = lp.implementation.map_or(String::new(), |i| format!(" ({})", i.name()));
                let action = match (&lp.not_applicable, lp.needs_run()) {
                    (Some(r), _) => format!("skip — {r}"),
                    (None, true) => format!("run{by}"),
                    (None, false) => format!("reuse{by}"),
                };
                println!("{:<11} {action}", lp.layer.name());
                for p in &lp.products {
                    println!("  {:<18} {:<9} {:?}", p.product.name(), p.source, p.state);
                }
            }
        }
        Cmd::Status { dir } => {
            let w = World::open(&dir)?;
            let plan = pipe::plan(&w)?;
            let edits = pipe::edit_status(&w);
            println!("{:<18} {:<10} {:<10} edit", "product", "source", "state");
            for lp in &plan.layers {
                for p in &lp.products {
                    let e = edits.get(p.product.name()).map_or(String::new(), |e| {
                        let mut s = Vec::new();
                        if e.rebased {
                            s.push("rebased".to_string());
                        }
                        if e.orphaned > 0 {
                            s.push(format!("{} orphaned", e.orphaned));
                        }
                        s.join(", ")
                    });
                    let state = match p.state {
                        State::Fresh => "fresh",
                        State::Stale => "stale",
                        State::Missing => "missing",
                        State::NotApplicable => "n/a",
                    };
                    println!("{:<18} {:<10} {:<10} {e}", p.product.name(), p.source, state);
                }
            }
        }
        Cmd::Export { dir, what, to, scale, caps } => {
            let w = World::open(&dir)?;
            let to = to.unwrap_or_else(|| pipe::edit_folder(&w));
            let (ps, vs) = parse_targets(&what)?;
            let ps: Vec<ProductId> = ps.into_iter().filter(|p| w.has_product(*p)).collect();
            for (p, files) in pipe::export(&w, &ps, &to, &geoforge_store::edit::ViewOptions { scale, caps })? {
                println!("{:<18} {}", p.name(), files.join(" "));
            }
            for (v, files) in pipe::export_views(&w, &vs, &to)? {
                println!("{:<18} {} (derived; not imported)", v.name(), files.join(" "));
            }
            if !ps.is_empty() {
                println!("edit views in {} — edit, then: geoforge-v2 import {} <product>", to.display(), dir.display());
            }
        }
        Cmd::Import { dir, product: name, from, authored } => {
            let mut w = World::open(&dir)?;
            let from = from.unwrap_or_else(|| pipe::edit_folder(&w));
            let r = pipe::import(&mut w, product(&name)?, &from, authored)?;
            println!("{}: imported as {} ({} cells changed, {} rows changed)", r.product.name(), r.mode, r.changed_cells, r.rows_changed);
            if let Some((iw, ih)) = r.resampled_from {
                println!("  the image was {iw}x{ih}; resampled to the grid");
            }
            if r.cap_cells > 0 {
                println!("  {} cells came from the polar cap views", r.cap_cells);
            }
            if r.snapped > 0 {
                println!("  {} pixels were not exact palette colours and were snapped to the nearest one", r.snapped);
            }
            println!("next: geoforge-v2 generate {} (downstream layers will update)", dir.display());
        }
        Cmd::Reset { dir, product: name } => {
            let mut w = World::open(&dir)?;
            pipe::reset(&mut w, product(&name)?)?;
            println!("{name}: edits removed; it will be generated");
        }
        Cmd::Render { dir, what, to } => {
            let w = World::open(&dir)?;
            let products = pipe::load_all(&w)?;
            let (ps, vs) = parse_targets(&what)?;
            let to = to.unwrap_or_else(|| w.root.join("renders"));
            for f in pipe::render_all(&w, &products, &to, &ps, &vs)? {
                println!("{}", f.display());
            }
        }
        Cmd::Show { dir } => show(&World::open(&dir)?)?,
        Cmd::Positions { dir, day, to } => positions(&World::open(&dir)?, day, to.as_deref())?,
        Cmd::Sky { dir, day, to, width } => sky(&World::open(&dir)?, day, to.as_deref(), width)?,
        Cmd::Weather { dir, lat, lon, day, days, json } => {
            let f = pipe::weather::forecast(&World::open(&dir)?, lat, lon, day, days)?;
            if json {
                println!("{}", serde_json::to_string_pretty(&f)?);
            } else {
                print_weather(&f);
            }
        }
        Cmd::Place { dir, lat, lon, radius_km, from_day, days, json } => {
            let w = World::open(&dir)?;
            let days = days.unwrap_or_else(|| w.load_body().map_or(365.25, |b| b.irradiance.year_days));
            let r = pipe::place::report(&w, lat, lon, radius_km, from_day, days)?;
            if json {
                println!("{}", serde_json::to_string_pretty(&r)?);
            } else {
                print_place(&r);
            }
        }
        Cmd::Seed { dir, product: name, owner } => {
            let w = World::open(&dir)?;
            let p = product(&name)?;
            let key = owner.map_or(w.manifest.focus, geoforge_core::Key);
            println!("{}", w.manifest.root_seed.for_product(p, key).0);
        }
        Cmd::Census { galaxy, count, from } => print_census(&geoforge_cosmos::census::census(&galaxy.galaxy(), from..from.saturating_add(count))),
        Cmd::Stats { dir } => {
            let w = World::open(&dir)?;
            let products = pipe::load_all(&w)?;
            let body = w.load_body()?;
            let grid = w.manifest.grid.build()?;
            let s = geoforge_pipeline::stats::summarize(&products, &grid, body.exterior.radius_km, body.interior.has_surface_liquid_water());
            println!("{}", serde_json::to_string_pretty(&s)?);
        }
    }
    Ok(())
}

fn show(w: &World) -> R<()> {
    let body = w.load_body()?;
    if let Some(sys) = w.load_system()? {
        println!("system {} (seed {}), {:.2} Gyr old", sys.root.name, w.manifest.root_seed.0, sys.root.epoch.system_age_gyr);
        if let Some(site) = &sys.root.site {
            println!("  drawn from the {}, {:.1} kpc from its centre", site.galaxy, site.radius_kpc);
        }
        for n in &sys.nodes {
            let a = n.elements.map_or(String::new(), |e| format!("  a={:.4} AU e={:.3}", e.a_au, e.e));
            let kind = match &n.node {
                NodeKind::Star(s) => describe_star(s),
                NodeKind::Body(b) => {
                    let what = geoforge_cosmos::body_context(&sys, n.key).map_or(String::new(), |c| geoforge_cosmos::describe(&c));
                    format!("{:?} {:.3} M⊕ {:.0} km — {what}", b.class, b.mass_kg / geoforge_core::body::EARTH_MASS_KG, b.radius_km)
                }
                NodeKind::Barycenter => "barycentre of the pair".into(),
                NodeKind::Belt(p) | NodeKind::Ring(p) => format!("belt {:.2}–{:.2} AU", p.inner_au, p.outer_au),
            };
            let mark = if n.key == w.manifest.focus { " ← focus" } else { "" };
            println!("  {:<14} {kind} {a}{mark}", n.name);
        }
    }
    println!(
        "\nfocus {}: radius {:.0} km, g {:.2} m/s², day {:.1} h, obliquity {:.1}°, flux {:.0} W/m² ({:.2}× Earth), ocean {:.2} km GEL, surface pressure {:.2} bar",
        body.name,
        body.exterior.radius_km,
        body.gravity_m_s2(),
        body.exterior.rotation_h,
        body.exterior.obliquity_deg,
        body.irradiance.flux_w_m2,
        body.irradiance.flux_w_m2 / 1361.0,
        body.interior.surface_ocean_gel_km(),
        body.interior.surface_pressure_bar()
    );
    if body.interior.has_atmosphere() {
        let rot = body.exterior.rotation_h;
        if let Some(s) = geoforge_surface::atmosphere::Substellar::for_body(&body) {
            println!("circulation: tidally locked; air rises under the star (0°N {:.0}°E) and returns along the ground from the night side; equatorial jet {:.1}", s.lon, s.superrotation);
        } else if body.interior.has_solid_surface() {
            let cells = geoforge_surface::atmosphere::Cells::for_body(rot, body.exterior.radius_km, 0.0, 0.0);
            println!("circulation: {} cells per hemisphere; eddy-driven jets: {}; Hadley cell to {:.0}°", cells.count(), cells.jets, cells.hadley);
        } else {
            let p = geoforge_surface::bands::Profile::for_body(&body, &geoforge_surface::bands::BandsParams::default());
            println!("bands: jets per hemisphere: {}; equatorial wind {:+.0} m/s", p.jets, p.wind(0.0));
        }
    }
    println!("magnetic field: {}", describe_field(&body));
    if body.interior.has_atmosphere() && body.interior.has_solid_surface() {
        println!("air: {}", describe_air(&body));
    }
    if !body.tides.is_empty() {
        println!("tides: {}", describe_tides(&body));
    }
    println!("calendar: {}", describe_calendar(&calendar_of(&body)));
    if let GridSpec::Equirect { width, height } = w.manifest.grid {
        println!("grid {width}×{height} equirectangular ({:.0} km cells at the equator)", 2.0 * std::f64::consts::PI * body.exterior.radius_km / f64::from(width));
    }
    Ok(())
}

/// "30 µT dipole tilted 7°; magnetopause 9.5 radii against a 2.0 nPa wind; auroras at 71°
/// magnetic latitude", or why there is none.
fn describe_field(body: &BodyContext) -> String {
    let f = MagneticField::of(body);
    let npa = |p: f64| if p < 0.1 { format!("{p:.1e} nPa") } else { format!("{p:.1} nPa") };
    let onto = if body.interior.has_atmosphere() { "air" } else { "surface" };
    if !f.dynamo {
        return match f.lifetime_gyr {
            Some(life) => format!(
                "none: the core dynamo ran {life:.1} Gyr and stopped {:.1} Gyr ago; the star's wind ({}) reaches the {onto}",
                body.epoch.system_age_gyr - life,
                npa(f.wind_npa)
            ),
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
fn calendar_of(body: &BodyContext) -> Calendar {
    let moons: Vec<(String, f64)> = body
        .tides
        .iter()
        .filter(|t| t.mass_kg < body.exterior.mass_kg)
        .map(|t| (t.name.clone(), geoforge_core::derived::orbital_period_days(t.distance_km / geoforge_core::body::AU_KM, body.exterior.mass_kg + t.mass_kg)))
        .collect();
    Calendar::new(body, &moons)
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
}

/// Positions on `day`, printed and optionally written as JSON.
fn positions(w: &World, day: f64, to: Option<&Path>) -> R<()> {
    let sys = w.load_system()?.ok_or("this world has no solar system (an archetype or hand-made body)")?;
    let body = w.load_body()?;
    let calendar = calendar_of(&body);
    let (year, local) = calendar.date(day);
    let unit = if calendar.solar_day_h.is_some() { "local day" } else { "day" };
    println!("day {day} (Earth days from the present): {}, year {year}, {unit} {local:.2}", body.name);
    let placed = geoforge_core::ephemeris::positions(&sys, day);
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
            let f = geoforge_core::ephemeris::lit_fraction(pos, focus, star);
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

fn print_weather(f: &pipe::weather::Forecast) {
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
            d.sky
        );
    }
}

/// "NNE" for 22.5°.
fn compass(deg: f64) -> &'static str {
    const POINTS: [&str; 16] = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE", "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"];
    POINTS[((deg.rem_euclid(360.0) / 22.5).round() as usize) % 16]
}

/// The sky on `day`, printed, and with `to` charted and catalogued.
fn sky(w: &World, day: f64, to: Option<&Path>, width: u32) -> R<()> {
    use pipe::sky::Kind;
    let sys = w.load_system()?.ok_or("this world has no solar system (an archetype or hand-made body)")?;
    let body = w.load_body()?;
    let calendar = calendar_of(&body);
    let s = pipe::sky::sky(&sys, w.manifest.focus, day, calendar.year_days)?;
    let (year, local) = calendar.date(day);
    println!("The sky of {} on day {day} (year {year}, day {local:.1}): {} naked-eye stars; {} {:.1} kpc from its centre", body.name, s.stars.len(), s.galaxy, s.radius_kpc);
    let radec = |o: &pipe::sky::Object| format!("RA {:>5.1}° Dec {:>+5.1}°", o.ra_deg, o.dec_deg);
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
                "  day {:>8.2} (year {y}, day {d:.1}): {:?} {:?} eclipse of {} by {}, magnitude {:.2}",
                e.day,
                e.coverage,
                e.kind,
                if e.kind == geoforge_core::ephemeris::EclipseKind::Solar { &e.star } else { &e.body },
                if e.kind == geoforge_core::ephemeris::EclipseKind::Solar { &e.body } else { &body.name },
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

fn print_census(c: &geoforge_cosmos::census::Census) {
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
