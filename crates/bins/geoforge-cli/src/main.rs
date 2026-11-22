//! `geoforge-v2`: make a world, generate it, look at it, edit any layer and regenerate.

use clap::{Parser, Subcommand};
use geoforge_core::grid::GridSpec;
use geoforge_core::product::{LayerId, ProductId};
use geoforge_core::body::{NodeKind, Star, StellarStage};
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
        return Ok((ProductId::ALL.to_vec(), View::ALL.to_vec()));
    }
    let (mut ps, mut vs) = (Vec::new(), Vec::new());
    for w in what {
        if let Some(p) = ProductId::from_name(w) {
            ps.push(p);
        } else if let Some(l) = LayerId::from_name(w) {
            ps.extend(l.products());
        } else if let Some(v) = View::from_name(w) {
            vs.push(v);
        } else {
            let views = View::ALL.map(View::name).join(", ");
            return Err(format!("unknown product, layer or view `{w}` (views: {views}; products: {})", ProductId::ALL.map(ProductId::name).join(", ")).into());
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
                let files = pipe::render_all(&w, &products, &w.root.join("renders"), &ProductId::ALL, &View::ALL)?;
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
        Cmd::Export { dir, what, to } => {
            let w = World::open(&dir)?;
            let to = to.unwrap_or_else(|| pipe::edit_folder(&w));
            let (ps, vs) = parse_targets(&what)?;
            let ps: Vec<ProductId> = ps.into_iter().filter(|p| w.has_product(*p)).collect();
            for (p, files) in pipe::export(&w, &ps, &to)? {
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
