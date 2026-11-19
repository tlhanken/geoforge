//! `geoforge-v2`: make a world, generate it, look at it, edit any layer and regenerate.

use clap::{Parser, Subcommand};
use geoforge_core::grid::GridSpec;
use geoforge_core::product::{LayerId, ProductId};
use geoforge_core::{BodyContext, Seed};
use geoforge_pipeline::{self as pipe, GenerateOptions, LayerAction, NewWorld, State};
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
    /// Create a world: a generated solar system and its focus body (or your own body file).
    New {
        dir: PathBuf,
        #[arg(long, default_value_t = 42)]
        seed: u64,
        /// Grid width in cells (height is half).
        #[arg(long, default_value_t = 2048)]
        width: u32,
        /// A body file (TOML or JSON BodyContext) instead of a generated system.
        #[arg(long)]
        body: Option<PathBuf>,
        /// Use an Earth-like reference body instead of a generated system (for calibration).
        #[arg(long, conflicts_with = "body")]
        earth_like: bool,
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
    /// Write edit views (PNG + world file + sidecar) of products for editing elsewhere.
    Export {
        dir: PathBuf,
        /// Product names, layer names, or "all".
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
    /// Write preview renders of stored products.
    Render {
        dir: PathBuf,
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

fn parse_products(what: &[String]) -> R<Vec<ProductId>> {
    if what.is_empty() || what.iter().any(|w| w == "all") {
        return Ok(ProductId::ALL.to_vec());
    }
    let mut out = Vec::new();
    for w in what {
        if let Some(p) = ProductId::from_name(w) {
            out.push(p);
        } else if let Some(l) = LayerId::from_name(w) {
            out.extend(l.products());
        } else {
            return Err(format!("unknown product or layer `{w}` (products: {})", ProductId::ALL.map(ProductId::name).join(", ")).into());
        }
    }
    out.sort();
    out.dedup();
    Ok(out)
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
        Cmd::New { dir, seed, width, body, earth_like } => {
            if seed > i64::MAX as u64 {
                return Err("seed must be at most 9223372036854775807".into());
            }
            let mut body = body.as_deref().map(read_body).transpose()?;
            if earth_like {
                body = Some(geoforge_core::Archetype::Earth.context(geoforge_core::Key(1), "Terra"));
            }
            let w = pipe::create_world(&dir, &NewWorld { seed: Seed(seed), grid: GridSpec::Equirect { width, height: width / 2 }, body })?;
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
                let files = pipe::render_all(&w, &products, &w.root.join("renders"), &[])?;
                println!("renders: {} files in {}", files.len(), w.root.join("renders").display());
            }
        }
        Cmd::Plan { dir } => {
            let w = World::open(&dir)?;
            let plan = pipe::plan(&w)?;
            for lp in &plan.layers {
                let action = match (&lp.not_applicable, lp.needs_run()) {
                    (Some(r), _) => format!("skip — {r}"),
                    (None, true) => "run".into(),
                    (None, false) => "reuse".into(),
                };
                println!("{:<10} {action}", lp.layer.name());
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
            let ps: Vec<ProductId> = parse_products(&what)?.into_iter().filter(|p| w.has_product(*p)).collect();
            for (p, files) in pipe::export(&w, &ps, &to)? {
                println!("{:<18} {}", p.name(), files.join(" "));
            }
            println!("edit views in {} — edit, then: geoforge-v2 import {} <product>", to.display(), dir.display());
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
            let only = if what.is_empty() { Vec::new() } else { parse_products(&what)? };
            let to = to.unwrap_or_else(|| w.root.join("renders"));
            for f in pipe::render_all(&w, &products, &to, &only)? {
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
        Cmd::Stats { dir } => {
            let w = World::open(&dir)?;
            let products = pipe::load_all(&w)?;
            let body = w.load_body()?;
            let grid = w.manifest.grid.build()?;
            let s = geoforge_pipeline::stats::summarize(&products, &grid, body.exterior.radius_km);
            println!("{}", serde_json::to_string_pretty(&s)?);
        }
    }
    Ok(())
}

fn show(w: &World) -> R<()> {
    let body = w.load_body()?;
    if let Some(sys) = w.load_system()? {
        println!("system {} (seed {})", sys.root.name, w.manifest.root_seed.0);
        for n in &sys.nodes {
            let a = n.elements.map_or(String::new(), |e| format!("a={:.3} AU e={:.3}", e.a_au, e.e));
            let kind = match &n.node {
                geoforge_core::body::NodeKind::Star(s) => format!("star {} {:.2} M☉ {:.0} K", s.spectral_class, s.mass_solar, s.temperature_k),
                geoforge_core::body::NodeKind::Body(b) => format!("{:?} {:.2} M⊕ {:.0} km", b.class, b.mass_kg / geoforge_core::body::EARTH_MASS_KG, b.radius_km),
                other => format!("{other:?}").split('(').next().unwrap_or("").split(' ').next().unwrap_or("").to_string(),
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
    if let GridSpec::Equirect { width, height } = w.manifest.grid {
        println!("grid {width}×{height} equirectangular ({:.0} km cells at the equator)", 2.0 * std::f64::consts::PI * body.exterior.radius_km / f64::from(width));
    }
    Ok(())
}
