//! `geoforge`: make a world, generate it, look at it, edit any layer and regenerate.

use clap::{Parser, Subcommand};
use geoforge_cosmos::GalaxyPreset;
use std::path::PathBuf;
use std::process::ExitCode;

mod almanac;
mod edit;
mod system;
mod world;

use world::parse_galaxy;

#[derive(Parser)]
#[command(
    name = "geoforge",
    version,
    about = "Procedural worlds you can edit at any layer",
    long_about = "Procedural worlds you can edit at any layer.\n\nA world is a directory: a solar system drawn from a galaxy, and its focus body generated layer by layer (tectonics, geology, elevation, circulation, climate, erosion, hydrology, biomes, resources). Any layer can be exported to images, edited or redrawn, and imported back; downstream layers follow. Places and days can then be asked about: what is there, its weather, its sky.",
    after_help = "Getting started:\n  geoforge new worlds/demo --seed 42        make a world (a system and its most Earth-like body)\n  geoforge generate worlds/demo --render    build every layer, with preview images in worlds/demo/renders\n  geoforge show worlds/demo                 the system and the world\n  geoforge place worlds/demo --lat 45 --lon 10     what is at a place\n  geoforge weather worlds/demo --lat 45 --lon 10   two weeks of weather there\n\nEditing:\n  geoforge export worlds/demo rough_elevation      write edit images to worlds/demo/edit\n  geoforge import worlds/demo rough_elevation      bring your edit back as a patch\n  geoforge generate worlds/demo                    the layers below follow\n\nRun `geoforge help <command>` for a command's options. Guide: docs/GETTING_STARTED.md"
)]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Create a world.
    ///
    /// A solar system drawn from a galaxy, with its most Earth-like body as the focus; or the
    /// habitable recipe (--habitable), a real body (--archetype) or your own body file (--body).
    New {
        /// The world directory.
        dir: PathBuf,
        /// The world's seed: the same seed makes the same world.
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
    /// Build every layer that is missing or stale.
    Generate {
        /// The world directory.
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
    Plan {
        /// The world directory.
        dir: PathBuf,
    },
    /// Where every product comes from (generated, patched, authored) and whether it is fresh.
    Status {
        /// The world directory.
        dir: PathBuf,
    },
    /// Write layers as images to edit elsewhere.
    ///
    /// Products become edit views (PNG + world file + sidecar) that GIMP, Krita, QGIS or
    /// Blender open; derived views (koppen, hazards, …) are written as PNG + legend + world
    /// file for reference.
    Export {
        /// The world directory.
        dir: PathBuf,
        /// Product names, layer names, view names, or "all".
        what: Vec<String>,
        /// Where to write (a folder, or a file for single outputs).
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
    /// Bring an edited image back.
    ///
    /// By default the edit becomes a patch on the generated product (only what you changed,
    /// kept when the layer regenerates); with --authored it replaces the product.
    Import {
        /// The world directory.
        dir: PathBuf,
        /// The product (run `geoforge export --help` for names).
        product: String,
        /// Folder holding the edited views (default: the world's edit folder).
        #[arg(long)]
        from: Option<PathBuf>,
        /// Replace the product with your data instead of patching it.
        #[arg(long)]
        authored: bool,
    },
    /// Discard a product's edits; it is generated again.
    Reset {
        /// The world directory.
        dir: PathBuf,
        /// The product (see `export` for names).
        product: String,
    },
    /// Write preview images of products and derived views.
    Render {
        /// The world directory.
        dir: PathBuf,
        /// Product names, layer names, view names, or "all" (the default).
        what: Vec<String>,
        /// Where to write (a folder, or a file for single outputs).
        #[arg(long)]
        to: Option<PathBuf>,
    },
    /// Describe the system and the focus: stars, planets, air, water, field, tides, calendar.
    Show {
        /// The world directory.
        dir: PathBuf,
    },
    /// What is at a place: climate, soil, deposits and hazards.
    ///
    /// Climate, soil and timber at the centre; the deposits likely within reach (and how many
    /// a local map holds); and each hazard's yearly rate and chance over a stretch of days.
    Place {
        /// The world directory.
        dir: PathBuf,
        /// Latitude, degrees (south negative).
        #[arg(long, allow_hyphen_values = true)]
        lat: f64,
        /// Longitude, degrees (west negative).
        #[arg(long, allow_hyphen_values = true)]
        lon: f64,
        /// Radius of the area, km.
        #[arg(long, default_value_t = 100.0)]
        radius_km: f64,
        /// The stretch of days: from this day (Earth days from the present), for this many
        /// (default: one of the world's years).
        #[arg(long, default_value_t = 0.0)]
        from_day: f64,
        /// How many days.
        #[arg(long)]
        days: Option<f64>,
        /// Print JSON instead.
        #[arg(long)]
        json: bool,
    },
    /// Where every star, planet and moon is on a day.
    ///
    /// Days are Earth days from the present; also the focus's date on its own calendar and
    /// how lit its moons are.
    Positions {
        /// The world directory.
        dir: PathBuf,
        /// Earth days from the present (the focus's northern spring equinox).
        #[arg(long, default_value_t = 0.0)]
        day: f64,
        /// Also write them as JSON.
        #[arg(long)]
        to: Option<PathBuf>,
    },
    /// The night sky on a day: stars, planets, moons and the year's eclipses.
    ///
    /// With --to, also a star chart (sky.png) and a catalogue of every naked-eye star
    /// (sky.json).
    Sky {
        /// The world directory.
        dir: PathBuf,
        /// Earth days from the present (the focus's northern spring equinox).
        #[arg(long, default_value_t = 0.0)]
        day: f64,
        /// Where to write (a folder, or a file for single outputs).
        #[arg(long)]
        to: Option<PathBuf>,
        /// Chart width in pixels (half as tall).
        #[arg(long, default_value_t = 2048)]
        width: u32,
    },
    /// Day-by-day weather at a place.
    ///
    /// Drawn from its climate, the season of the date and the storms passing; with the day's
    /// phenomena (thunderstorms, tornadoes, quakes, auroras, …). The same every time for the
    /// same date. Days are Earth days from the present.
    Weather {
        /// The world directory.
        dir: PathBuf,
        /// Latitude, degrees (south negative).
        #[arg(long, allow_hyphen_values = true)]
        lat: f64,
        /// Longitude, degrees (west negative).
        #[arg(long, allow_hyphen_values = true)]
        lon: f64,
        /// Earth days from the present (the focus's northern spring equinox).
        #[arg(long, default_value_t = 0.0)]
        day: f64,
        /// How many days (local days on a world with days).
        #[arg(long, default_value_t = 14)]
        days: usize,
        /// Print JSON instead of a table.
        #[arg(long)]
        json: bool,
    },
    /// A place's year day by day: sunrise and sunset, the moons, the weather and the day's
    /// events; then the year's months of rain, drought index (SPI-3) and soil water, and the
    /// place's droughts over 30 years.
    Almanac {
        /// The world directory.
        dir: PathBuf,
        /// Latitude, degrees (south negative).
        #[arg(long, allow_hyphen_values = true)]
        lat: f64,
        /// Longitude, degrees (west negative).
        #[arg(long, allow_hyphen_values = true)]
        lon: f64,
        /// The local year (0 the present one).
        #[arg(long, default_value_t = 0, allow_hyphen_values = true)]
        year: i64,
        /// Print JSON instead.
        #[arg(long)]
        json: bool,
    },
    /// The weather map of a day: isobars, fronts, rain and storm centres.
    WeatherMap {
        /// The world directory.
        dir: PathBuf,
        /// Earth days from the present (the focus's northern spring equinox).
        #[arg(long, default_value_t = 0.0)]
        day: f64,
        /// Where to write (a folder, or a file for single outputs).
        #[arg(long)]
        to: Option<PathBuf>,
        /// Image width in pixels.
        #[arg(long, default_value_t = 1440)]
        width: u32,
    },
    /// Print the seed a product is generated from.
    Seed {
        /// The world directory.
        dir: PathBuf,
        /// The product (run `geoforge export --help` for names).
        product: String,
        /// Owner key (defaults to the focus body).
        #[arg(long)]
        owner: Option<u64>,
    },
    /// Summary statistics of the products, with Earth's for comparison.
    Stats {
        /// The world directory.
        dir: PathBuf,
    },
    /// The A/B page: every product an observed world holds as authored data, beside the same
    /// product in each model world, with one renderer and key and the numbers between them.
    Compare {
        /// The observed world (its authored products are the observations).
        real: PathBuf,
        /// The model worlds.
        #[arg(required = true)]
        models: Vec<PathBuf>,
        /// A label for each model, in order; default its folder's name.
        #[arg(long = "label")]
        labels: Vec<String>,
        /// A place or box whose values are listed, "name:lat,lon" or
        /// "name:lat0..lat1,lon0..lon1" (repeat for more).
        #[arg(long = "probe")]
        probes: Vec<String>,
        /// Where an observed product came from, "product=text" (repeat for more), for its
        /// caption.
        #[arg(long = "source")]
        sources: Vec<String>,
        /// Folder to write (index.html and img/); default <real>/compare.
        #[arg(long)]
        to: Option<PathBuf>,
        /// Width of the maps, px.
        #[arg(long, default_value_t = 1200)]
        width: u32,
    },
    /// A report page to take the world's pulse: every map with its key, the numbers against
    /// Earth's, places and their weather, the weather map and the night sky.
    Report {
        /// The world directory.
        dir: PathBuf,
        /// Folder to write (index.html and img/); default <dir>/report.
        #[arg(long)]
        to: Option<PathBuf>,
        /// A place to report on, "name:lat,lon" (repeat for more); default a temperate, a
        /// tropical and a dry one.
        #[arg(long = "place")]
        places: Vec<String>,
        /// Day of the weather and sky, Earth days from the present.
        #[arg(long, default_value_t = 30.0)]
        day: f64,
        /// Days of weather at each place.
        #[arg(long, default_value_t = 14)]
        days: usize,
        /// Width of the maps, px.
        #[arg(long, default_value_t = 1200)]
        width: u32,
    },
    /// Draw many systems from a galaxy and count them.
    ///
    /// Kinds of star, planets and wet worlds, with example seeds of each for `new`.
    Census {
        /// The galaxy: milky-way, elliptical, dwarf, starburst.
        #[arg(long, default_value = "milky-way", value_parser = parse_galaxy)]
        galaxy: GalaxyPreset,
        /// How many seeds, and the first.
        #[arg(long, default_value_t = 10_000)]
        count: u64,
        /// Folder holding the edited views (default: the world's edit folder).
        #[arg(long, default_value_t = 0)]
        from: u64,
    },
}

fn main() -> ExitCode {
    // Piped into `head` and the like, stop quietly when the reader goes away.
    #[cfg(unix)]
    // SAFETY: restoring the default SIGPIPE action before any other thread starts.
    unsafe {
        libc::signal(libc::SIGPIPE, libc::SIG_DFL);
    }
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

type R<T> = Result<T, Box<dyn std::error::Error>>;

fn run(cli: Cli) -> R<()> {
    use geoforge_store::World;
    match cli.cmd {
        Cmd::New { dir, seed, galaxy, habitable, width, body, archetype, focus } => world::new_world(&dir, seed, galaxy, habitable, width, body, archetype, focus),
        Cmd::Generate { dir, force, until, from, render } => world::generate(&dir, force, until, from, render),
        Cmd::Plan { dir } => world::plan(&dir),
        Cmd::Status { dir } => world::status(&dir),
        Cmd::Export { dir, what, to, scale, caps } => edit::export(&dir, &what, to, scale, caps),
        Cmd::Import { dir, product, from, authored } => edit::import(&dir, &product, from, authored),
        Cmd::Reset { dir, product } => edit::reset(&dir, &product),
        Cmd::Render { dir, what, to } => edit::render(&dir, &what, to),
        Cmd::Show { dir } => system::show(&World::open(&dir)?),
        Cmd::Positions { dir, day, to } => system::positions(&World::open(&dir)?, day, to.as_deref()),
        Cmd::Census { galaxy, count, from } => {
            system::print_census(&geoforge_cosmos::census::census(&galaxy.galaxy(), from..from.saturating_add(count)));
            Ok(())
        }
        Cmd::Sky { dir, day, to, width } => almanac::sky(&World::open(&dir)?, day, to.as_deref(), width),
        Cmd::Place { dir, lat, lon, radius_km, from_day, days, json } => almanac::place(&dir, lat, lon, radius_km, from_day, days, json),
        Cmd::Weather { dir, lat, lon, day, days, json } => almanac::weather(&dir, lat, lon, day, days, json),
        Cmd::Almanac { dir, lat, lon, year, json } => almanac::almanac(&dir, lat, lon, year, json),
        Cmd::WeatherMap { dir, day, to, width } => almanac::weather_map(&dir, day, to, width),
        Cmd::Seed { dir, product, owner } => world::seed(&dir, &product, owner),
        Cmd::Stats { dir } => world::stats(&dir),
        Cmd::Report { dir, to, places, day, days, width } => world::report(&dir, to, &places, day, days, width),
        Cmd::Compare { real, models, labels, probes, sources, to, width } => world::compare(&real, &models, &labels, &probes, &sources, to, width),
    }
}
