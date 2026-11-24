//! `geoforge`: make a world, generate it, look at it, edit any layer and regenerate.

use clap::{Parser, Subcommand};
use geoforge_cosmos::GalaxyPreset;
use std::path::PathBuf;
use std::process::ExitCode;

mod almanac;
mod edit;
mod stars;
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
    /// Create a world: a body of a universe, made from the universe's seed and its address.
    ///
    /// By default the universe's default body: in its galaxy most like the Milky Way, the
    /// wet world of a Sun-like star near that galaxy's reference circle. --at makes any other
    /// galaxy's, system's or body's instead; in a directory that already holds the universe it
    /// is added there and becomes the focus. A real body on its own (--archetype) or your own
    /// body file (--body) make a world outside any universe.
    New {
        /// The universe (or world) directory.
        dir: PathBuf,
        /// The universe's seed: the same seed makes the same universe.
        #[arg(long, default_value_t = 42)]
        seed: u64,
        /// Where in the universe: an address, e.g. `42/h6:0,0,-1:17.0/f12:3,-4,0:5/4`, or
        /// `real` for the real universe (the Earth). A galaxy is taken to its default star, a
        /// system to its most Earth-like body. In a universe's directory the universe part may
        /// be left out.
        #[arg(long, conflicts_with_all = ["body", "archetype"])]
        at: Option<String>,
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
        /// Focus on another body of the system: its name (see `show`), or "<belt name> #<n>"
        /// for a belt's n-th largest member.
        #[arg(long, conflicts_with_all = ["body", "archetype"])]
        focus: Option<String>,
    },
    /// The bodies made in a universe's directory, or make one of them its focus.
    Focus {
        /// The universe directory.
        dir: PathBuf,
        /// The body to work on: its address (the universe part may be left out).
        at: Option<String>,
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
        /// The calendar dates are shown in: `world` (its own `calendar.json`, else
        /// twelfths), `twelfths`, `lunar` (the heaviest moon) or `lunar:<moon>`.
        #[arg(long, default_value = "world")]
        calendar: geoforge_almanac::calendar::Choice,
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
    /// Edit the galaxy's stars around the world: add one, move, change or remove one, add a
    /// cluster. Positions are in the focus's sky (RA and Dec in degrees, distance in light
    /// years); stars by identity (from `stars near` or `stars export`) or name.
    Star {
        #[command(subcommand)]
        action: StarCmd,
    },
    /// The stars around the world: the nearest of every brightness, or the sky's to a
    /// magnitude written out to edit in a spreadsheet and read back.
    Stars {
        #[command(subcommand)]
        action: StarsCmd,
    },
    /// The world's own constellations, drawn on its sky chart.
    Constellation {
        #[command(subcommand)]
        action: ConstellationCmd,
    },
    /// The galaxies around the world's own: its place in the universe, its satellites and
    /// neighbours, nearest first.
    Universe {
        dir: PathBuf,
        /// How far to look, Mpc.
        #[arg(long, default_value_t = 5.0)]
        radius_mpc: f64,
        /// Only galaxies brighter than this (apparent V from the home galaxy's centre).
        #[arg(long, default_value_t = 12.0)]
        limit: f64,
    },
    /// Go to a star of this world's sky: its system's most Earth-like body, made in the same
    /// universe (and made its focus). A world on its own makes a new one in --to, its galaxy
    /// and star edits carried over.
    Visit {
        /// The world whose sky holds the star.
        dir: PathBuf,
        /// The star: its identity or name.
        star: String,
        /// For a world on its own: the new world's directory.
        #[arg(long)]
        to: Option<PathBuf>,
    },
    /// Draw many systems from a galaxy and count them.
    ///
    /// Kinds of star, planets and wet worlds, with example seeds of each for `new`.
    Census {
        /// The galaxy (as for `new`): milky-way, sa … sd, sba … sbd, lenticular, elliptical,
        /// dwarf-spheroidal, dwarf, starburst.
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

#[derive(Subcommand)]
enum StarCmd {
    /// Add a star.
    Add {
        dir: PathBuf,
        #[arg(long)]
        name: String,
        #[arg(long, allow_hyphen_values = true)]
        ra: f64,
        #[arg(long, allow_hyphen_values = true)]
        dec: f64,
        #[arg(long)]
        distance_ly: f64,
        /// Initial mass, M☉.
        #[arg(long, default_value_t = 1.0)]
        mass: f64,
        #[arg(long, default_value_t = 4.6)]
        age_gyr: f64,
        /// [Fe/H], dex.
        #[arg(long, default_value_t = 0.0, allow_hyphen_values = true)]
        feh: f64,
        /// An unresolved companion's initial mass, M☉.
        #[arg(long)]
        companion: Option<f64>,
    },
    /// Move a star (any of RA, Dec, distance).
    Move {
        dir: PathBuf,
        star: String,
        #[arg(long, allow_hyphen_values = true)]
        ra: Option<f64>,
        #[arg(long, allow_hyphen_values = true)]
        dec: Option<f64>,
        #[arg(long)]
        distance_ly: Option<f64>,
    },
    /// Rename a star or make it something else.
    Change {
        dir: PathBuf,
        star: String,
        #[arg(long)]
        name: Option<String>,
        #[arg(long)]
        mass: Option<f64>,
        #[arg(long)]
        age_gyr: Option<f64>,
        #[arg(long, allow_hyphen_values = true)]
        feh: Option<f64>,
    },
    /// Take a star away.
    Remove { dir: PathBuf, star: String },
    /// Add a cluster of stars born together.
    Cluster {
        dir: PathBuf,
        #[arg(long)]
        name: String,
        #[arg(long, allow_hyphen_values = true)]
        ra: f64,
        #[arg(long, allow_hyphen_values = true)]
        dec: f64,
        #[arg(long)]
        distance_ly: f64,
        /// Plummer radius, pc.
        #[arg(long, default_value_t = 2.0)]
        radius_pc: f64,
        /// Stars of every brightness.
        #[arg(long, default_value_t = 500)]
        members: u32,
        #[arg(long, default_value_t = 0.1)]
        age_gyr: f64,
        #[arg(long, default_value_t = 0.0, allow_hyphen_values = true)]
        feh: f64,
    },
}

#[derive(Subcommand)]
enum StarsCmd {
    /// The systems within a distance, nearest first, faint ones and all.
    Near {
        dir: PathBuf,
        #[arg(long, default_value_t = 20.0)]
        radius_ly: f64,
    },
    /// Write the sky's stars to a magnitude (CSV, or JSON for a .json file).
    Export {
        dir: PathBuf,
        #[arg(long)]
        to: PathBuf,
        #[arg(long, default_value_t = 6.5)]
        limit: f64,
    },
    /// Read an edited export back: changed rows become changes, rows without an identity
    /// added stars.
    Import { dir: PathBuf, from: PathBuf },
}

#[derive(Subcommand)]
enum ConstellationCmd {
    /// Add (or replace) a constellation from lines between stars, each `A-B` (identities or
    /// names).
    Add {
        dir: PathBuf,
        name: String,
        #[arg(long = "line", required = true)]
        lines: Vec<String>,
    },
    Remove {
        dir: PathBuf,
        name: String,
    },
    List {
        dir: PathBuf,
    },
    /// Save the figures drawn from the brightest stars as the world's own, to rename or
    /// edit (the sky chart draws them anyway while the world has none).
    Auto {
        dir: PathBuf,
        /// Replace the world's own constellations.
        #[arg(long)]
        replace: bool,
    },
}

type R<T> = Result<T, Box<dyn std::error::Error>>;

fn run(cli: Cli) -> R<()> {
    use geoforge_store::World;
    match cli.cmd {
        Cmd::New { dir, seed, at, width, body, archetype, focus } => world::new_world(&dir, seed, at.as_deref(), width, body, archetype, focus.as_deref()),
        Cmd::Focus { dir, at } => world::focus(&dir, at.as_deref()),
        Cmd::Star { action } => match action {
            StarCmd::Add { dir, name, ra, dec, distance_ly, mass, age_gyr, feh, companion } => stars::add(&dir, &name, ra, dec, distance_ly, mass, age_gyr, feh, companion),
            StarCmd::Move { dir, star, ra, dec, distance_ly } => stars::move_star(&dir, &star, ra, dec, distance_ly),
            StarCmd::Change { dir, star, name, mass, age_gyr, feh } => stars::change(&dir, &star, name, mass, age_gyr, feh),
            StarCmd::Remove { dir, star } => stars::remove(&dir, &star),
            StarCmd::Cluster { dir, name, ra, dec, distance_ly, radius_pc, members, age_gyr, feh } => stars::cluster(&dir, &name, ra, dec, distance_ly, radius_pc, members, age_gyr, feh),
        },
        Cmd::Stars { action } => match action {
            StarsCmd::Near { dir, radius_ly } => stars::near(&dir, radius_ly),
            StarsCmd::Export { dir, to, limit } => stars::export(&dir, &to, limit),
            StarsCmd::Import { dir, from } => stars::import(&dir, &from),
        },
        Cmd::Constellation { action } => match action {
            ConstellationCmd::Add { dir, name, lines } => stars::constellation_add(&dir, &name, &lines),
            ConstellationCmd::Remove { dir, name } => stars::constellation_remove(&dir, &name),
            ConstellationCmd::List { dir } => stars::constellation_list(&dir),
            ConstellationCmd::Auto { dir, replace } => stars::constellation_auto(&dir, replace),
        },
        Cmd::Universe { dir, radius_mpc, limit } => stars::universe(&dir, radius_mpc, limit),
        Cmd::Visit { dir, star, to } => stars::visit(&dir, &star, to.as_deref()),
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
        Cmd::Almanac { dir, lat, lon, year, calendar, json } => almanac::almanac(&dir, lat, lon, year, &calendar, json),
        Cmd::WeatherMap { dir, day, to, width } => almanac::weather_map(&dir, day, to, width),
        Cmd::Seed { dir, product, owner } => world::seed(&dir, &product, owner),
        Cmd::Stats { dir } => world::stats(&dir),
        Cmd::Report { dir, to, places, day, days, width } => world::report(&dir, to, &places, day, days, width),
        Cmd::Compare { real, models, labels, probes, sources, to, width } => world::compare(&real, &models, &labels, &probes, &sources, to, width),
    }
}
