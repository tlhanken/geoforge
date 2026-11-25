//! `geoforge-viewer DIR`: look around a universe in a browser.

use clap::Parser;
use std::path::PathBuf;

/// Look around a universe in a browser: zoom from its galaxies to a galaxy's stars, a system
/// and a body's maps, search, and make or generate any body you find.
#[derive(Parser)]
#[command(name = "geoforge-viewer", version)]
struct Args {
    /// The universe directory (`geoforge new` makes one).
    dir: PathBuf,
    /// Where to listen: 127.0.0.1 for this computer only, 0.0.0.0 to reach it from another
    /// device on your network.
    #[arg(long, default_value = "127.0.0.1")]
    host: String,
    #[arg(long, default_value_t = 8137)]
    port: u16,
}

fn main() {
    let a = Args::parse();
    let run = || -> geoforge_viewer::Result<()> {
        let viewer = geoforge_viewer::Viewer::open(&a.dir)?;
        println!("geoforge viewer: http://{}:{}/ (Ctrl-C to stop)", if a.host == "0.0.0.0" { "localhost" } else { &a.host }, a.port);
        geoforge_viewer::server::serve(viewer, &format!("{}:{}", a.host, a.port))
    };
    if let Err(e) = run() {
        eprintln!("geoforge-viewer: {e}");
        std::process::exit(1);
    }
}
