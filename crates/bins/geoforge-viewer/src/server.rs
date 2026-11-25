//! The viewer behind a local web server: the page and its script, and an API by address.
//!
//! ```text
//! GET  /                                   the page
//! GET  /api/start                          the universe, its focus body and home galaxy
//! GET  /api/scene?at=ADDRESS[&stars=1]     a scene (scene::Scene as JSON)
//! GET  /api/search?q=TEXT                  hits (search::Hit)
//! GET  /api/habitable?at=ADDRESS          the systems around it with a habitable or
//!                                         temperate world: [{address, habitability}]
//! GET  /api/texture?at=ADDRESS&layer=NAME&width=W    a map as PNG
//! POST /api/make?at=ADDRESS                make the body's world (and focus it)
//! POST /api/generate?at=ADDRESS            generate its layers, in the background
//! ```

use crate::{scene, Result, ViewError, Viewer};
use geoforge_cosmos::worlds::Habitability;
use std::sync::Arc;
use tiny_http::{Header, Method, Request, Response, Server};

const PAGE: &str = include_str!("../web/index.html");
const SCRIPT: &str = include_str!("../web/app.js");
const THREE: &str = include_str!("../web/three.js");

/// Grid width of a body made from the viewer, cells.
const MAKE_WIDTH: u32 = 1024;
/// Requests handled at once.
const WORKERS: usize = 4;

/// Serve `viewer` at `address` (e.g. `127.0.0.1:8137`) until the process ends.
///
/// # Errors
/// If the address cannot be listened on.
pub fn serve(viewer: Viewer, address: &str) -> Result<()> {
    let server = Arc::new(Server::http(address).map_err(|e| ViewError(format!("cannot listen on {address}: {e}")))?);
    let viewer = Arc::new(viewer);
    let workers: Vec<_> = (0..WORKERS)
        .map(|_| {
            let (server, viewer) = (Arc::clone(&server), Arc::clone(&viewer));
            std::thread::spawn(move || {
                while let Ok(rq) = server.recv() {
                    handle(&viewer, rq);
                }
            })
        })
        .collect();
    for w in workers {
        let _ = w.join();
    }
    Ok(())
}

fn handle(viewer: &Arc<Viewer>, rq: Request) {
    let url = rq.url().to_string();
    let (path, query) = url.split_once('?').unwrap_or((&url, ""));
    let arg = |k: &str| query.split('&').filter_map(|kv| kv.split_once('=')).find(|(key, _)| *key == k).map(|(_, v)| decode(v));
    let answer = match (rq.method(), path) {
        (Method::Get, "/") => Ok(text(PAGE, "text/html; charset=utf-8")),
        (Method::Get, "/app.js") => Ok(text(SCRIPT, "text/javascript; charset=utf-8")),
        (Method::Get, "/three.js") => Ok(text(THREE, "text/javascript; charset=utf-8")),
        (Method::Get, "/api/start") => {
            let (universe, focus) = viewer.start();
            Ok(json(&serde_json::json!({ "universe": universe.to_string(), "focus": focus.map(|f| f.to_string()), "home": viewer.home().to_string() })))
        }
        (Method::Get, "/api/scene") => at(viewer, arg("at")).and_then(|a| {
            let mut v = serde_json::to_value(viewer.scene(&a, arg("stars").is_some())?).map_err(|e| ViewError(e.to_string()))?;
            v["crumbs"] = serde_json::to_value(viewer.crumbs(&a)?).map_err(|e| ViewError(e.to_string()))?;
            Ok(json(&v))
        }),
        (Method::Get, "/api/habitable") => at(viewer, arg("at")).and_then(|a| {
            let scene::Scene::Stars(st) = viewer.scene(&a, true)? else { return Err(ViewError("no stars there".into())) };
            let systems: Vec<geoforge_core::Address> = st.stars.iter().filter_map(|s| s.address.parse().ok()).collect();
            let found: Vec<serde_json::Value> = systems
                .iter()
                .zip(viewer.habitability(&systems))
                .filter(|(_, h)| *h != Habitability::None)
                .map(|(a, h)| serde_json::json!({ "address": a.to_string(), "habitability": if h == Habitability::Temperate { "temperate" } else { "habitable" } }))
                .collect();
            Ok(json(&found))
        }),
        (Method::Get, "/api/search") => viewer.search(&arg("q").unwrap_or_default()).map(|h| json(&h)),
        (Method::Get, "/api/texture") => at(viewer, arg("at")).and_then(|a| {
            let width = arg("width").and_then(|w| w.parse().ok()).unwrap_or(1024);
            viewer.texture(&a, &arg("layer").unwrap_or_else(|| "appearance".into()), width).map(|png| Response::from_data(png).with_header(header("Content-Type", "image/png")))
        }),
        (Method::Post, "/api/make") => at(viewer, arg("at")).and_then(|a| viewer.make(&a, MAKE_WIDTH)).map(|w| json(&serde_json::json!({ "address": w.manifest.address.map(|a| a.to_string()) }))),
        (Method::Post, "/api/generate") => at(viewer, arg("at")).map(|a| {
            let v = Arc::clone(viewer);
            v.jobs.lock().expect("jobs").insert(a, "queued".into());
            std::thread::spawn(move || {
                let _ = v.generate(&a, MAKE_WIDTH);
            });
            json(&serde_json::json!({ "address": a.to_string(), "job": "queued" }))
        }),
        _ => Err(ViewError(format!("nothing at {path}"))),
    };
    let response = answer.unwrap_or_else(|e| Response::from_string(serde_json::json!({ "error": e.0 }).to_string()).with_status_code(400).with_header(header("Content-Type", "application/json")));
    let _ = rq.respond(response);
}

fn at(viewer: &Viewer, text: Option<String>) -> Result<geoforge_core::Address> {
    viewer.parse(&text.ok_or_else(|| ViewError("give an address: at=…".into()))?)
}

fn text(s: &str, kind: &str) -> Response<std::io::Cursor<Vec<u8>>> {
    Response::from_string(s).with_header(header("Content-Type", kind))
}

fn json<T: serde::Serialize>(v: &T) -> Response<std::io::Cursor<Vec<u8>>> {
    Response::from_string(serde_json::to_string(v).unwrap_or_default()).with_header(header("Content-Type", "application/json"))
}

fn header(k: &str, v: &str) -> Header {
    Header::from_bytes(k.as_bytes(), v.as_bytes()).expect("a valid header")
}

/// A URL query value: `+` and `%XX` decoded.
fn decode(v: &str) -> String {
    let bytes = v.as_bytes();
    let mut out = Vec::with_capacity(bytes.len());
    let mut i = 0;
    while i < bytes.len() {
        match bytes[i] {
            b'+' => out.push(b' '),
            b'%' if i + 2 < bytes.len() => {
                if let Ok(b) = u8::from_str_radix(std::str::from_utf8(&bytes[i + 1..i + 3]).unwrap_or("zz"), 16) {
                    out.push(b);
                    i += 2;
                } else {
                    out.push(b'%');
                }
            }
            b => out.push(b),
        }
        i += 1;
    }
    String::from_utf8_lossy(&out).into_owned()
}

#[cfg(test)]
mod tests {
    #[test]
    fn queries_decode() {
        assert_eq!(super::decode("7%2Fh5%3A0%2C-1%2C0%3A37.0+x"), "7/h5:0,-1,0:37.0 x");
        assert_eq!(super::decode("100%"), "100%");
    }
}
