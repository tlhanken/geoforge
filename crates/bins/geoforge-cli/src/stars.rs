//! The galaxy's stars around a world: `star` (add, move, change, remove, cluster), `stars`
//! (near, export, import), `constellation` (add, remove, list) and `visit` (a new world at a
//! star). Positions are given in the focus's sky: right ascension and declination (degrees)
//! and distance (light-years); stars by identity or by name.

use crate::R;
use geoforge_almanac::sky::{starfield, Frame};
use geoforge_core::galactic::{ClusterEdit, Constellation, StarId, StarKind, StarRecord};
use geoforge_cosmos::sky::{neighbourhood, NAKED_EYE};
use geoforge_cosmos::Starfield;
use geoforge_store::World;
use std::path::Path;

const LY_PER_PC: f64 = geoforge_cosmos::sky::LY_PER_PC;

/// A world with its stars and sky frame.
struct Here {
    world: World,
    field: Starfield,
    frame: Frame,
}

impl Here {
    fn open(dir: &Path) -> R<Self> {
        let world = World::open(dir)?;
        let system = world.load_system()?.ok_or("this world has no solar system (a hand-made body)")?;
        let field = starfield(&world)?;
        let frame = Frame::of(&system, world.manifest.focus, &field)?;
        Ok(Self { world, field, frame })
    }

    /// The star named or identified by `text`: an identity, or the name of an edited star,
    /// a naked-eye star or one within 50 light-years.
    fn find(&self, text: &str) -> R<StarId> {
        if let Ok(id) = text.parse::<StarId>() {
            return Ok(id);
        }
        let e = self.field.edits();
        if let Some(r) = e.add.iter().chain(&e.change).find(|r| r.name.as_deref().is_some_and(|n| n.eq_ignore_ascii_case(text))) {
            return Ok(r.id);
        }
        let named = |id: &StarId| self.field.name_of(id).eq_ignore_ascii_case(text);
        let n = neighbourhood(&self.field, &self.frame.observer, NAKED_EYE);
        if let Some(s) = n.stars.iter().find(|s| named(&s.id)) {
            return Ok(s.id);
        }
        let near = self.field.near(self.frame.observer.position_pc, 50.0 / LY_PER_PC, 40.0);
        near.iter()
            .find(|s| named(&s.id))
            .map(|s| s.id)
            .ok_or_else(|| format!("no star `{text}` (give its identity, or the name of an edited star, a naked-eye star or one within 50 light-years)").into())
    }

    fn save(&self, edits: &geoforge_core::galactic::StarEdits) -> R<()> {
        self.world.save_star_edits(edits)?;
        Ok(())
    }

    fn describe(&self, id: &StarId) -> String {
        match self.field.star(id) {
            Some(s) => {
                let (ra, dec, d) = self.frame.radec_of(s.position_pc);
                format!("{} ({id}): {}, RA {ra:.2}° Dec {dec:+.2}°, {:.1} ly", self.field.name_of(id), s.star.spectral_class, d * LY_PER_PC)
            }
            None => format!("{id}: gone"),
        }
    }
}

/// A star's kind from its mass, age and metallicity.
fn kind(mass: f64, age_gyr: f64, metallicity_dex: f64, companion_mass: Option<f64>) -> R<StarKind> {
    if !(0.08..=150.0).contains(&mass) || age_gyr < 0.0 {
        return Err("a star's mass is 0.08 to 150 M☉, its age at least 0".into());
    }
    Ok(StarKind { initial_mass: mass, age_gyr, metallicity_dex, companion_mass })
}

#[allow(clippy::too_many_arguments)]
pub fn add(dir: &Path, name: &str, ra: f64, dec: f64, distance_ly: f64, mass: f64, age_gyr: f64, feh: f64, companion: Option<f64>) -> R<()> {
    let h = Here::open(dir)?;
    let mut e = h.field.edits().clone();
    let id = e.next_added();
    e.put(StarRecord { id, name: Some(name.into()), position_pc: h.frame.galactic_of(ra, dec, distance_ly / LY_PER_PC), kind: kind(mass, age_gyr, feh, companion)? });
    h.save(&e)?;
    println!("added {name} as {id}");
    Ok(())
}

/// Move `star`: any of RA, Dec and distance (the rest kept).
pub fn move_star(dir: &Path, star: &str, ra: Option<f64>, dec: Option<f64>, distance_ly: Option<f64>) -> R<()> {
    let h = Here::open(dir)?;
    let id = h.find(star)?;
    let s = h.field.star(&id).ok_or("that star is gone")?;
    let (r0, d0, dist0) = h.frame.radec_of(s.position_pc);
    let to = h.frame.galactic_of(ra.unwrap_or(r0), dec.unwrap_or(d0), distance_ly.map_or(dist0, |l| l / LY_PER_PC));
    let mut e = h.field.edits().clone();
    let name = e.record(&id).and_then(|r| r.name.clone());
    e.put(StarRecord { id, name, position_pc: to, kind: s.kind });
    h.save(&e)?;
    println!("moved {}", Here::open(dir)?.describe(&id));
    Ok(())
}

/// Change `star`'s name or what it is.
pub fn change(dir: &Path, star: &str, name: Option<String>, mass: Option<f64>, age_gyr: Option<f64>, feh: Option<f64>) -> R<()> {
    let h = Here::open(dir)?;
    let id = h.find(star)?;
    let s = h.field.star(&id).ok_or("that star is gone")?;
    let k = kind(mass.unwrap_or(s.kind.initial_mass), age_gyr.unwrap_or(s.kind.age_gyr), feh.unwrap_or(s.kind.metallicity_dex), s.kind.companion_mass)?;
    let mut e = h.field.edits().clone();
    let name = name.or_else(|| e.record(&id).and_then(|r| r.name.clone())).or_else(|| Some(h.field.name_of(&id)));
    e.put(StarRecord { id, name, position_pc: s.position_pc, kind: k });
    h.save(&e)?;
    println!("changed {}", Here::open(dir)?.describe(&id));
    Ok(())
}

pub fn remove(dir: &Path, star: &str) -> R<()> {
    let h = Here::open(dir)?;
    let id = h.find(star)?;
    let mut e = h.field.edits().clone();
    e.take_away(id);
    h.save(&e)?;
    println!("removed {id}");
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub fn cluster(dir: &Path, name: &str, ra: f64, dec: f64, distance_ly: f64, radius_pc: f64, members: u32, age_gyr: f64, feh: f64) -> R<()> {
    let h = Here::open(dir)?;
    let mut e = h.field.edits().clone();
    let n = e.next_cluster();
    e.clusters.push(ClusterEdit { n, name: name.into(), centre_pc: h.frame.galactic_of(ra, dec, distance_ly / LY_PER_PC), radius_pc, members, age_gyr, metallicity_dex: feh });
    h.save(&e)?;
    println!("added the cluster {name} ({members} stars within ~{radius_pc} pc), its members c{n}:<n>");
    Ok(())
}

/// The stars within `radius_ly`, nearest first: the census, faint stars and all.
pub fn near(dir: &Path, radius_ly: f64) -> R<()> {
    let h = Here::open(dir)?;
    let here = h.frame.observer.position_pc;
    let mut stars = h.field.near(here, radius_ly / LY_PER_PC, 40.0);
    stars.retain(|s| Some(s.id) != h.frame.observer.exclude && s.distance_pc(here) > h.frame.observer.clear_pc);
    stars.sort_by(|a, b| a.distance_pc(here).total_cmp(&b.distance_pc(here)));
    println!("{} systems within {radius_ly} light-years:", stars.len());
    for s in &stars {
        let (ra, dec, d) = h.frame.radec_of(s.position_pc);
        let pair = if s.kind.companion_mass.is_some() { " + companion" } else { "" };
        let v = s.absolute_v + 5.0 * (d / 10.0).log10();
        println!("  {:>7.2} ly  {:<14} {:<8}{pair:<12} V {v:>5.1}  RA {ra:>6.2}° Dec {dec:>+6.2}°  {}", d * LY_PER_PC, h.field.name_of(&s.id), s.star.spectral_class, s.id);
    }
    Ok(())
}

const COLUMNS: [&str; 12] = ["id", "name", "ra_deg", "dec_deg", "distance_ly", "v", "absolute_v", "spectral_class", "initial_mass", "age_gyr", "metallicity_dex", "companion_mass"];

/// The sky's stars to magnitude `limit`, with identities, as CSV (or JSON for a `.json`
/// path), to edit in a spreadsheet and import.
pub fn export(dir: &Path, to: &Path, limit: f64) -> R<()> {
    let h = Here::open(dir)?;
    let n = neighbourhood(&h.field, &h.frame.observer, limit);
    let rows: Vec<Vec<String>> = n
        .stars
        .iter()
        .map(|s| {
            let (ra, dec, d) = h.frame.radec_of(s.position_pc);
            let k = h.field.star(&s.id).map_or(StarKind { initial_mass: s.star.initial_mass(), age_gyr: s.star.age_gyr, metallicity_dex: s.star.metallicity_dex, companion_mass: None }, |g| g.kind);
            vec![
                s.id.to_string(),
                h.field.name_of(&s.id),
                format!("{ra:.5}"),
                format!("{dec:.5}"),
                format!("{:.4}", d * LY_PER_PC),
                format!("{:.3}", s.v),
                format!("{:.3}", s.absolute_v),
                s.star.spectral_class.clone(),
                format!("{:.5}", k.initial_mass),
                format!("{:.5}", k.age_gyr),
                format!("{:.4}", k.metallicity_dex),
                k.companion_mass.map_or(String::new(), |m| format!("{m:.5}")),
            ]
        })
        .collect();
    if to.extension().is_some_and(|e| e == "json") {
        let objs: Vec<serde_json::Map<String, serde_json::Value>> =
            rows.iter().map(|r| COLUMNS.iter().zip(r).map(|(c, v)| ((*c).to_string(), serde_json::Value::String(v.clone()))).collect()).collect();
        std::fs::write(to, serde_json::to_string_pretty(&objs)?)?;
    } else {
        let mut out = COLUMNS.join(",") + "\n";
        for r in &rows {
            out += &r.iter().map(|v| csv_field(v)).collect::<Vec<_>>().join(",");
            out.push('\n');
        }
        std::fs::write(to, out)?;
    }
    println!("{} stars to magnitude {:.1} written to {}", rows.len(), n.limit, to.display());
    Ok(())
}

/// Read an export back: a row whose star moved, was renamed or made something else becomes a
/// change; a row without an identity, an added star. Rows left out are not removed (use
/// `geoforge star remove`).
pub fn import(dir: &Path, from: &Path) -> R<()> {
    let h = Here::open(dir)?;
    let text = std::fs::read_to_string(from)?;
    let rows: Vec<std::collections::HashMap<String, String>> = if from.extension().is_some_and(|e| e == "json") {
        serde_json::from_str::<Vec<std::collections::HashMap<String, serde_json::Value>>>(&text)?
            .into_iter()
            .map(|m| m.into_iter().map(|(k, v)| (k, v.as_str().map_or_else(|| v.to_string(), str::to_string))).collect())
            .collect()
    } else {
        let mut lines = text.lines().filter(|l| !l.trim().is_empty());
        let header = csv_split(lines.next().ok_or("an empty file")?);
        lines.map(|l| header.iter().cloned().zip(csv_split(l)).collect()).collect()
    };
    let mut e = h.field.edits().clone();
    let (mut changed, mut added) = (0, 0);
    for (i, row) in rows.iter().enumerate() {
        let get = |c: &str| row.get(c).map(|v| v.trim()).filter(|v| !v.is_empty());
        let num = |c: &str| -> R<Option<f64>> { get(c).map(|v| v.parse::<f64>().map_err(|_| format!("row {}: {c} `{v}` is not a number", i + 1).into())).transpose() };
        let (ra, dec, ly) = (num("ra_deg")?, num("dec_deg")?, num("distance_ly")?);
        let name = get("name").map(str::to_string);
        match get("id") {
            Some(id) => {
                let id: StarId = id.parse()?;
                let Some(s) = h.field.star(&id) else { continue };
                let (r0, d0, dist0) = h.frame.radec_of(s.position_pc);
                let position_pc = h.frame.galactic_of(ra.unwrap_or(r0), dec.unwrap_or(d0), ly.map_or(dist0, |l| l / LY_PER_PC));
                let k = kind(
                    num("initial_mass")?.unwrap_or(s.kind.initial_mass),
                    num("age_gyr")?.unwrap_or(s.kind.age_gyr),
                    num("metallicity_dex")?.unwrap_or(s.kind.metallicity_dex),
                    num("companion_mass")?.or(s.kind.companion_mass),
                )?;
                let moved = s.distance_pc(position_pc) > 1e-3 * (1.0 + dist0);
                let renamed = name.as_ref().is_some_and(|n| *n != h.field.name_of(&id));
                let remade = (k.initial_mass - s.kind.initial_mass).abs() > 1e-4
                    || (k.age_gyr - s.kind.age_gyr).abs() > 1e-4
                    || (k.metallicity_dex - s.kind.metallicity_dex).abs() > 1e-3
                    || k.companion_mass.is_some() != s.kind.companion_mass.is_some();
                if moved || renamed || remade {
                    let keep = e.record(&id).and_then(|r| r.name.clone());
                    e.put(StarRecord { id, name: if renamed { name } else { keep }, position_pc, kind: k });
                    changed += 1;
                }
            }
            None => {
                let (Some(ra), Some(dec), Some(ly), Some(m)) = (ra, dec, ly, num("initial_mass")?) else {
                    return Err(format!("row {}: a new star needs ra_deg, dec_deg, distance_ly and initial_mass", i + 1).into());
                };
                let id = e.next_added();
                e.put(StarRecord {
                    id,
                    name,
                    position_pc: h.frame.galactic_of(ra, dec, ly / LY_PER_PC),
                    kind: kind(m, num("age_gyr")?.unwrap_or(1.0), num("metallicity_dex")?.unwrap_or(0.0), num("companion_mass")?)?,
                });
                added += 1;
            }
        }
    }
    h.save(&e)?;
    println!("imported {}: {changed} stars changed, {added} added", from.display());
    Ok(())
}

fn csv_field(v: &str) -> String {
    if v.contains([',', '"', '\n']) {
        format!("\"{}\"", v.replace('"', "\"\""))
    } else {
        v.to_string()
    }
}

/// One CSV line's fields, with quotes.
fn csv_split(line: &str) -> Vec<String> {
    let (mut out, mut cur, mut quoted) = (Vec::new(), String::new(), false);
    let mut chars = line.chars().peekable();
    while let Some(c) = chars.next() {
        match (c, quoted) {
            ('"', true) if chars.peek() == Some(&'"') => {
                cur.push('"');
                chars.next();
            }
            ('"', _) => quoted = !quoted,
            (',', false) => out.push(std::mem::take(&mut cur)),
            _ => cur.push(c),
        }
    }
    out.push(cur);
    out
}

/// Add a constellation from lines `A-B` between stars (identities or names).
pub fn constellation_add(dir: &Path, name: &str, lines: &[String]) -> R<()> {
    let h = Here::open(dir)?;
    let mut figures = h.world.load_constellations()?;
    let mut c = Constellation { name: name.into(), lines: Vec::new() };
    for l in lines {
        let (a, b) = l.split_once('-').filter(|(a, b)| !a.is_empty() && !b.is_empty()).ok_or_else(|| format!("a line is `A-B`, not `{l}`"))?;
        c.lines.push([h.find(a)?, h.find(b)?]);
    }
    figures.retain(|f| !f.name.eq_ignore_ascii_case(name));
    figures.push(c);
    h.world.save_constellations(&figures)?;
    println!("constellation {name}: {} lines", lines.len());
    Ok(())
}

pub fn constellation_remove(dir: &Path, name: &str) -> R<()> {
    let w = World::open(dir)?;
    let mut figures = w.load_constellations()?;
    let before = figures.len();
    figures.retain(|f| !f.name.eq_ignore_ascii_case(name));
    if figures.len() == before {
        return Err(format!("no constellation `{name}`").into());
    }
    w.save_constellations(&figures)?;
    println!("removed {name}");
    Ok(())
}

pub fn constellation_list(dir: &Path) -> R<()> {
    let h = Here::open(dir)?;
    for f in h.world.load_constellations()? {
        println!("{}:", f.name);
        for [a, b] in &f.lines {
            println!("  {} — {}", h.describe(a), h.describe(b));
        }
    }
    Ok(())
}

/// Save the automatic figures as the world's constellations.
pub fn constellation_auto(dir: &Path, replace: bool) -> R<()> {
    let w = World::open(dir)?;
    if !replace && !w.load_constellations()?.is_empty() {
        return Err("the world has constellations of its own (--replace to replace them)".into());
    }
    let system = w.load_system()?.ok_or("this world has no solar system (a hand-made body)")?;
    let field = starfield(&w)?;
    let s = geoforge_almanac::sky::sky(&system, w.manifest.focus, 0.0, 0.0, &field)?;
    let figures = s.automatic_constellations();
    for f in &figures {
        println!("{}: {} lines", f.name, f.lines.len());
    }
    w.save_constellations(&figures)?;
    println!("saved {} constellations", figures.len());
    Ok(())
}

/// The galaxies within `radius_mpc` of the world's own, brighter than `limit` from it.
pub fn universe(dir: &Path, radius_mpc: f64, limit: f64) -> R<()> {
    use geoforge_cosmos::universe::Universe;
    let w = World::open(dir)?;
    let system = w.load_system()?.ok_or("this world has no solar system (a hand-made body)")?;
    let site = system.root.site.as_ref().ok_or("this system has no place in a galaxy")?;
    let field = starfield(&w)?;
    let Some(u) = Universe::of_site(site, system.root.seed) else {
        println!("This world is in the real universe: its sky has the Local Group and the bright galaxies near it (`geoforge sky`).");
        return Ok(());
    };
    let home = u.home_at(system.root.placement.as_ref(), field.galaxy());
    println!(
        "Home: {} ({}, {}), {:.1e} M☉ of stars in a {:.1e} M☉ halo, at ({:.1}, {:.1}, {:.1}) Mpc in universe {}",
        field.galaxy().name,
        home.kind(),
        home.id,
        home.stellar_msun,
        home.halo_msun,
        home.position_mpc[0],
        home.position_mpc[1],
        home.position_mpc[2],
        u.seed().0
    );
    let mut near: Vec<_> = u
        .galaxies_seen(home.position_mpc, radius_mpc, limit)
        .into_iter()
        .filter(|g| g.id != home.id)
        .map(|g| {
            let d = (0..3).map(|k| (g.position_mpc[k] - home.position_mpc[k]).powi(2)).sum::<f64>().sqrt();
            let v = g.absolute_v() + 5.0 * (d * 1e5).max(1e-9).log10();
            (d, v, g)
        })
        .filter(|(_, v, _)| *v <= limit)
        .collect();
    near.sort_by(|a, b| a.0.total_cmp(&b.0));
    let with_home = |g: &geoforge_cosmos::universe::UniverseGalaxy| {
        if g.id.bin == home.id.bin && g.id.cube == home.id.cube && g.id.n == home.id.n {
            "satellite"
        } else if g.is_central() {
            ""
        } else {
            "satellite of another"
        }
    };
    println!("{} galaxies within {radius_mpc} Mpc brighter than V {limit}:", near.len());
    println!("  {:<18} {:<20} {:<8} {:>8} {:>10} {:>6}", "name", "identity", "type", "Mpc", "stars M☉", "V");
    for (d, v, g) in &near {
        println!("  {:<18} {:<20} {:<8} {:>8.3} {:>10.1e} {:>6.1} {}", g.name(), g.id.to_string(), g.kind(), d, g.stellar_msun, v, with_home(g));
    }
    Ok(())
}

/// The system at `star` of `dir`'s galaxy: in `dir`'s universe, its focus body made there;
/// for a world on its own, a new world in `to`, the galaxy and its edits carried over.
pub fn visit(dir: &Path, star: &str, to: Option<&Path>) -> R<()> {
    let h = Here::open(dir)?;
    let id = h.find(star)?;
    let generated = geoforge_cosmos::galactic::system_at(&h.field, &id).ok_or("that star is gone")?.map_err(|b| format!("{b}: no world to stand on there"))?;
    // A system placed on its own is no star of the galaxy: make its star one, in both worlds,
    // so the new world sees it.
    let mut edits = h.field.edits().clone();
    let system = h.world.load_system()?.ok_or("no system")?;
    let here = h.frame.observer.position_pc;
    let home_star = system.root.site.as_ref().and_then(|s| s.star);
    if home_star.is_none() && !edits.add.iter().any(|r| r.position_pc == here) {
        let stars: Vec<(&str, &geoforge_core::body::Star)> =
            system.nodes.iter().filter_map(|n| if let geoforge_core::body::NodeKind::Star(s) = &n.node { Some((n.name.as_str(), s)) } else { None }).collect();
        if let Some((name, primary)) = stars.first() {
            let id = edits.next_added();
            let kind = StarKind { initial_mass: primary.initial_mass(), age_gyr: primary.age_gyr, metallicity_dex: primary.metallicity_dex, companion_mass: stars.get(1).map(|s| s.1.initial_mass()) };
            edits.put(StarRecord { id, name: Some((*name).to_string()), position_pc: here, kind });
            h.save(&edits)?;
        }
    }
    let grid = h.world.manifest.grid;
    if let (Some(root), Some(here)) = (&h.world.universe, h.world.manifest.address) {
        let mut u = geoforge_store::UniverseDir::open(root)?;
        let w = geoforge_pipeline::create_at(&mut u, &here.up_to(geoforge_core::Level::Galaxy).star(geoforge_core::StarRef::Star(id)), None, grid)?;
        crate::system::show(&w)?;
        println!("\n{} at {} — now the focus of {}", w.manifest.focus_name, h.describe(&id), root.display());
        return Ok(());
    }
    let to = to.ok_or("this world is on its own: give the new world's directory with --to")?;
    let body = geoforge_cosmos::any_body(&generated.system, generated.focus).ok_or("the system's focus has no body")?;
    let seed = id.seed(h.field.galaxy_seed());
    let w =
        geoforge_pipeline::create_world(to, &geoforge_pipeline::NewWorld { seed, recipe: geoforge_cosmos::Recipe::default(), grid, body: Some(body), system: Some(generated.system), focus: None })?;
    if let Some(g) = h.world.load_galaxy_model::<geoforge_cosmos::Galaxy>()? {
        w.save_galaxy_model(&g)?;
    }
    w.save_star_edits(&edits)?;
    crate::system::show(&w)?;
    println!("\ncreated {} at {} — next: geoforge generate {}", to.display(), h.describe(&id), to.display());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn csv_quotes() {
        assert_eq!(csv_split(r#"a,"b, c","d ""e""""#), vec!["a", "b, c", r#"d "e""#]);
        assert_eq!(csv_split(&["x", &csv_field("y, \"z\"")].join(",")), vec!["x", "y, \"z\""]);
    }
}
