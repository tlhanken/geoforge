//! Worlds set against an observed one, layer by layer: the A/B page.
//!
//! The reference world holds observed data as authored products: Earth's relief, winds, sea
//! temperatures, climate. Each model world ran some layers on real inputs. For every product
//! the reference holds, the page shows the observed map and then each model's, with the same
//! renderer and key, and gives the numbers that say how close the model is. A model's copy of
//! an observed input, which is the same edit file, is not compared: it is an input, not a
//! result. The observed data is resampled onto each model's grid, so the grids may differ.

use crate::report::{caption, esc, fit, key_html, lat_lon, map, num, CSS, GROUPS, SCRIPT};
use crate::{views, PipelineError, Result, View};
use geoforge_core::grid::{AnyGrid, Cell, Grid};
use geoforge_core::product::{ProductData, ProductId, Products, Source};
use geoforge_core::BodyContext;
use geoforge_store::{codec, World};
use image::RgbImage;
use serde::Serialize;
use std::fmt::Write;

/// A model world and its label on the page.
pub struct Model<'a> {
    pub label: String,
    pub world: &'a World,
}

/// A named point or latitude–longitude box whose values are listed, observed beside each
/// model's (a station, a mountain range, a boundary current).
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Probe {
    pub name: String,
    pub lat: (f64, f64),
    pub lon: (f64, f64),
}

/// What goes into a comparison page.
#[derive(Debug, Clone, PartialEq)]
pub struct Options {
    /// Width of the maps, px.
    pub width: u32,
    pub probes: Vec<Probe>,
    /// Where each observed product came from, by product name ("temperature" → "ERA-Interim
    /// 2 m, 2005"), for the observed map's caption.
    pub sources: Vec<(String, String)>,
}

impl Default for Options {
    fn default() -> Self {
        Self { width: 1200, probes: Vec::new(), sources: Vec::new() }
    }
}

/// How close a model's map is to the observed one.
#[derive(Debug, Clone, PartialEq, Serialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum Metrics {
    /// A quantity: area-weighted model − observed (`bias`, `rmse`) and pattern correlation
    /// (of the logarithm for rain), over `over`; `extra` lists (measure, observed, model).
    Scalar { over: &'static str, bias: f64, rmse: f64, corr: f64, log: bool, extra: Vec<(String, f64, f64)> },
    /// A flow: area share where the direction is within 90° of the observed, over the sea
    /// and over land below 1.5 km, and the model's mean speed over the observed.
    Vector { agree_sea: f64, agree_low_land: f64, speed_ratio: f64 },
    /// Classes: area share of land where they agree (for Köppen, also by group), and each
    /// class's share of the land (observed, model).
    Classes { agree: f64, agree_group: Option<f64>, shares: Vec<(String, f64, f64)> },
}

/// One product or view: the metrics of each model, and the probes' values (observed first).
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Section {
    pub name: String,
    pub title: String,
    pub models: Vec<(String, Metrics)>,
    pub probes: Vec<(String, Vec<Option<f64>>)>,
}

/// A comparison page, its images (`img/<name>.jpg`) and its numbers.
#[derive(Debug, Clone)]
pub struct Comparison {
    pub title: String,
    pub html: String,
    pub images: Vec<(String, RgbImage)>,
    pub sections: Vec<Section>,
}

/// Ground below this is "low land" for the winds: the observed winds are at 850 hPa, which
/// lies under higher ground.
const LOW_LAND_M: f64 = 1500.0;

/// A product's values, per cell.
enum Values {
    Scalar(Vec<f64>),
    Vector(Vec<[f64; 2]>),
    Classes(Vec<Option<String>>),
}

impl Values {
    fn of(p: ProductId, products: &Products) -> Option<Self> {
        let scalar = |v: &[f32]| Some(Values::Scalar(v.iter().map(|&x| f64::from(x)).collect()));
        let vector = |v: &[[f32; 2]]| Some(Values::Vector(v.iter().map(|w| [f64::from(w[0]), f64::from(w[1])]).collect()));
        match products.get(p)? {
            ProductData::CrustAge(f)
            | ProductData::RockAge(f)
            | ProductData::RoughElevation(f)
            | ProductData::SstAnomaly(f)
            | ProductData::Temperature(f)
            | ProductData::TemperatureRange(f)
            | ProductData::Precipitation(f)
            | ProductData::PrecipitationWarmShare(f)
            | ProductData::Elevation(f)
            | ProductData::ZonalWind(f)
            | ProductData::CloudBands(f) => scalar(&f.values),
            ProductData::WindJan(f) | ProductData::WindJul(f) | ProductData::OceanCurrents(f) => vector(&f.values),
            ProductData::Biome(b) => Some(Values::Classes(b.owner.values.iter().map(|&k| b.table.get(geoforge_core::Key(u64::from(k))).map(|c| c.slug.clone())).collect())),
            _ => None,
        }
    }

    /// The values at each of `cells` (indices into these values).
    fn pick(&self, cells: &[usize]) -> Self {
        match self {
            Values::Scalar(v) => Values::Scalar(cells.iter().map(|&c| v[c]).collect()),
            Values::Vector(v) => Values::Vector(cells.iter().map(|&c| v[c]).collect()),
            Values::Classes(v) => Values::Classes(cells.iter().map(|&c| v[c].clone()).collect()),
        }
    }
}

/// Where a product's numbers are read.
#[derive(Clone, Copy, PartialEq)]
enum Over {
    Sea,
    Land,
    All,
}

fn over(p: ProductId) -> Over {
    match p {
        ProductId::SstAnomaly | ProductId::OceanCurrents => Over::Sea,
        ProductId::Temperature | ProductId::TemperatureRange | ProductId::Precipitation | ProductId::PrecipitationWarmShare | ProductId::Biome => Over::Land,
        _ => Over::All,
    }
}

/// One world's maps and the observed values on its grid.
struct Side<'a> {
    world: &'a World,
    products: Products,
    grid: AnyGrid,
    body: BodyContext,
}

impl<'a> Side<'a> {
    fn open(world: &'a World) -> Result<Self> {
        Ok(Self { world, products: world.load_products()?, grid: world.manifest.grid.build()?, body: world.load_body()? })
    }

    fn render(&self, name: &str) -> Result<Option<(RgbImage, Option<geoforge_store::render::Legend>)>> {
        let source = views::Source::new(&self.products, &self.grid, &self.body);
        map(name, &self.products, &self.grid, &self.body, &source, &source.sea())
    }

    fn koppen(&self) -> Option<Vec<Option<String>>> {
        let source = views::Source::new(&self.products, &self.grid, &self.body);
        match views::compute(View::Koppen, &source)? {
            views::Map::Classes(c) => Some(c.cells.iter().map(|k| k.map(|k| c.legend[k].code.clone())).collect()),
            _ => None,
        }
    }

    /// The cell of this world nearest each of `other`'s cells.
    fn cells_under(&self, other: &AnyGrid) -> Vec<usize> {
        if self.grid.spec() == other.spec() {
            return (0..other.len()).collect();
        }
        (0..other.len() as Cell).map(|c| self.grid.locate(other.center(c)) as usize).collect()
    }

    /// Ground height per cell: the elevation after erosion, or before.
    fn height(&self) -> Option<Vec<f64>> {
        [ProductId::Elevation, ProductId::RoughElevation].into_iter().find_map(|p| match Values::of(p, &self.products) {
            Some(Values::Scalar(v)) => Some(v),
            _ => None,
        })
    }
}

/// The comparison page: each product `reference` holds as observed data, beside each model's.
pub fn build(reference: &World, models: &[Model], opts: &Options) -> Result<Comparison> {
    if models.is_empty() {
        return Err(PipelineError::Invalid("compare needs at least one model world".into()));
    }
    let real = Side::open(reference)?;
    let sides: Vec<Side> = models.iter().map(|m| Side::open(m.world)).collect::<Result<_>>()?;
    let observed: Vec<ProductId> = ProductId::ALL.into_iter().filter(|&p| matches!(reference.manifest.source(p), Source::Authored { .. }) && real.products.get(p).is_some()).collect();
    let real_height = real.height();
    let mut sections = Vec::new();
    let mut images = Vec::new();
    let mut body_html = String::new();
    let mut nav = Vec::new();

    for (anchor, heading, names) in GROUPS {
        let mut group = String::new();
        for &name in names {
            let section = if name == "koppen" {
                // Köppen read from the observed climate, when the reference holds one.
                let climate = [ProductId::Temperature, ProductId::TemperatureRange, ProductId::Precipitation].iter().all(|p| observed.contains(p));
                if !climate {
                    continue;
                }
                let Some(real_k) = real.koppen() else { continue };
                let results: Vec<(usize, Metrics)> = sides
                    .iter()
                    .enumerate()
                    .filter_map(|(i, s)| s.koppen().map(|k| (i, classes(&Values::Classes(real_k.clone()).pick(&real.cells_under(&s.grid)), &k, &s.grid, None, true))))
                    .collect();
                Some(("koppen", "Köppen–Geiger climates; observed: read from the observed temperature and rain".to_string(), results))
            } else {
                let Some(p) = ProductId::from_name(name) else { continue };
                if !observed.contains(&p) {
                    continue;
                }
                let Some(real_v) = Values::of(p, &real.products) else { continue };
                let real_hash = codec::read_content_hash(&reference.edits_dir(), p);
                let mut results = Vec::new();
                for (i, s) in sides.iter().enumerate() {
                    let input = matches!(s.world.manifest.source(p), Source::Authored { .. }) && codec::read_content_hash(&s.world.edits_dir(), p) == real_hash;
                    let Some(model_v) = Values::of(p, &s.products) else { continue };
                    if input {
                        continue;
                    }
                    let under = real.cells_under(&s.grid);
                    let height: Option<Vec<f64>> = real_height.as_ref().map(|h| under.iter().map(|&c| h[c]).collect());
                    let r = real_v.pick(&under);
                    let m = match (&r, &model_v) {
                        (Values::Scalar(r), Values::Scalar(m)) => scalar(p, r, m, &s.grid, height.as_deref()),
                        (Values::Vector(r), Values::Vector(m)) => vector(r, m, &s.grid, height.as_deref()),
                        (Values::Classes(_), Values::Classes(m)) => classes(&r, m, &s.grid, height.as_deref(), false),
                        _ => continue,
                    };
                    results.push((i, m));
                }
                Some((name, caption(name), results))
            };
            let Some((name, title, results)) = section else { continue };
            if results.is_empty() {
                continue;
            }
            // Maps: observed, then each compared model, with the observed map's key.
            let Some((img, key)) = real.render(name)? else { continue };
            let source = opts.sources.iter().find(|(p, _)| p == name).map_or(String::new(), |(_, s)| s.clone());
            let mut figs = figure(&format!("real_{name}"), "Observed", &source);
            images.push((format!("real_{name}"), fit(&img, opts.width)));
            for (i, _) in &results {
                if let Some((img, _)) = sides[*i].render(name)? {
                    let id = format!("{}_{name}", slug(&models[*i].label));
                    figs += &figure(&id, &format!("Model: {}", models[*i].label), "");
                    images.push((id, fit(&img, opts.width)));
                }
            }
            let probes = probe_rows(name, &real, &sides, &results, &opts.probes);
            let labels: Vec<&str> = results.iter().map(|(i, _)| models[*i].label.as_str()).collect();
            let _ = write!(
                group,
                r#"<h3>{}</h3><div class="pair">{figs}</div>{}<div class="cols">{}{}</div>"#,
                esc(&title),
                key.as_ref().map_or(String::new(), |k| format!(r#"<div class="card keycard">{}</div>"#, key_html(k))),
                metrics_card(&labels, &results.iter().map(|(_, m)| m).collect::<Vec<_>>()),
                probes_card(name, &labels, &probes, &opts.probes)
            );
            sections.push(Section {
                name: name.to_string(),
                title,
                models: results.into_iter().map(|(i, m)| (models[i].label.clone(), m)).collect(),
                probes: opts.probes.iter().map(|p| p.name.clone()).zip(probes).filter(|(_, v)| v.iter().any(Option::is_some)).collect(),
            });
        }
        if !group.is_empty() {
            nav.push((anchor, heading));
            let _ = write!(body_html, r#"<section class="band" id="{anchor}"><div class="wrap"><h2>{}</h2>{group}</div></section>"#, esc(heading));
        }
    }
    if sections.is_empty() {
        return Err(PipelineError::Invalid(format!("{} holds no observed (authored) product that a model world also has", real.body.name)));
    }

    let title = format!("{} A/B, layer by layer", real.body.name);
    let model_list: String = models.iter().zip(&sides).map(|(m, s)| format!("<li><b>{}</b>: {} cells</li>", esc(&m.label), s.grid.len())).collect();
    let nav_html: String = nav.iter().map(|(a, t)| format!(r##"<a href="#{a}">{}</a>"##, esc(t))).collect();
    let json = serde_json::to_string(&sections).unwrap_or_default().replace("</", "<\\/");
    let html = format!(
        r#"<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{t}</title>
<link rel="preconnect" href="https://fonts.googleapis.com"><link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;600&family=Spectral:wght@500&display=swap" rel="stylesheet">
<style>{CSS}
.pair{{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,420px),1fr));gap:14px;margin-top:12px}}.keycard{{margin-top:8px}}</style></head><body>
<header class="band cover"><div class="wrap"><h1>{t}</h1><p class="lede">Each layer ran on observed inputs and is set beside the observed result, so what differs belongs to that layer. The observed map comes first; the models follow, drawn by the same renderer with the same key. Observed data is resampled onto each model's grid for the numbers.</p><ul>{model_list}</ul></div></header>
<nav class="nav">{nav_html}</nav>
{body_html}
<script type="application/json" id="metrics">{json}</script>
<script>{SCRIPT}</script></body></html>"#,
        t = esc(&title),
    );
    Ok(Comparison { title, html, images, sections })
}

fn slug(s: &str) -> String {
    s.chars().map(|c| if c.is_ascii_alphanumeric() { c.to_ascii_lowercase() } else { '_' }).collect()
}

fn figure(id: &str, title: &str, note: &str) -> String {
    format!(
        r#"<figure><img loading="lazy" src="img/{id}.jpg" alt="{t}"><figcaption><b>{t}</b>{}</figcaption></figure>"#,
        if note.is_empty() { String::new() } else { format!(" {}", esc(note)) },
        t = esc(title)
    )
}

/// Whether a cell counts, by the observed ground height (sea below 0 m).
fn counts(o: Over, height: Option<&[f64]>, c: usize) -> bool {
    match (o, height) {
        (Over::All, _) | (_, None) => true,
        (Over::Sea, Some(h)) => h[c] < 0.0,
        (Over::Land, Some(h)) => h[c] >= 0.0,
    }
}

/// Area-weighted bias, RMSE and correlation of `m` against `r`.
fn scalar(p: ProductId, r: &[f64], m: &[f64], grid: &AnyGrid, height: Option<&[f64]>) -> Metrics {
    let o = over(p);
    let log = p == ProductId::Precipitation;
    let tf = |x: f64| if log { (x.max(0.0) + 50.0).ln() } else { x };
    let cells: Vec<usize> = (0..r.len()).filter(|&c| counts(o, height, c) && r[c].is_finite() && m[c].is_finite()).collect();
    let a: Vec<f64> = cells.iter().map(|&c| grid.area_sr(c as Cell)).collect();
    let sum_a: f64 = a.iter().sum::<f64>().max(1e-300);
    let mean = |f: &dyn Fn(usize) -> f64| cells.iter().zip(&a).map(|(&c, w)| w * f(c)).sum::<f64>() / sum_a;
    let bias = mean(&|c| m[c] - r[c]);
    let rmse = mean(&|c| (m[c] - r[c]).powi(2)).sqrt();
    let (rt, mt): (Vec<f64>, Vec<f64>) = r.iter().zip(m).map(|(&a, &b)| (tf(a), tf(b))).unzip();
    let corr = weighted_corr(&cells, &a, &rt, &mt);
    let mut extra = Vec::new();
    if matches!(p, ProductId::Elevation | ProductId::RoughElevation) {
        // Hypsometry: how much is land, and how high.
        let all: f64 = (0..r.len()).map(|c| grid.area_sr(c as Cell)).sum::<f64>().max(1e-300);
        let share = |v: &[f64], lo: f64| (0..v.len()).filter(|&c| v[c] >= lo).map(|c| grid.area_sr(c as Cell)).sum::<f64>() / all;
        let land_mean = |v: &[f64]| {
            let (s, w) = (0..v.len()).filter(|&c| v[c] >= 0.0).fold((0.0, 0.0), |(s, w), c| (s + grid.area_sr(c as Cell) * v[c], w + grid.area_sr(c as Cell)));
            s / w.max(1e-300)
        };
        extra.push(("land, share of the surface".to_string(), share(r, 0.0), share(m, 0.0)));
        extra.push(("mean height of land, m".to_string(), land_mean(r), land_mean(m)));
        extra.push(("above 1 km, share of the surface".to_string(), share(r, 1000.0), share(m, 1000.0)));
        extra.push(("above 2 km, share of the surface".to_string(), share(r, 2000.0), share(m, 2000.0)));
        let land: Vec<usize> = (0..r.len()).filter(|&c| r[c] >= 0.0 && m[c] >= 0.0).collect();
        let w: Vec<f64> = land.iter().map(|&c| grid.area_sr(c as Cell)).collect();
        extra.push(("height correlation where both are land".to_string(), f64::NAN, weighted_corr(&land, &w, r, m)));
    }
    Metrics::Scalar { over: ["sea", "land", "everywhere"][o as usize], bias, rmse, corr, log, extra }
}

/// Area-weighted Pearson correlation of `r` and `m` over `cells` (weights `w`).
fn weighted_corr(cells: &[usize], w: &[f64], r: &[f64], m: &[f64]) -> f64 {
    let sw: f64 = w.iter().sum::<f64>().max(1e-300);
    let mean = |f: &dyn Fn(usize) -> f64| cells.iter().zip(w).map(|(&c, w)| w * f(c)).sum::<f64>() / sw;
    let (mr, mm) = (mean(&|c| r[c]), mean(&|c| m[c]));
    let cov = mean(&|c| (r[c] - mr) * (m[c] - mm));
    let (vr, vm) = (mean(&|c| (r[c] - mr).powi(2)), mean(&|c| (m[c] - mm).powi(2)));
    cov / (vr * vm).sqrt().max(1e-300)
}

/// Direction agreement over the sea and low land, and the speed ratio.
fn vector(r: &[[f64; 2]], m: &[[f64; 2]], grid: &AnyGrid, height: Option<&[f64]>) -> Metrics {
    let agree = |keep: &dyn Fn(usize) -> bool| {
        let (mut yes, mut all) = (0.0, 0.0);
        for c in (0..r.len()).filter(|&c| keep(c) && r[c][0].hypot(r[c][1]) > 0.0 && m[c][0].hypot(m[c][1]) > 0.0) {
            let a = grid.area_sr(c as Cell);
            all += a;
            if r[c][0] * m[c][0] + r[c][1] * m[c][1] > 0.0 {
                yes += a;
            }
        }
        if all > 0.0 {
            yes / all
        } else {
            f64::NAN
        }
    };
    let h = |c: usize| height.map_or(0.0, |h| h[c]);
    let (mut sr, mut sm) = (0.0, 0.0);
    for c in 0..r.len() {
        let a = grid.area_sr(c as Cell);
        sr += a * r[c][0].hypot(r[c][1]);
        sm += a * m[c][0].hypot(m[c][1]);
    }
    Metrics::Vector { agree_sea: agree(&|c| h(c) < 0.0), agree_low_land: agree(&|c| (0.0..LOW_LAND_M).contains(&h(c))), speed_ratio: sm / sr.max(1e-300) }
}

/// Agreement and each class's share of land. Köppen codes are compared by group (A, BW, BS,
/// C, D, E) for the shares, and in full for the agreement.
fn classes(r: &Values, m: &[Option<String>], grid: &AnyGrid, height: Option<&[f64]>, koppen: bool) -> Metrics {
    let Values::Classes(r) = r else { unreachable!("classes of a class map") };
    let group = |k: &str| {
        if koppen {
            if k.starts_with('B') {
                k.chars().take(2).collect()
            } else {
                k.chars().take(1).collect()
            }
        } else {
            k.to_string()
        }
    };
    let land = |c: usize| counts(Over::Land, height, c) && (r[c].is_some() || m[c].is_some());
    let (mut yes, mut yes_group, mut all) = (0.0, 0.0, 0.0);
    let mut shares: Vec<(String, f64, f64)> = Vec::new();
    for c in (0..r.len()).filter(|&c| land(c)) {
        let a = grid.area_sr(c as Cell);
        all += a;
        if r[c].is_some() && r[c] == m[c] {
            yes += a;
        }
        if let (Some(x), Some(y)) = (&r[c], &m[c]) {
            if group(x) == group(y) {
                yes_group += a;
            }
        }
        for (k, observed) in [(&r[c], true), (&m[c], false)] {
            let Some(k) = k else { continue };
            let g = group(k);
            let i = shares.iter().position(|s| s.0 == g).unwrap_or_else(|| {
                shares.push((g, 0.0, 0.0));
                shares.len() - 1
            });
            if observed {
                shares[i].1 += a;
            } else {
                shares[i].2 += a;
            }
        }
    }
    let all = all.max(1e-300);
    for s in &mut shares {
        s.1 /= all;
        s.2 /= all;
    }
    shares.sort_by(|a, b| a.0.cmp(&b.0));
    Metrics::Classes { agree: yes / all, agree_group: koppen.then(|| yes_group / all), shares }
}

/// Each probe's value of a quantity: observed, then each compared model.
fn probe_rows(name: &str, real: &Side, sides: &[Side], results: &[(usize, Metrics)], probes: &[Probe]) -> Vec<Vec<Option<f64>>> {
    let Some(p) = ProductId::from_name(name) else { return vec![Vec::new(); probes.len()] };
    let value = |s: &Side, pr: &Probe| -> Option<f64> {
        let Some(Values::Scalar(v)) = Values::of(p, &s.products) else { return None };
        let height = s.height();
        // Heights are compared everywhere but listed for land.
        let o = if matches!(p, ProductId::Elevation | ProductId::RoughElevation) { Over::Land } else { over(p) };
        let at = |lat: f64, lon: f64| s.grid.locate(geoforge_core::UnitVec3::from_lat_lon(lat, lon)) as usize;
        // A place counts if its cell is in the product's domain (land or sea), a box if most of
        // it is.
        let inside: Vec<usize> = if pr.lat.0 == pr.lat.1 && pr.lon.0 == pr.lon.1 {
            vec![at(pr.lat.0, pr.lon.0)]
        } else {
            (0..s.grid.len())
                .filter(|&c| {
                    let ll = s.grid.center(c as Cell).to_lat_lon();
                    (pr.lat.0..=pr.lat.1).contains(&ll.lat_deg) && (pr.lon.0..=pr.lon.1).contains(&ll.lon_deg)
                })
                .collect()
        };
        let cells: Vec<usize> = inside.iter().copied().filter(|&c| counts(o, height.as_deref(), c)).collect();
        if 2 * cells.len() < inside.len().max(1) {
            return None;
        }
        let (sum, w) = cells.iter().filter(|&&c| v[c].is_finite()).fold((0.0, 0.0), |(s_, w), &c| (s_ + s.grid.area_sr(c as Cell) * v[c], w + s.grid.area_sr(c as Cell)));
        (w > 0.0).then(|| sum / w)
    };
    probes.iter().map(|pr| std::iter::once(value(real, pr)).chain(results.iter().map(|(i, _)| value(&sides[*i], pr))).collect()).collect()
}

fn pct(x: f64) -> String {
    if x.is_finite() {
        format!("{:.0}%", 100.0 * x)
    } else {
        "—".into()
    }
}

fn metrics_card(labels: &[&str], metrics: &[&Metrics]) -> String {
    let head: String = labels.iter().map(|l| format!("<th>{}</th>", esc(l))).collect();
    let row = |name: &str, f: &dyn Fn(&Metrics) -> String| format!("<tr><td>{}</td>{}</tr>", esc(name), metrics.iter().map(|m| format!(r#"<td class="n">{}</td>"#, f(m))).collect::<String>());
    let rows = match metrics[0] {
        Metrics::Scalar { over, log, extra, .. } => {
            let mut rows = format!(
                "{}{}{}",
                row(&format!("bias, over {over}"), &|m| sc(m).0),
                row("RMSE", &|m| sc(m).1),
                row(if *log { "pattern correlation (log)" } else { "pattern correlation" }, &|m| sc(m).2)
            );
            for (k, (name, observed, _)) in extra.iter().enumerate() {
                let share = name.contains("share");
                let corr = name.contains("correlation");
                let show = |x: f64| {
                    if share {
                        pct(x)
                    } else if corr {
                        format!("{x:.2}")
                    } else {
                        num(x)
                    }
                };
                let observed = if observed.is_finite() { format!(" (observed {})", show(*observed)) } else { String::new() };
                rows += &format!(
                    "<tr><td>{}{observed}</td>{}</tr>",
                    esc(name),
                    metrics.iter().map(|m| format!(r#"<td class="n">{}</td>"#, if let Metrics::Scalar { extra, .. } = m { show(extra[k].2) } else { String::new() })).collect::<String>()
                );
            }
            rows
        }
        Metrics::Vector { .. } => {
            let v = |m: &Metrics| if let Metrics::Vector { agree_sea, agree_low_land, speed_ratio } = m { (*agree_sea, *agree_low_land, *speed_ratio) } else { (f64::NAN, f64::NAN, f64::NAN) };
            format!(
                "{}{}{}",
                row("direction agrees, sea", &|m| pct(v(m).0)),
                row("direction agrees, land below 1.5 km", &|m| pct(v(m).1)),
                row("mean speed, model ÷ observed", &|m| format!("{:.2}", v(m).2))
            )
        }
        Metrics::Classes { shares, .. } => {
            let mut rows = row("land agreeing", &|m| if let Metrics::Classes { agree, .. } = m { pct(*agree) } else { String::new() });
            if let Metrics::Classes { agree_group: Some(_), .. } = metrics[0] {
                rows += &row("land agreeing in group (A, BW, BS, C, D, E)", &|m| if let Metrics::Classes { agree_group: Some(g), .. } = m { pct(*g) } else { String::new() });
            }
            for (name, observed, _) in shares {
                rows += &format!(
                    "<tr><td>{} (observed {})</td>{}</tr>",
                    esc(name),
                    pct(*observed),
                    metrics
                        .iter()
                        .map(|m| format!(
                            r#"<td class="n">{}</td>"#,
                            if let Metrics::Classes { shares, .. } = m { shares.iter().find(|s| &s.0 == name).map_or("0%".into(), |s| pct(s.2)) } else { String::new() }
                        ))
                        .collect::<String>()
                );
            }
            rows
        }
    };
    format!(r#"<div class="card"><h3>Against the observed map</h3><div class="scroll"><table><thead><tr><th></th>{head}</tr></thead><tbody>{rows}</tbody></table></div></div>"#)
}

fn sc(m: &Metrics) -> (String, String, String) {
    match m {
        Metrics::Scalar { bias, rmse, corr, .. } => (num(*bias), num(*rmse), format!("{corr:.2}")),
        _ => Default::default(),
    }
}

fn probes_card(name: &str, labels: &[&str], rows: &[Vec<Option<f64>>], probes: &[Probe]) -> String {
    let body: String = probes
        .iter()
        .zip(rows)
        .filter(|(_, v)| v.first().is_some_and(Option::is_some))
        .map(|(p, v)| {
            let at = if p.lat.0 == p.lat.1 { lat_lon(p.lat.0, p.lon.0) } else { format!("{}–{}", lat_lon(p.lat.0, p.lon.0), lat_lon(p.lat.1, p.lon.1)) };
            format!(r#"<tr><td>{} <span class="note">{at}</span></td>{}</tr>"#, esc(&p.name), v.iter().map(|x| format!(r#"<td class="n">{}</td>"#, x.map_or("—".into(), num))).collect::<String>())
        })
        .collect();
    if body.is_empty() {
        return String::new();
    }
    let head: String = labels.iter().map(|l| format!("<th>{}</th>", esc(l))).collect();
    format!(
        r#"<div class="card"><h3>At named places: {}</h3><div class="scroll"><table><thead><tr><th></th><th>Observed</th>{head}</tr></thead><tbody>{body}</tbody></table></div></div>"#,
        esc(&caption(name))
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use geoforge_core::grid::GridSpec;

    fn grid() -> AnyGrid {
        GridSpec::Equirect { width: 36, height: 18 }.build().unwrap()
    }

    #[test]
    fn a_map_against_itself_is_exact_and_an_offset_shows_as_bias() {
        let g = grid();
        let r: Vec<f64> = (0..g.len()).map(|c| g.center(c as Cell).to_lat_lon().lat_deg.to_radians().cos() * 30.0).collect();
        let Metrics::Scalar { bias, rmse, corr, .. } = scalar(ProductId::Temperature, &r, &r, &g, None) else { panic!() };
        assert!(bias.abs() < 1e-12 && rmse.abs() < 1e-12 && (corr - 1.0).abs() < 1e-9);
        let m: Vec<f64> = r.iter().map(|v| v - 2.0).collect();
        let Metrics::Scalar { bias, rmse, corr, .. } = scalar(ProductId::Temperature, &r, &m, &g, None) else { panic!() };
        assert!((bias + 2.0).abs() < 1e-9 && (rmse - 2.0).abs() < 1e-9 && (corr - 1.0).abs() < 1e-9);
    }

    #[test]
    fn sea_and_land_are_read_from_the_observed_heights() {
        let g = grid();
        // Sea west of the meridian, land east; the model is 5 °C off over the sea only.
        let height: Vec<f64> = (0..g.len()).map(|c| if g.center(c as Cell).to_lat_lon().lon_deg < 0.0 { -1000.0 } else { 500.0 }).collect();
        let r = vec![10.0; g.len()];
        let m: Vec<f64> = (0..g.len()).map(|c| if height[c] < 0.0 { 15.0 } else { 10.0 }).collect();
        let Metrics::Scalar { bias, over, .. } = scalar(ProductId::Temperature, &r, &m, &g, Some(&height)) else { panic!() };
        assert!(over == "land" && bias.abs() < 1e-9, "{over} {bias}");
        let Metrics::Scalar { bias, over, .. } = scalar(ProductId::SstAnomaly, &r, &m, &g, Some(&height)) else { panic!() };
        assert!(over == "sea" && (bias - 5.0).abs() < 1e-9, "{over} {bias}");
    }

    #[test]
    fn winds_agree_by_direction_and_compare_by_speed() {
        let g = grid();
        let height = vec![-100.0; g.len()];
        let r = vec![[4.0, 1.0]; g.len()];
        let Metrics::Vector { agree_sea, speed_ratio, .. } = vector(&r, &vec![[8.0, 2.0]; g.len()], &g, Some(&height)) else { panic!() };
        assert!((agree_sea - 1.0).abs() < 1e-12 && (speed_ratio - 2.0).abs() < 1e-9);
        let Metrics::Vector { agree_sea, agree_low_land, .. } = vector(&r, &vec![[-4.0, 0.0]; g.len()], &g, Some(&height)) else { panic!() };
        assert!(agree_sea.abs() < 1e-12 && agree_low_land.is_nan());
    }

    #[test]
    fn koppen_classes_agree_in_full_and_share_by_group() {
        let g = grid();
        let r: Vec<Option<String>> = (0..g.len()).map(|c| Some(if c % 2 == 0 { "BWh" } else { "Cfb" }.to_string())).collect();
        let m: Vec<Option<String>> = (0..g.len()).map(|c| Some(if c % 2 == 0 { "BWk" } else { "Cfb" }.to_string())).collect();
        let Metrics::Classes { agree, agree_group, shares } = classes(&Values::Classes(r), &m, &g, None, true) else { panic!() };
        // Half the cells (by area) agree in full, all of them by group; both maps are half BW,
        // half C.
        assert!((agree - 0.5).abs() < 0.02 && agree_group == Some(1.0), "{agree} {agree_group:?}");
        assert_eq!(shares.iter().map(|s| s.0.as_str()).collect::<Vec<_>>(), ["BW", "C"]);
        assert!(shares.iter().all(|s| (s.1 - s.2).abs() < 1e-12));
    }
}
