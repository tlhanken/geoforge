//! Known cases for worlds without plates or water: generated Mars and Moon, checked against
//! what MOLA, LOLA and the landers measured. Each statistic is the median of three seeds,
//! so the cases test the generator's habits rather than one draw.
//!
//! `cargo test --test known_planets -- --nocapture` prints the scorecard;
//! `cargo test --test known_planets -- --ignored` asserts the known gaps.

mod known;

use geoforge_core::grid::AnyGrid;
use geoforge_core::product::{LayerId, Products};
use geoforge_core::surface::{CrustType, GeologicProvince};
use geoforge_core::{Archetype, Cell, Grid, Key, Vec3};
use known::{cells_in, mean, Scorecard};
use std::collections::BTreeMap;
use std::sync::OnceLock;

const SEEDS: [u64; 3] = [1, 2, 3];

/// What one generated world measures.
#[derive(Default)]
struct Stats {
    /// Highest and lowest point, m from the mean surface.
    max_m: f64,
    min_m: f64,
    /// Share of the surface that is lowland (lava plains, maria).
    lowland_share: f64,
    /// Share of the lowland in the hemisphere centred on it.
    lowland_one_side: f64,
    /// Commonest highland elevation minus commonest lowland elevation, m.
    mode_gap_m: f64,
    /// Equivalent diameter of the largest impact basin, km.
    largest_basin_km: f64,
    /// Mars only, °C: global and equatorial annual means; the coldest month at 75–90°,
    /// the colder hemisphere.
    t_global: f64,
    t_equator: f64,
    t_polar_winter: f64,
}

/// The commonest elevation of `cells`, in 250 m bins.
fn mode(g: &AnyGrid, e: &[f64], cells: &[Cell]) -> f64 {
    let mut h: BTreeMap<i64, f64> = BTreeMap::new();
    for &c in cells {
        *h.entry((e[c as usize] / 250.0).floor() as i64).or_default() += g.area_sr(c);
    }
    h.into_iter()
        .max_by(|a, b| a.1.total_cmp(&b.1))
        .map_or(0.0, |(k, _)| (k as f64 + 0.5) * 250.0)
}

fn measure(a: Archetype, seed: u64) -> Stats {
    let g = known::grid();
    let body = a.context(Key(1), a.name());
    let climate = a == Archetype::Mars;
    let layers: &[LayerId] = if climate {
        &[
            LayerId::Tectonics,
            LayerId::Geology,
            LayerId::Elevation,
            LayerId::Circulation,
            LayerId::Climate,
            LayerId::Erosion,
        ]
    } else {
        &[LayerId::Tectonics, LayerId::Geology, LayerId::Elevation, LayerId::Erosion]
    };
    let p = known::run_seeded(seed, &body, &g, Products::new(), layers);
    let e = known::scalar(&p, Products::elevation);
    let all: Vec<Cell> = g.cells().collect();
    let datum = mean(&g, &all, |c| e[c as usize]);
    let crust = p.crust_type().unwrap();
    let (low, high): (Vec<Cell>, Vec<Cell>) = all
        .iter()
        .partition(|&&c| *crust.get(c) == CrustType::Oceanic);
    let area = |cells: &[Cell]| cells.iter().map(|&c| g.area_sr(c)).sum::<f64>();
    let centroid = low.iter().fold(Vec3::new(0.0, 0.0, 0.0), |s, &c| {
        s + g.center(c).v() * g.area_sr(c)
    });
    let axis = centroid.normalized().unwrap();
    let one_side: Vec<Cell> = low
        .iter()
        .copied()
        .filter(|&c| g.center(c).dot(axis) > 0.0)
        .collect();

    let prov = p.provinces().unwrap();
    let mut basins: BTreeMap<Key, f64> = BTreeMap::new();
    for &c in &all {
        if let Some(k) = prov.owner_of(c).filter(|&k| {
            prov.table
                .get(k)
                .is_some_and(|r| r.kind == GeologicProvince::ImpactBasin)
        }) {
            *basins.entry(k).or_default() += g.area_sr(c);
        }
    }
    let r = body.exterior.radius_km;
    let largest = basins.values().copied().fold(0.0, f64::max);

    let mut s = Stats {
        max_m: e.iter().copied().fold(f64::MIN, f64::max) - datum,
        min_m: e.iter().copied().fold(f64::MAX, f64::min) - datum,
        lowland_share: area(&low) / area(&all),
        lowland_one_side: area(&one_side) / area(&low),
        mode_gap_m: mode(&g, &e, &high) - mode(&g, &e, &low),
        largest_basin_km: 2.0 * (largest * r * r / std::f64::consts::PI).sqrt(),
        ..Stats::default()
    };
    if climate {
        let t = known::scalar(&p, Products::temperature);
        let tr = known::scalar(&p, Products::temperature_range);
        let band = |lo: f64, hi: f64| cells_in(&g, lo..=hi, -180.0, 180.0);
        s.t_global = mean(&g, &all, |c| t[c as usize]);
        s.t_equator = mean(&g, &band(-10.0, 10.0), |c| t[c as usize]);
        let winter = |cells: Vec<Cell>| mean(&g, &cells, |c| t[c as usize] - tr[c as usize] / 2.0);
        s.t_polar_winter = winter(band(75.0, 90.0)).min(winter(band(-90.0, -75.0)));
    }
    s
}

/// The median over the seeds of one statistic.
fn median(runs: &[Stats], f: impl Fn(&Stats) -> f64) -> f64 {
    let mut v: Vec<f64> = runs.iter().map(f).collect();
    v.sort_by(f64::total_cmp);
    v[v.len() / 2]
}

fn scorecard() -> &'static Scorecard {
    static S: OnceLock<Scorecard> = OnceLock::new();
    S.get_or_init(|| {
        let runs = |a: Archetype| -> Vec<Stats> {
            std::thread::scope(|sc| {
                let hs: Vec<_> = SEEDS.iter().map(|&seed| sc.spawn(move || measure(a, seed))).collect();
                hs.into_iter().map(|h| h.join().unwrap()).collect()
            })
        };
        let (mars, moon) = (runs(Archetype::Mars), runs(Archetype::Moon));
        let km = |v: f64| v / 1000.0;
        let mut s = Scorecard::default();

        // Mars: MOLA (Smith et al. 1999, 2001; Aharonson et al. 2001).
        s.range("mars", "highest point (Olympus Mons, 21.2 km)", km(median(&mars, |x| x.max_m)), 17.0..=26.0, "km");
        s.range("mars", "lowest point (Hellas floor, −8.2 km)", km(median(&mars, |x| x.min_m)), -10.5..=-6.0, "km");
        s.range("mars", "full range (~30 km)", km(median(&mars, |x| x.max_m - x.min_m)), 25.0..=35.0, "km")
            .gap("the extremes overshoot: 29, 35 and 39 km over seeds 1–3");
        s.range("mars", "lowland share of the surface (about a third)", median(&mars, |x| x.lowland_share), 0.28..=0.42, "");
        s.range("mars", "lowland in one hemisphere", median(&mars, |x| x.lowland_one_side), 0.85..=1.0, "");
        s.range("mars", "highland peak above the lowland peak (5.5 km)", km(median(&mars, |x| x.mode_gap_m)), 4.0..=7.0, "km");
        s.range("mars", "largest impact basin (Hellas, ~2300 km)", median(&mars, |x| x.largest_basin_km), 1500.0..=3200.0, "km");
        // Mars climate: −63 °C global mean; InSight at 4.5°N measured 217.5 K (−55.6 °C) below
        // the daily swing; CO₂ frost holds the winter caps near 148 K (−125 °C).
        s.range("mars", "global mean temperature (−63)", median(&mars, |x| x.t_global), -71.0..=-55.0, "°C");
        s.range("mars", "equatorial annual mean (InSight −55.6)", median(&mars, |x| x.t_equator), -66.0..=-45.0, "°C");
        s.range("mars", "polar winter (CO₂ frost point −125)", median(&mars, |x| x.t_polar_winter), -135.0..=-115.0, "°C")
            .gap("polar winters fall below the CO₂ frost point: condensing CO₂ should release latent heat and hold the cap near 148 K");

        // Moon: LOLA (Smith et al. 2010); maria 16% of the surface, 31% of the near side and
        // ~1–2% of the far side; South Pole–Aitken, 2500 km across.
        s.range("moon", "highest point (SPA's north rim, 10.8 km)", km(median(&moon, |x| x.max_m)), 7.5..=13.0, "km")
            .gap("no far-side highland bulge: the highest point is a crater rim, not a regional rise");
        s.range("moon", "lowest point (Antoniadi in SPA, −9.1 km)", km(median(&moon, |x| x.min_m)), -11.0..=-7.0, "km");
        s.range("moon", "maria share of the surface (16%)", median(&moon, |x| x.lowland_share), 0.13..=0.20, "");
        s.range("moon", "maria on one side (31% near vs ~2% far)", median(&moon, |x| x.lowland_one_side), 0.85..=1.0, "");
        s.range("moon", "largest impact basin (South Pole–Aitken, 2500 km)", median(&moon, |x| x.largest_basin_km), 1800.0..=3200.0, "km")
            .gap("nothing the size of South Pole–Aitken: the largest basin is ~1000 km across");
        s
    })
}

#[test]
fn planets_known_cases() {
    let s = scorecard();
    println!("{}", s.report());
    s.assert_cases("planets");
}

#[test]
#[ignore = "known gaps: run with --ignored"]
fn planets_known_gaps() {
    scorecard().assert_gaps();
}
