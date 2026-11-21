//! Known cases, circulation → climate → hydrology → biomes: Earth's real topography (ETOPO1,
//! one degree) goes in as the elevation; winds, currents, sea-surface anomalies,
//! temperature, rain, rivers and biomes come out and are checked against what Earth does.
//! The biome classifier is also fed station climates directly, so its own errors show
//! apart from the climate model's.
//!
//! `cargo test --test known_earth_climate -- --nocapture` prints the scorecard;
//! `cargo test --test known_earth_climate -- --ignored` asserts the known gaps.

mod known;

use geoforge_core::grid::AnyGrid;
use geoforge_core::product::{LayerId, Products};
use geoforge_core::{Archetype, Cell, Field, Grid, Key};
use geoforge_surface::biomes::{classify, BiomesParams, CellClimate};
use known::{cells_in, mean, near, Scorecard};
use std::collections::BTreeMap;
use std::ops::RangeInclusive;
use std::sync::OnceLock;

/// Climate normals of a station (approximate, rounded).
struct Station {
    name: &'static str,
    lat: f64,
    lon: f64,
    elev_m: f64,
    /// Annual mean, °C.
    t: f64,
    /// Warmest minus coldest month, °C.
    range: f64,
    /// mm/yr.
    precip: f64,
    /// Share of the year's rain in the warm half-year.
    warm_share: f64,
    koppen: &'static str,
    /// Biomes the station's surroundings belong to (WWF terrestrial biomes, in this model's
    /// classes).
    biomes: &'static [&'static str],
}

#[allow(clippy::too_many_arguments)]
const fn st(
    name: &'static str,
    lat: f64,
    lon: f64,
    elev_m: f64,
    t: f64,
    range: f64,
    precip: f64,
    warm_share: f64,
    koppen: &'static str,
    biomes: &'static [&'static str],
) -> Station {
    Station {
        name,
        lat,
        lon,
        elev_m,
        t,
        range,
        precip,
        warm_share,
        koppen,
        biomes,
    }
}

const TF: &str = "temperate_deciduous_forest";
#[rustfmt::skip]
const STATIONS: &[Station] = &[
    st("London", 51.5, -0.1, 25.0, 11.3, 13.5, 620.0, 0.48, "Cfb", &[TF]),
    st("Moscow", 55.75, 37.6, 156.0, 5.8, 25.7, 700.0, 0.62, "Dfb", &[TF, "boreal_forest"]),
    st("Yakutsk", 62.0, 129.7, 98.0, -8.8, 58.0, 240.0, 0.72, "Dfd", &["boreal_forest"]),
    st("Chicago", 41.9, -87.6, 180.0, 10.0, 27.9, 950.0, 0.6, "Dfa", &[TF, "temperate_grassland"]),
    st("Beijing", 39.9, 116.4, 44.0, 12.9, 30.0, 570.0, 0.85, "Dwa", &[TF, "temperate_grassland"]),
    st("Vancouver", 49.25, -123.1, 70.0, 10.5, 14.5, 1450.0, 0.28, "Cfb", &["temperate_rainforest"]),
    st("Reykjavik", 64.1, -21.9, 50.0, 4.7, 11.0, 850.0, 0.42, "Cfc", &["tundra", "boreal_forest"]),
    st("Los Angeles", 34.05, -118.25, 70.0, 18.6, 9.0, 380.0, 0.08, "Csb", &["mediterranean_woodland"]),
    st("Perth", -31.95, 115.85, 20.0, 18.7, 11.0, 730.0, 0.12, "Csa", &["mediterranean_woodland"]),
    st("Cape Town", -33.9, 18.45, 40.0, 17.0, 8.5, 515.0, 0.2, "Csb", &["mediterranean_woodland"]),
    st("Buenos Aires", -34.6, -58.4, 25.0, 17.9, 13.5, 1200.0, 0.58, "Cfa", &["temperate_grassland", TF]),
    st("Cairo", 30.05, 31.25, 75.0, 22.0, 14.5, 25.0, 0.15, "BWh", &["hot_desert"]),
    st("Riyadh", 24.7, 46.7, 610.0, 26.5, 22.0, 100.0, 0.15, "BWh", &["hot_desert"]),
    st("In Salah", 27.2, 2.5, 270.0, 25.5, 25.0, 15.0, 0.3, "BWh", &["hot_desert"]),
    st("Timbuktu", 16.8, -3.0, 260.0, 29.0, 12.0, 200.0, 0.95, "BWh", &["hot_desert", "savanna"]),
    st("Phoenix", 33.45, -112.07, 340.0, 24.0, 21.0, 200.0, 0.45, "BWh", &["hot_desert"]),
    st("Alice Springs", -23.7, 133.9, 545.0, 21.0, 16.0, 285.0, 0.68, "BWh", &["hot_desert"]),
    st("Lima", -12.05, -77.05, 150.0, 19.4, 7.0, 15.0, 0.6, "BWh", &["hot_desert"]),
    st("Tehran", 35.7, 51.4, 1190.0, 18.5, 27.0, 230.0, 0.15, "BSk", &["temperate_grassland", "cold_desert", "hot_desert"]),
    st("Dodge City", 37.75, -100.0, 790.0, 13.0, 26.0, 560.0, 0.75, "BSk", &["temperate_grassland"]),
    st("Ulaanbaatar", 47.9, 106.9, 1350.0, -0.4, 40.0, 270.0, 0.88, "Dwc", &["temperate_grassland", "cold_desert"]),
    st("Astana", 51.15, 71.45, 350.0, 3.5, 37.0, 330.0, 0.6, "Dfb", &["temperate_grassland"]),
    st("Manaus", -3.1, -60.0, 70.0, 27.0, 1.8, 2300.0, 0.4, "Af", &["tropical_rainforest"]),
    st("Singapore", 1.35, 103.8, 15.0, 27.5, 1.8, 2500.0, 0.45, "Af", &["tropical_rainforest"]),
    st("Kinshasa", -4.3, 15.3, 310.0, 25.3, 4.0, 1400.0, 0.6, "Aw", &["tropical_seasonal_forest", "savanna"]),
    st("Mumbai", 19.1, 72.9, 10.0, 27.2, 6.0, 2200.0, 0.97, "Am", &["tropical_seasonal_forest"]),
    st("Nairobi", -1.3, 36.8, 1660.0, 17.8, 4.0, 870.0, 0.45, "Cwb", &["savanna", "tropical_seasonal_forest"]),
    // Lhasa's valley is steppe, below the treeline; the plateau above is alpine.
    st("Lhasa", 29.65, 91.1, 3650.0, 8.0, 17.0, 430.0, 0.92, "Cwb", &["temperate_grassland", "alpine"]),
    st("Utqiagvik", 71.3, -156.8, 10.0, -11.0, 31.0, 140.0, 0.65, "ET", &["tundra"]),
    st("Summit, Greenland", 72.6, -38.5, 3200.0, -30.0, 30.0, 250.0, 0.5, "EF", &["ice_sheet"]),
    st("Vostok", -78.45, 106.8, 3490.0, -55.0, 35.0, 20.0, 0.5, "EF", &["ice_sheet"]),
];

// Known gaps, by station or place: what the model gets wrong there.
const EUROPE: &str = "maritime air fades too fast across Europe: Moscow's range is 31 °C (observed 26) and its mean 5 °C cold";
const SIBERIA: &str = "Siberian winters too mild: Yakutsk's range is 45 °C (observed 58)";
const SOUTHWEST: &str = "the interior Southwest is drenched by winter storms (Great Basin 1600 mm at 40°N 117°W, observed ~250)";
const WEST_AFRICA: &str = "the West African monsoon reaches too far north (2500 mm at 13°N 0°E, observed ~600)";
const ARABIA: &str = "the summer rain belt reaches into Arabia (600 mm at 22°N 47°E, all of it in summer), which cools it";
const LIMA: &str = "Lima runs 7 °C warm: the Humboldt's cold (now −3.9 °C offshore) does not reach the air ashore, which has no stratus deck";
const ICELAND: &str = "cold East Greenland water chills the air reaching Iceland: Reykjavik is 7 °C cold";
const KAZAKHSTAN: &str = "the Kazakh steppe is three times too wet and less continental than observed (990 mm, range 22 °C at 50°N 70°E; observed ~330 mm, 37 °C)";
const TIBET: &str = "Tibetan summers too warm (warmest month 16 °C at 4900 m, observed ~7) and the monsoon does not climb the plateau";
const ANTARCTICA: &str = "the Antarctic plateau is ~17 °C too warm";
const TEMPERATURE_GAPS: &[(&str, &str)] = &[
    ("Moscow", EUROPE),
    ("Yakutsk", SIBERIA),
    ("Reykjavik", ICELAND),
    ("Beijing", "north China runs 4 °C warm"),
    ("Tehran", "the Iranian plateau runs 6 °C cool"),
    ("Riyadh", ARABIA),
    ("Timbuktu", WEST_AFRICA),
    ("Phoenix", SOUTHWEST),
    ("Lima", LIMA),
    ("Utqiagvik", "the Arctic coast runs 4 °C warm"),
    ("Vostok", ANTARCTICA),
];
const RANGE_GAPS: &[(&str, &str)] = &[];
const RAIN_GAPS: &[(&str, &str)] = &[
    ("London", "no rain shadow in the lee of Britain's hills: London gets west-coast rain"),
    ("Yakutsk", "the Siberian interior gets an eighth of its rain"),
    ("Phoenix", SOUTHWEST),
    ("Timbuktu", WEST_AFRICA),
    ("Tehran", "the Alborz's rain spreads to its dry southern foot (485–700 mm in the cells around Tehran)"),
    ("Astana", KAZAKHSTAN),
    ("Nairobi", "the East African highlands are drenched (2800 mm, observed 870)"),
    ("Lhasa", TIBET),
];
const PLACE_GAPS: &[(&str, &str)] = &[
    ("Siberian taiga", "central Siberia is too dry for taiga (230 mm at 60°N 100°E, observed ~400–500)"),
    ("Amazon", "every one-degree lowland cell on a big river turns wetland: the wetland rule keys on discharge, not on how much of the cell floods"),
    ("Congo", "the Congo basin gets ~1200 mm (observed 1600–2000), too dry for forest"),
    ("Sahel", WEST_AFRICA),
    ("Eastern US", "the eastern US is too dry for forest (640 mm at 38°N 82°W, observed ~1100)"),
    ("Kazakh steppe", KAZAKHSTAN),
    ("Greece", "Greece is drenched by winter storms (1450 mm at 38°N 22°E, observed 400–1100), too wet for Mediterranean woodland"),
    ("Tibet", TIBET),
];
const CLASSIFIER_GAPS: &[(&str, &str)] = &[
    ("Vancouver", "temperate rainforest needs 1500 mm and a range under 14 °C; Vancouver has 1450 mm and 14.5 °C"),
    ("Reykjavik", "mild-winter subpolar climates fall through to temperate forest: boreal needs winters below −3 °C"),
    ("Tehran", "a semi-arid 18 °C mean reads as savanna even with frosty winters"),
    ("Nairobi", "tropical highlands read as temperate: the tropical classes need an 18 °C coldest month"),
];

/// Köppen major group (B split into BW and BS), Kottek et al. (2006) thresholds.
fn koppen_group(t: f64, range: f64, precip: f64, warm_share: f64) -> &'static str {
    let (warm, cold) = (t + range / 2.0, t - range / 2.0);
    let pth = 10.0
        * (2.0 * t
            + if warm_share >= 2.0 / 3.0 {
                28.0
            } else if warm_share <= 1.0 / 3.0 {
                0.0
            } else {
                14.0
            });
    if warm < 10.0 {
        "E"
    } else if precip < pth / 2.0 {
        "BW"
    } else if precip < pth {
        "BS"
    } else if cold >= 18.0 {
        "A"
    } else if cold > -3.0 {
        "C"
    } else {
        "D"
    }
}

fn observed_group(k: &str) -> &str {
    if k.starts_with('B') {
        &k[..2]
    } else {
        &k[..1]
    }
}

/// Earth's climate on its real topography, with the fields the cases read.
struct Earth {
    grid: AnyGrid,
    p: Products,
    elev: Vec<f64>,
    t: Vec<f64>,
    t_range: Vec<f64>,
    precip: Vec<f64>,
    warm_share: Vec<f64>,
}

fn earth() -> &'static Earth {
    static E: OnceLock<Earth> = OnceLock::new();
    E.get_or_init(|| {
        let g = known::grid();
        let e = known::etopo(&g);
        let mut inputs = Products::new();
        inputs.insert(known::elevation_product(&g, e.clone()));
        let body = Archetype::Earth.context(Key(1), "Earth");
        let p = known::run(
            &body,
            &g,
            inputs,
            &[
                LayerId::Circulation,
                LayerId::Climate,
                LayerId::Hydrology,
                LayerId::Biomes,
            ],
        );
        Earth {
            t: known::scalar(&p, Products::temperature),
            t_range: known::scalar(&p, Products::temperature_range),
            precip: known::scalar(&p, Products::precipitation),
            warm_share: known::scalar(&p, Products::precipitation_warm_share),
            elev: e.iter().map(|&v| f64::from(v)).collect(),
            grid: g,
            p,
        }
    })
}

impl Earth {
    fn sea(&self, lat: RangeInclusive<f64>, lon0: f64, lon1: f64) -> Vec<Cell> {
        cells_in(&self.grid, lat, lon0, lon1)
            .into_iter()
            .filter(|&c| self.elev[c as usize] < 0.0)
            .collect()
    }

    fn land(&self) -> Vec<Cell> {
        self.grid
            .cells()
            .filter(|&c| self.elev[c as usize] >= 0.0)
            .collect()
    }

    /// The cell at a place, or the nearest land cell (a coastal station whose cell is sea).
    fn at_land(&self, lat: f64, lon: f64) -> usize {
        near(&self.grid, lat, lon, 2, |c| self.elev[c as usize] >= 0.0) as usize
    }

    /// The land cell within a degree of a station whose height is nearest the station's: a
    /// one-degree cell averages over valleys and ranges, and Lima's own cell is mostly Andes.
    fn station_cell(&self, st: &Station) -> usize {
        let near = cells_in(&self.grid, (st.lat - 1.0)..=(st.lat + 1.0), st.lon - 1.0, st.lon + 1.0);
        near.into_iter()
            .filter(|&c| self.elev[c as usize] >= 0.0)
            .min_by(|&a, &b| (self.elev[a as usize] - st.elev_m).abs().total_cmp(&(self.elev[b as usize] - st.elev_m).abs()))
            .map_or_else(|| self.at_land(st.lat, st.lon), |c| c as usize)
    }

    fn station(&self, name: &str) -> (&'static Station, usize) {
        let st = STATIONS.iter().find(|st| st.name == name).unwrap();
        (st, self.station_cell(st))
    }

    /// Area-weighted mean (east, north) of a vector field.
    fn vector_mean(&self, f: &Field<[f32; 2]>, cells: &[Cell]) -> (f64, f64) {
        (
            mean(&self.grid, cells, |c| f64::from(f.get(c)[0])),
            mean(&self.grid, cells, |c| f64::from(f.get(c)[1])),
        )
    }

    /// The commonest land biome, by area, within 1.5° of a place (lakes left out: they are
    /// hydrology's).
    fn dominant_biome(&self, lat: f64, lon: f64) -> String {
        let b = self.p.biome().unwrap();
        let mut area: BTreeMap<String, f64> = BTreeMap::new();
        for c in cells_in(&self.grid, (lat - 1.5)..=(lat + 1.5), lon - 1.5, lon + 1.5) {
            if let Some(r) = b
                .row_at(c)
                .filter(|r| self.elev[c as usize] >= 0.0 && r.slug != "lake")
            {
                *area.entry(r.slug.clone()).or_default() += self.grid.area_sr(c);
            }
        }
        area.into_iter()
            .max_by(|a, b| a.1.total_cmp(&b.1))
            .map_or("none".into(), |(k, _)| k)
    }
}

fn winds(s: &mut Scorecard, e: &Earth) {
    let (jan, jul) = (e.p.wind_jan().unwrap(), e.p.wind_jul().unwrap());
    let east = |w, cells: &[Cell]| e.vector_mean(w, cells).0;
    let north = |w, cells: &[Cell]| e.vector_mean(w, cells).1;
    for (name, w) in [("January", jan), ("July", jul)] {
        s.range(
            "winds",
            format!("{name} trades, 10–20°N, east–west"),
            east(w, &e.sea(10.0..=20.0, -180.0, 180.0)),
            -10.0..=-1.0,
            "m/s",
        );
        s.range(
            "winds",
            format!("{name} trades, 10–20°S, east–west"),
            east(w, &e.sea(-20.0..=-10.0, -180.0, 180.0)),
            -10.0..=-1.0,
            "m/s",
        );
        s.range(
            "winds",
            format!("{name} westerlies, 40–55°S, east–west"),
            east(w, &e.sea(-55.0..=-40.0, -180.0, 180.0)),
            2.0..=15.0,
            "m/s",
        );
    }
    s.range(
        "winds",
        "January westerlies, 40–55°N, east–west",
        east(jan, &e.sea(40.0..=55.0, -180.0, 180.0)),
        2.0..=15.0,
        "m/s",
    );
    let arabian = e.sea(8.0..=18.0, 55.0, 72.0);
    s.range(
        "winds",
        "July southwest monsoon, Arabian Sea, east–west",
        east(jul, &arabian),
        2.0..=15.0,
        "m/s",
    );
    s.range(
        "winds",
        "January northeast monsoon, Arabian Sea, east–west",
        east(jan, &arabian),
        -10.0..=0.0,
        "m/s",
    );
    s.range(
        "winds",
        "July West African monsoon, 5–12°N 10°W–10°E, northward",
        north(jul, &cells_in(&e.grid, 5.0..=12.0, -10.0, 10.0)),
        0.5..=10.0,
        "m/s",
    );
    s.range(
        "winds",
        "January South China Sea, northward",
        north(jan, &e.sea(10.0..=20.0, 110.0, 120.0)),
        -10.0..=-0.5,
        "m/s",
    );
}

fn ocean(s: &mut Scorecard, e: &Earth) {
    let cur = e.p.ocean_currents().unwrap();
    let flow = |lat, lon0, lon1| e.vector_mean(cur, &e.sea(lat, lon0, lon1));
    s.range(
        "currents",
        "Gulf Stream, 30–38°N 80–72°W, northward",
        flow(30.0..=38.0, -80.0, -72.0).1,
        0.1..=2.0,
        "m/s",
    );
    s.range(
        "currents",
        "Kuroshio, 25–35°N 125–140°E, northward",
        flow(25.0..=35.0, 125.0, 140.0).1,
        0.05..=2.0,
        "m/s",
    );
    s.range(
        "currents",
        "California Current, 25–40°N 130–120°W, northward",
        flow(25.0..=40.0, -130.0, -120.0).1,
        -1.0..=-0.01,
        "m/s",
    );
    s.range(
        "currents",
        "Humboldt Current, 30–10°S 85–75°W, northward",
        flow(-30.0..=-10.0, -85.0, -75.0).1,
        0.01..=1.0,
        "m/s",
    );
    s.range(
        "currents",
        "Benguela Current, 30–15°S 5–15°E, northward",
        flow(-30.0..=-15.0, 5.0, 15.0).1,
        0.01..=1.0,
        "m/s",
    );
    s.range(
        "currents",
        "Antarctic Circumpolar, 50–60°S, eastward",
        flow(-60.0..=-50.0, -180.0, 180.0).0,
        0.05..=1.0,
        "m/s",
    );
    s.range(
        "currents",
        "South Equatorial, 3°S–3°N 150–100°W, eastward",
        flow(-3.0..=3.0, -150.0, -100.0).0,
        -1.0..=-0.02,
        "m/s",
    );

    let sst = e.p.sst_anomaly().unwrap();
    let anom = |lat, lon0, lon1| mean(&e.grid, &e.sea(lat, lon0, lon1), |c| f64::from(*sst.get(c)));
    s.range(
        "sea surface",
        "Gulf Stream, 35–42°N 72–55°W, warmer than the zonal mean",
        anom(35.0..=42.0, -72.0, -55.0),
        0.0..=10.0,
        "°C",
    );
    // Observed: annual SST departure from the zonal mean (scripts/real-earth/check.py).
    #[rustfmt::skip]
    let observed = [
        ("Kuroshio, 30–40°N 130–150°E", 30.0..=40.0, 130.0, 150.0, 0.2),
        ("California Current, 25–40°N 130–115°W", 25.0..=40.0, -130.0, -115.0, -4.2),
        ("Humboldt, 25–5°S 85–72°W", -25.0..=-5.0, -85.0, -72.0, -4.9),
        ("Benguela, 30–15°S 5–14°E", -30.0..=-15.0, 5.0, 14.0, -4.4),
        ("Canary, 15–30°N 25–12°W", 15.0..=30.0, -25.0, -12.0, -3.0),
        ("NE Atlantic, 45–65°N 30°W–10°E", 45.0..=65.0, -30.0, 10.0, 4.4),
    ];
    for (name, lat, lon0, lon1, obs) in observed {
        // Within 2.5 °C of observed, on the observed side of zero.
        let (lo, hi) = if obs > 0.0 { (f64::max(0.0, obs - 2.5), obs + 2.5) } else { (obs - 2.5, f64::min(0.0, obs + 2.5)) };
        s.range("sea surface", format!("{name} (observed {obs:+.1})"), anom(lat, lon0, lon1), lo..=hi, "°C");
    }
}

fn temperature(s: &mut Scorecard, e: &Earth) {
    let g = &e.grid;
    let all: Vec<Cell> = g.cells().collect();
    let zone =
        |lo: f64, hi: f64| mean(g, &cells_in(g, lo..=hi, -180.0, 180.0), |c| e.t[c as usize]);
    s.range(
        "temperature",
        "global mean",
        mean(g, &all, |c| e.t[c as usize]),
        13.5..=15.5,
        "°C",
    );
    s.range(
        "temperature",
        "equator, 5°S–5°N",
        zone(-5.0, 5.0),
        25.0..=28.0,
        "°C",
    );
    s.range(
        "temperature",
        "75–90°N",
        zone(75.0, 90.0),
        -25.0..=-10.0,
        "°C",
    );
    s.range(
        "temperature",
        "75–90°S",
        zone(-90.0, -75.0),
        -50.0..=-35.0,
        "°C",
    )
    .gap(ANTARCTICA);
    for st in STATIONS {
        let c = e.station_cell(st);
        // Compare at the station's height: the cell can be a few hundred metres off.
        let at_station = e.t[c] + 6.5 * (e.elev[c].max(0.0) - st.elev_m) / 1000.0;
        s.range(
            "temperature",
            format!("{} mean (observed {:.1})", st.name, st.t),
            at_station,
            (st.t - 4.0)..=(st.t + 4.0),
            "°C",
        )
        .gap_if_listed(st.name, TEMPERATURE_GAPS);
    }
    for (name, lo, hi) in [
        ("Singapore", 0.0, 4.0),
        ("Manaus", 0.0, 4.0),
        ("London", 9.0, 18.0),
        ("Perth", 7.0, 15.0),
        ("Moscow", 18.0, 32.0),
        ("Beijing", 22.0, 36.0),
        ("Yakutsk", 45.0, 65.0),
    ] {
        let (st, c) = e.station(name);
        s.range(
            "temperature",
            format!("{name} annual range (observed {:.0})", st.range),
            e.t_range[c],
            lo..=hi,
            "°C",
        )
        .gap_if_listed(name, RANGE_GAPS);
    }
}

fn precipitation(s: &mut Scorecard, e: &Earth) {
    let g = &e.grid;
    for st in STATIONS {
        let v = e.precip[e.station_cell(st)];
        // Within a factor of two or so; a hyperarid station may read zero.
        let lo = if st.precip < 50.0 {
            0.0
        } else {
            st.precip * 0.5
        };
        s.range(
            "precipitation",
            format!("{} (observed {:.0})", st.name, st.precip),
            v,
            lo..=(st.precip * 1.8 + 60.0),
            "mm",
        )
        .gap_if_listed(st.name, RAIN_GAPS);
    }
    for name in ["Los Angeles", "Perth", "Cape Town"] {
        s.range(
            "precipitation",
            format!("{name}: dry summers, warm-half share"),
            e.warm_share[e.station(name).1],
            0.0..=0.35,
            "",
        );
    }
    for name in ["Mumbai", "Beijing", "Timbuktu"] {
        s.range(
            "precipitation",
            format!("{name}: summer monsoon, warm-half share"),
            e.warm_share[e.station(name).1],
            0.7..=1.0,
            "",
        );
    }
    s.range(
        "precipitation",
        "land mean (Earth 780)",
        mean(g, &e.land(), |c| e.precip[c as usize]),
        650.0..=950.0,
        "mm",
    );
    let both = |lo: f64, hi: f64| {
        [
            e.sea(lo..=hi, -180.0, 180.0),
            e.sea(-hi..=-lo, -180.0, 180.0),
        ]
        .concat()
    };
    s.range(
        "precipitation",
        "subtropical oceans, 20–30° (Earth ~700)",
        mean(g, &both(20.0, 30.0), |c| e.precip[c as usize]),
        450.0..=1000.0,
        "mm",
    );
    s.range(
        "precipitation",
        "equatorial oceans, 5°S–10°N (Earth ~2000)",
        mean(g, &e.sea(-5.0..=10.0, -180.0, 180.0), |c| {
            e.precip[c as usize]
        }),
        1500.0..=3000.0,
        "mm",
    );

    let matches = STATIONS
        .iter()
        .filter(|st| {
            let c = e.station_cell(st);
            koppen_group(e.t[c], e.t_range[c], e.precip[c], e.warm_share[c])
                == observed_group(st.koppen)
        })
        .count();
    s.range(
        "climate class",
        "stations in their Köppen group",
        matches as f64 / STATIONS.len() as f64,
        0.65..=1.0,
        "",
    )
    .gap("follows the station misses above");
}

fn hydrology(s: &mut Scorecard, e: &Earth) {
    let g = &e.grid;
    let rivers = e.p.rivers().unwrap();
    let biggest = rivers
        .nodes
        .rows
        .values()
        .max_by(|a, b| a.discharge_m3_s.total_cmp(&b.discharge_m3_s))
        .unwrap();
    let (lat, lon) = known::lat_lon(g, biggest.cell);
    let at_amazon = (lat + 0.5).abs() < 6.0 && (lon + 50.0).abs() < 8.0;
    s.one_of(
        "hydrology",
        "largest river reaches the sea near the Amazon mouth",
        if at_amazon { "Amazon" } else { "elsewhere" },
        &["Amazon"],
    );
    let sheds = e.p.watersheds().unwrap();
    let shed = |lat: f64, lon: f64| sheds.row_at(e.at_land(lat, lon) as Cell);
    let to_med =
        shed(15.6, 32.5).is_some_and(|w| w.reaches_sea && known::lat_lon(g, w.outlet).0 > 28.0);
    s.one_of(
        "hydrology",
        "Khartoum drains to the Mediterranean",
        if to_med { "yes" } else { "no" },
        &["yes"],
    );
    for (name, lat, lon) in [
        ("Great Basin", 40.0, -117.0),
        ("Tarim", 40.0, 83.0),
        ("Lake Eyre", -28.0, 137.5),
        ("Lake Chad", 14.0, 16.0),
    ] {
        let closed = shed(lat, lon).is_some_and(|w| !w.reaches_sea);
        s.one_of(
            "hydrology",
            format!("{name} is a closed basin"),
            if closed { "closed" } else { "drains to sea" },
            &["closed"],
        )
        .gap_if_listed(name, &[("Great Basin", SOUTHWEST)]);
    }
}

fn biomes(s: &mut Scorecard, e: &Earth) {
    let biome = e.p.biome().unwrap();
    let land = e.land();
    let share = |slug: &str| {
        mean(&e.grid, &land, |c| {
            f64::from(u8::from(biome.row_at(c).is_some_and(|r| r.slug == slug)))
        })
    };
    for (slug, lo, hi) in [
        ("hot_desert", 0.10, 0.20),
        ("temperate_grassland", 0.04, 0.13),
        ("ice_sheet", 0.08, 0.13),
        ("tropical_rainforest", 0.05, 0.12),
        ("boreal_forest", 0.08, 0.15),
        ("tundra", 0.04, 0.10),
        ("savanna", 0.10, 0.20),
        ("mediterranean_woodland", 0.01, 0.04),
    ] {
        s.range("biomes", format!("{slug}, share of land"), share(slug), lo..=hi, "")
            .gap_if_listed(slug, &[("savanna", "the Sahel comes out rainforest (the West African monsoon reaches too far north), so savanna runs short")]);
    }
    for (name, lat, lon, want) in [
        ("Sahara", 23.0, 10.0, &["hot_desert"][..]),
        ("Arabia", 22.0, 47.0, &["hot_desert"]),
        ("Australian outback", -25.0, 133.0, &["hot_desert"]),
        ("Amazon", -3.0, -62.0, &["tropical_rainforest"]),
        (
            "Congo",
            0.0,
            22.0,
            &["tropical_rainforest", "tropical_seasonal_forest"],
        ),
        (
            "Cerrado",
            -15.0,
            -50.0,
            &["savanna", "tropical_seasonal_forest"],
        ),
        ("Sahel", 13.0, 0.0, &["savanna"]),
        ("Siberian taiga", 60.0, 100.0, &["boreal_forest"]),
        ("Canadian taiga", 55.0, -100.0, &["boreal_forest"]),
        ("Canadian Arctic", 68.0, -100.0, &["tundra"]),
        ("Greenland", 75.0, -40.0, &["ice_sheet"]),
        ("Antarctica", -80.0, 0.0, &["ice_sheet"]),
        ("Great Plains", 40.0, -100.0, &["temperate_grassland"]),
        ("Kazakh steppe", 50.0, 70.0, &["temperate_grassland"]),
        ("Pampas", -35.0, -62.0, &["temperate_grassland", TF]),
        ("Central Europe", 50.0, 10.0, &[TF]),
        ("Eastern US", 38.0, -82.0, &[TF]),
        ("Greece", 38.0, 22.0, &["mediterranean_woodland"]),
        ("California", 37.0, -121.0, &["mediterranean_woodland"]),
        ("Tibet", 33.0, 88.0, &["alpine", "cold_desert", "tundra"]),
    ] {
        s.one_of("biomes", name, &e.dominant_biome(lat, lon), want)
            .gap_if_listed(name, PLACE_GAPS);
    }
}

/// The classifier alone, fed station normals: its own errors, apart from the climate's.
fn classifier(s: &mut Scorecard) {
    let p = BiomesParams::default();
    for st in STATIONS {
        let c = CellClimate {
            t_mean: st.t,
            t_range: st.range,
            precip_mm: st.precip,
            elevation_m: st.elev_m,
            lat_deg: st.lat,
            warm_share: st.warm_share,
        };
        s.one_of(
            "classifier",
            format!("{} from its normals", st.name),
            &classify(&c, &p).class().slug,
            st.biomes,
        )
        .gap_if_listed(st.name, CLASSIFIER_GAPS);
    }
}

fn scorecard() -> &'static Scorecard {
    static S: OnceLock<Scorecard> = OnceLock::new();
    S.get_or_init(|| {
        let e = earth();
        known::dump("climate_temperature", e.t.iter().copied());
        known::dump("climate_precipitation", e.precip.iter().copied());
        known::dump("climate_warm_share", e.warm_share.iter().copied());
        known::dump("climate_t_range", e.t_range.iter().copied());
        known::dump("circulation_sst_anomaly", e.p.sst_anomaly().unwrap().values.iter().map(|&v| f64::from(v)));
        let (jul, jan) = (e.p.wind_jul().unwrap(), e.p.wind_jan().unwrap());
        known::dump("circulation_wind_jul", jul.values.iter().flat_map(|w| [f64::from(w[0]), f64::from(w[1])]));
        known::dump("circulation_wind_jan", jan.values.iter().flat_map(|w| [f64::from(w[0]), f64::from(w[1])]));
        let mut s = Scorecard::default();
        winds(&mut s, e);
        ocean(&mut s, e);
        temperature(&mut s, e);
        precipitation(&mut s, e);
        hydrology(&mut s, e);
        biomes(&mut s, e);
        classifier(&mut s);
        s
    })
}

#[test]
fn climate_known_cases() {
    let s = scorecard();
    println!("{}", s.report());
    s.assert_cases("earth_climate");
}

#[test]
#[ignore = "known gaps: run with --ignored"]
fn climate_known_gaps() {
    scorecard().assert_gaps();
}
