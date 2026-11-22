//! Where resources are likely: prospectivity, soils and timber, worked out from the stored
//! geology, climate and rivers (derived views, never stored).
//!
//! A global map cannot place individual mines, so each deposit kind gets λ, the expected
//! number of deposits per 10,000 km² (a 100 × 100 km region), summed over the settings it
//! forms in. Exported as P = 1 − e^(−λ), the chance such a region holds at least one; a
//! region zoomed into for local work draws its deposits as a Poisson process with that λ,
//! seeded by place ([`draw_count`]). The rates are orders of magnitude, scaled to Earth's
//! known significant deposits; the settings follow deposit geology (porphyry copper in arcs,
//! banded iron in Archean–Paleoproterozoic cratons, tin in collisional granites, orogenic gold
//! in greenstone belts) and Artifexian's Logs 45–48 (peat, bog iron, salt flats, supergene
//! enrichment where winters stay above 6 °C, the tree line).
//!
//! Coal needs land plants: on Earth, swamps of the last ~360 Myr. Other worlds are given the
//! same window.
//!
//! The one stored part is the **giant deposits** ([`GiantDeposits`], the resources layer): a
//! few per kind, as many as one in 50 of the kind's expected deposits, placed where the ground
//! favours the kind, with a size (log-normal about the kind's typical giant), a depth and the
//! host rock's age. Edit them like any table.

use crate::ground::{self, Ground};
use geoforge_core::grid::Cell;
use geoforge_core::key::Key;
use geoforge_core::product::{Generator, LayerId, Precondition, ProductData, ProductId, Products, RunContext};
use geoforge_core::seed::{unit, Seed};
use geoforge_core::shapes::Table;
use geoforge_core::surface::{GeologicProvince as P, GiantDeposit};
use geoforge_core::{BodyContext, Error};
use serde::{Deserialize, Serialize};

/// A kind of deposit.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Deposit {
    Iron,
    Copper,
    Tin,
    LeadZincSilver,
    Gold,
    Coal,
    Salt,
    Peat,
    Gemstones,
}

impl Deposit {
    pub const ALL: [Self; 9] = [Self::Iron, Self::Copper, Self::Tin, Self::LeadZincSilver, Self::Gold, Self::Coal, Self::Salt, Self::Peat, Self::Gemstones];

    #[must_use]
    pub fn slug(self) -> &'static str {
        match self {
            Self::Iron => "iron",
            Self::Copper => "copper",
            Self::Tin => "tin",
            Self::LeadZincSilver => "lead_zinc_silver",
            Self::Gold => "gold",
            Self::Coal => "coal",
            Self::Salt => "salt",
            Self::Peat => "peat",
            Self::Gemstones => "gemstones",
        }
    }

    #[must_use]
    pub fn from_slug(s: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|d| d.slug() == s)
    }

    /// The settings it forms in.
    #[must_use]
    pub fn settings(self) -> &'static [Setting] {
        match self {
            Self::Iron => IRON,
            Self::Copper => COPPER,
            Self::Tin => TIN,
            Self::LeadZincSilver => LEAD_ZINC_SILVER,
            Self::Gold => GOLD,
            Self::Coal => COAL,
            Self::Salt => SALT,
            Self::Peat => PEAT,
            Self::Gemstones => GEMSTONES,
        }
    }
}

/// One way a deposit forms: expected deposits per 10,000 km² where it applies, times how
/// much it applies at a cell (0 none, 1 fully, more where conditions enrich it).
pub struct Setting {
    pub name: &'static str,
    pub rate: f64,
    applies: fn(&Ground, usize) -> f64,
}

/// Placers lie in rivers up to this far from their source, km.
const PLACER_KM: f64 = 300.0;

fn on(b: bool) -> f64 {
    if b {
        1.0
    } else {
        0.0
    }
}

fn age(g: &Ground, c: usize) -> f64 {
    f64::from(g.rock_age[c])
}

fn land(g: &Ground, c: usize) -> bool {
    !g.sea[c]
}

const OROGENS: &[P] = &[P::CollisionOrogen, P::PaleoOrogen, P::OrogenicPlateau, P::BasementUplift];
const BASINS: &[P] = &[P::Platform, P::IntracratonicBasin, P::ForelandBasin];

/// Supergene enrichment: weathering concentrates copper where water percolates seasonally and
/// winters stay above 6 °C (snow stops it): arid and seasonally dry climates.
fn supergene(g: &Ground, c: usize) -> f64 {
    let dry_season = matches!(g.group(c), 'B') || matches!(g.code(c), "Aw" | "Csa" | "Csb" | "Cwa" | "Cwb");
    if dry_season && g.coldest(c) > 6.0 {
        2.0
    } else {
        1.0
    }
}

const IRON: &[Setting] = &[
    Setting { name: "banded iron formation", rate: 0.4, applies: |g, c| on(land(g, c) && matches!(g.rock[c], "greenstone" | "cratonic_basement") && (1800.0..2600.0).contains(&age(g, c))) },
    Setting { name: "skarn", rate: 0.3, applies: |g, c| on(land(g, c) && g.rock[c] == "granitoid" && g.touches_rock(c, "carbonate")) },
    Setting { name: "laterite", rate: 0.4, applies: |g, c| on(land(g, c) && matches!(g.rock[c], "flood_basalt" | "ophiolite") && matches!(g.code(c), "Aw" | "Am")) },
    Setting {
        name: "oolitic ironstone",
        rate: 0.1,
        applies: |g, c| on(land(g, c) && matches!(g.rock[c], "carbonate" | "clastic_sedimentary") && g.province_is(c, BASINS) && (50.0..500.0).contains(&age(g, c))),
    },
    Setting {
        name: "bog iron",
        rate: 0.3,
        applies: |g, c| on(land(g, c) && (g.biome[c] == "wetland" || (g.elevation[c] < 300.0 && g.slope[c] < 1.0 && matches!(g.code(c), "Dfb" | "Dfc" | "Cfb" | "Cfc" | "Af")))),
    },
    Setting { name: "meteoric iron", rate: 0.002, applies: |g, c| on(matches!(g.code(c), "BWh" | "BWk" | "ET")) },
];

const COPPER: &[Setting] = &[
    Setting { name: "porphyry", rate: 0.6, applies: |g, c| on(land(g, c) && (g.province_is(c, &[P::VolcanicArc]) || g.rock[c] == "arc_volcanic") && age(g, c) < 200.0) * supergene(g, c) },
    Setting { name: "volcanogenic massive sulfide", rate: 0.25, applies: |g, c| on(land(g, c) && matches!(g.rock[c], "ophiolite" | "accretionary_melange" | "greenstone")) * supergene(g, c) },
    Setting { name: "native copper in flood basalt", rate: 0.1, applies: |g, c| on(land(g, c) && g.rock[c] == "flood_basalt") },
    Setting {
        name: "sediment-hosted",
        rate: 0.15,
        applies: |g, c| on(land(g, c) && (g.rock[c] == "rift_fill" || (g.rock[c] == "clastic_sedimentary" && g.province_is(c, &[P::IntracratonicBasin, P::ContinentalRift])))),
    },
    Setting { name: "skarn", rate: 0.2, applies: |g, c| on(land(g, c) && g.rock[c] == "granitoid" && g.touches_rock(c, "carbonate")) },
    Setting { name: "iron oxide–copper–gold", rate: 0.05, applies: |g, c| on(land(g, c) && g.rock[c] == "cratonic_basement" && age(g, c) > 1000.0) },
];

const TIN: &[Setting] = &[Setting { name: "collisional granite", rate: 0.4, applies: |g, c| on(land(g, c) && g.rock[c] == "granitoid" && g.province_is(c, OROGENS)) }];

const LEAD_ZINC_SILVER: &[Setting] = &[
    Setting { name: "Mississippi Valley-type", rate: 0.2, applies: |g, c| on(land(g, c) && g.rock[c] == "carbonate" && g.province_is(c, BASINS)) },
    Setting {
        name: "sedimentary exhalative",
        rate: 0.1,
        applies: |g, c| {
            on(land(g, c) && matches!(g.rock[c], "rift_fill" | "clastic_sedimentary") && g.province_is(c, &[P::ContinentalRift, P::BackarcBasin, P::IntracratonicBasin]) && age(g, c) > 500.0)
        },
    },
    Setting { name: "epithermal silver", rate: 0.3, applies: |g, c| on(land(g, c) && g.rock[c] == "arc_volcanic" && age(g, c) < 50.0) },
    Setting { name: "volcanogenic massive sulfide", rate: 0.15, applies: |g, c| on(land(g, c) && matches!(g.rock[c], "ophiolite" | "accretionary_melange" | "greenstone")) },
];

const GOLD: &[Setting] = &[
    Setting { name: "orogenic gold in greenstone", rate: 1.0, applies: |g, c| on(land(g, c) && g.rock[c] == "greenstone") * if g.group(c) == 'A' { 1.3 } else { 1.0 } },
    Setting { name: "orogenic gold in old belts", rate: 0.2, applies: |g, c| on(land(g, c) && g.rock[c] == "metamorphic" && g.province_is(c, OROGENS)) },
    Setting { name: "epithermal and porphyry gold", rate: 0.3, applies: |g, c| on(land(g, c) && (g.province_is(c, &[P::VolcanicArc]) || g.rock[c] == "arc_volcanic") && age(g, c) < 100.0) },
];

const COAL: &[Setting] = &[
    Setting {
        name: "foreland basin coal (high grade)",
        rate: 0.8,
        applies: |g, c| on(land(g, c) && g.rock[c] == "clastic_sedimentary" && g.province_is(c, &[P::ForelandBasin]) && (2.0..360.0).contains(&age(g, c))),
    },
    Setting {
        name: "basin coal and lignite",
        rate: 0.2,
        applies: |g, c| on(land(g, c) && g.rock[c] == "clastic_sedimentary" && g.province_is(c, &[P::IntracratonicBasin, P::Platform]) && (2.0..360.0).contains(&age(g, c))),
    },
];

const SALT: &[Setting] = &[
    Setting { name: "rock salt", rate: 1.5, applies: |g, c| on(land(g, c) && g.rock[c] == "evaporite") },
    Setting { name: "salt flats", rate: 2.0, applies: |g, c| on(g.salt_lake[c] || (land(g, c) && g.closed_basin[c] && g.group(c) == 'B' && g.slope[c] < 0.5)) },
    Setting { name: "sea salt by the sun", rate: 1.0, applies: |g, c| on(g.coast(c) && (g.group(c) == 'B' || g.code(c).starts_with("Cs"))) },
];

/// Peat: cold or cool, wet, flat ground (Artifexian keys it to non-dry-summer temperate and
/// continental climates and to ground glaciated in the last ice age, which we do not store).
fn peaty(g: &Ground, c: usize) -> bool {
    let wet_enough = if g.group(c) == 'C' { 700.0 } else { 400.0 };
    land(g, c) && g.elevation[c] < 600.0 && g.slope[c] < 2.0 && f64::from(g.precipitation[c]) > wet_enough && matches!(g.code(c), "Dfb" | "Dfc" | "Dfd" | "Dwc" | "Dwd" | "Cfb" | "Cfc" | "ET")
}

const PEAT: &[Setting] = &[
    Setting { name: "peatland", rate: 0.8, applies: |g, c| on(peaty(g, c)) },
    Setting { name: "fen and bog in wetlands", rate: 2.0, applies: |g, c| on(land(g, c) && g.biome[c] == "wetland" && g.warmest(c) < 25.0) },
];

const GEMSTONES: &[Setting] = &[
    Setting { name: "kimberlite diamonds", rate: 0.05, applies: |g, c| on(land(g, c) && g.province_is(c, &[P::Craton]) && age(g, c) > 2500.0) },
    Setting { name: "pegmatite gems", rate: 0.05, applies: |g, c| on(land(g, c) && g.rock[c] == "granitoid" && g.province_is(c, OROGENS)) },
];

/// λ for `kind` at every cell: expected deposits per 10,000 km².
#[must_use]
pub fn prospectivity(g: &Ground, kind: Deposit) -> Vec<f64> {
    let mut lambda: Vec<f64> = (0..g.len()).map(|c| kind.settings().iter().map(|s| s.rate * (s.applies)(g, c)).sum()).collect();
    // Tin and gold wash out of their hills into the rivers below: placers.
    if matches!(kind, Deposit::Tin | Deposit::Gold) {
        let sources: Vec<Cell> = (0..g.len()).filter(|&c| lambda[c] > 0.0).map(|c| c as Cell).collect();
        let near = g.distance_over_land_km(sources, PLACER_KM);
        for c in 0..g.len() {
            if g.river_order[c] >= 3 && near[c].is_finite() {
                lambda[c] += PLACER_RATE * (1.0 - near[c] / PLACER_KM);
            }
        }
    }
    lambda
}

/// Placers right below their source, fading to nothing [`PLACER_KM`] downstream.
pub const PLACER_RATE: f64 = 0.4;

/// The settings that apply at a cell and their λ (placers not included).
#[must_use]
pub fn settings_at(g: &Ground, kind: Deposit, c: usize) -> Vec<(&'static str, f64)> {
    kind.settings().iter().map(|s| (s.name, s.rate * (s.applies)(g, c))).filter(|s| s.1 > 0.0).collect()
}

/// How many deposits a region holds: a Poisson draw with mean `mean`, seeded by `seed` (the
/// world's seed, the kind and the place), so every zoom into the same place agrees.
#[must_use]
pub fn draw_count(mean: f64, seed: Seed) -> u32 {
    // Knuth's method: fine for the small means of a region.
    let limit = (-mean.max(0.0)).exp();
    let (mut k, mut p) = (0u32, 1.0);
    loop {
        p *= geoforge_core::seed::unit(seed.child_u64(b"poisson", u64::from(k)));
        if p <= limit || k > 1000 {
            return k;
        }
        k += 1;
    }
}

/// Typical ore of a giant deposit of each kind, million tonnes (the median; sizes spread
/// log-normally about it, a factor e either way).
fn giant_ore_mt(kind: Deposit) -> f64 {
    match kind {
        Deposit::Iron => 2000.0,
        Deposit::Copper => 1500.0,
        Deposit::Tin => 50.0,
        Deposit::LeadZincSilver => 150.0,
        Deposit::Gold => 80.0,
        Deposit::Coal => 10_000.0,
        Deposit::Salt => 20_000.0,
        Deposit::Peat => 500.0,
        Deposit::Gemstones => 20.0,
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct DepositsParams {
    /// One giant per this many expected deposits of a kind.
    pub deposits_per_giant: f64,
    /// At most this many giants per kind.
    pub most_per_kind: u32,
}

impl Default for DepositsParams {
    fn default() -> Self {
        Self { deposits_per_giant: 50.0, most_per_kind: 12 }
    }
}

/// The resources layer: the giant deposits.
pub struct GiantDeposits;

impl Generator for GiantDeposits {
    type Params = DepositsParams;
    const NAME: &'static str = "giant_deposits";
    const VERSION: u32 = 1;
    const LAYER: LayerId = LayerId::Resources;

    fn needs() -> &'static [ProductId] {
        &ground::READS
    }
    fn applies(body: &BodyContext) -> Result<(), Precondition> {
        geoforge_surface::preconditions::surface_water(body)
    }

    fn run(ctx: &RunContext<'_>, inputs: &Products, p: &DepositsParams) -> Result<Products, Error> {
        let mut out = Products::new();
        if inputs.deposits().is_some() {
            return Ok(out);
        }
        let g = Ground::read(inputs, ctx.grid, ctx.radius_km(), ctx.body.interior.has_surface_liquid_water()).ok_or(Error::MissingInput(ProductId::Biome))?;
        let mut table = Table::new();
        let seed = ctx.seed(ProductId::Deposits);
        for kind in Deposit::ALL {
            for d in giants(&g, kind, seed.child(kind.slug().as_bytes()), p) {
                table.insert(Key(table.rows.len() as u64 + 1), d);
            }
        }
        out.insert(ProductData::Deposits(table));
        Ok(out)
    }
}

/// The giant deposits of `kind`: placed cell by cell in proportion to λ × area.
fn giants(g: &Ground, kind: Deposit, seed: Seed, p: &DepositsParams) -> Vec<GiantDeposit> {
    let lambda = prospectivity(g, kind);
    let weight: Vec<f64> = (0..g.len()).map(|c| lambda[c] * g.area_km2(c) / 1e4).collect();
    let total: f64 = weight.iter().sum();
    let count = draw_count(total / p.deposits_per_giant.max(1.0), seed.child(b"count")).min(p.most_per_kind);
    let mut cumulative = Vec::with_capacity(weight.len());
    let mut run = 0.0;
    for w in &weight {
        run += w;
        cumulative.push(run);
    }
    (0..count)
        .map(|i| {
            let s = seed.child_u64(b"giant", u64::from(i));
            let target = unit(s.child(b"where")) * total;
            let c = cumulative.partition_point(|&x| x < target).min(g.len() - 1);
            let ll = g.grid.center(c as Cell).to_lat_lon();
            let setting = settings_at(g, kind, c).into_iter().max_by(|a, b| a.1.total_cmp(&b.1)).map_or("placer", |s| s.0);
            // Box–Muller for the log-normal size.
            let normal = (-2.0 * unit(s.child(b"size-a")).max(1e-12).ln()).sqrt() * (std::f64::consts::TAU * unit(s.child(b"size-b"))).cos();
            let surface = matches!(setting, "placer" | "salt flats" | "sea salt by the sun" | "peatland" | "fen and bog in wetlands" | "laterite" | "bog iron");
            GiantDeposit {
                kind: kind.slug().into(),
                setting: setting.into(),
                lat_deg: ll.lat_deg,
                lon_deg: ll.lon_deg,
                ore_mt: giant_ore_mt(kind) * normal.exp(),
                depth_m: if surface { 0.0 } else { -300.0 * unit(s.child(b"depth")).max(1e-12).ln() },
                age_ma: f64::from(g.rock_age[c]),
            }
        })
        .collect()
}

/// Soil classes (after the WRB's reference groups), with fertility for farming (0–1).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Soil {
    BareRockAndIce,
    Cryosol,
    Histosol,
    Andosol,
    Fluvisol,
    Aridisol,
    Chernozem,
    Ferralsol,
    Podzol,
    ForestSoil,
}

impl Soil {
    pub const ALL: [Self; 10] = [Self::BareRockAndIce, Self::Cryosol, Self::Histosol, Self::Andosol, Self::Fluvisol, Self::Aridisol, Self::Chernozem, Self::Ferralsol, Self::Podzol, Self::ForestSoil];

    #[must_use]
    pub fn code(self) -> &'static str {
        self.info().0
    }

    #[must_use]
    pub fn name(self) -> &'static str {
        self.info().1
    }

    #[must_use]
    pub fn color(self) -> [u8; 3] {
        self.info().2
    }

    /// Fertility for farming, 0–1.
    #[must_use]
    pub fn fertility(self) -> f64 {
        self.info().3
    }

    fn info(self) -> (&'static str, &'static str, [u8; 3], f64) {
        match self {
            Self::BareRockAndIce => ("rock_ice", "Bare rock and ice", [225, 225, 230], 0.0),
            Self::Cryosol => ("cryosol", "Cryosol (permafrost)", [170, 185, 200], 0.05),
            Self::Histosol => ("histosol", "Histosol (peat)", [90, 70, 60], 0.3),
            Self::Andosol => ("andosol", "Andosol (young volcanic)", [140, 40, 40], 0.9),
            Self::Fluvisol => ("fluvisol", "Fluvisol (floodplain alluvium)", [60, 150, 200], 1.0),
            Self::Aridisol => ("aridisol", "Aridisol (desert)", [235, 210, 160], 0.1),
            Self::Chernozem => ("chernozem", "Chernozem (black grassland soil)", [40, 35, 30], 1.0),
            Self::Ferralsol => ("ferralsol", "Ferralsol (leached tropical)", [200, 90, 40], 0.3),
            Self::Podzol => ("podzol", "Podzol (acid forest soil)", [150, 150, 130], 0.3),
            Self::ForestSoil => ("luvisol", "Luvisol and cambisol (temperate forest soils)", [150, 110, 70], 0.7),
        }
    }
}

/// The soil at a cell, `None` at sea: from climate, rock, relief and rivers.
#[must_use]
pub fn soil(g: &Ground, c: usize) -> Option<Soil> {
    if g.sea[c] {
        return None;
    }
    let code = g.code(c);
    let volcanic = (matches!(g.rock[c], "arc_volcanic" | "flood_basalt") && age(g, c) < 10.0) || (g.province_is(c, &[P::VolcanicArc, P::HotspotTrack, P::ContinentalHotspotTrack]) && age(g, c) < 5.0);
    Some(match code {
        "EF" => Soil::BareRockAndIce,
        _ if g.biome[c] == "ice_sheet" || g.slope[c] > 60.0 => Soil::BareRockAndIce,
        "ET" => Soil::Cryosol,
        _ if g.temperature[c] < -2.0 => Soil::Cryosol,
        _ if peaty(g, c) && g.biome[c] == "wetland" => Soil::Histosol,
        _ if volcanic => Soil::Andosol,
        _ if g.river_order[c] >= 4 && g.slope[c] < 1.0 => Soil::Fluvisol,
        _ if g.group(c) == 'B' => Soil::Aridisol,
        _ if g.biome[c] == "temperate_grassland" => Soil::Chernozem,
        _ if g.group(c) == 'A' => Soil::Ferralsol,
        _ if g.biome[c] == "boreal_forest" || matches!(code, "Dfc" | "Dfd" | "Dwc" | "Dwd" | "Dsc" | "Dsd" | "Cfc") => Soil::Podzol,
        _ => Soil::ForestSoil,
    })
}

/// Tree line, m (Artifexian, Log 45): 4000 m to 30° latitude, falling 130 m a degree to 50°,
/// then 75 m a degree to sea level near 69°.
#[must_use]
pub fn tree_line_m(lat_deg: f64) -> f64 {
    let l = lat_deg.abs();
    if l <= 30.0 {
        4000.0
    } else if l <= 50.0 {
        4000.0 - 130.0 * (l - 30.0)
    } else {
        (1400.0 - 75.0 * (l - 50.0)).max(0.0)
    }
}

/// Standing timber, 0–1, by biome, below the tree line and outside dry and polar climates.
#[must_use]
pub fn timber(g: &Ground, c: usize) -> Option<f64> {
    if g.sea[c] {
        return None;
    }
    let treeless = matches!(g.group(c), 'B' | 'E') || g.elevation[c] > tree_line_m(g.lat[c]);
    let density = match g.biome[c] {
        "tropical_rainforest" | "temperate_rainforest" => 1.0,
        "tropical_seasonal_forest" | "temperate_deciduous_forest" => 0.8,
        "boreal_forest" => 0.7,
        "mediterranean_woodland" => 0.4,
        "wetland" => 0.3,
        "savanna" => 0.2,
        "temperate_grassland" => 0.1,
        _ => 0.0,
    };
    Some(if treeless { 0.0 } else { density })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tree_line() {
        assert_eq!(tree_line_m(10.0), 4000.0);
        assert!((tree_line_m(46.0) - 1920.0).abs() < 1e-9, "the Alps: ~2000 m");
        assert_eq!(tree_line_m(-75.0), 0.0);
    }

    #[test]
    fn poisson_draws_are_seeded_and_average_right() {
        assert_eq!(draw_count(2.0, Seed(9)), draw_count(2.0, Seed(9)));
        let mean: f64 = (0..4000u64).map(|s| f64::from(draw_count(1.5, Seed(s)))).sum::<f64>() / 4000.0;
        assert!((mean - 1.5).abs() < 0.08, "{mean}");
        assert_eq!(draw_count(0.0, Seed(1)), 0);
    }

    #[test]
    fn kinds_round_trip() {
        for d in Deposit::ALL {
            assert_eq!(Deposit::from_slug(d.slug()), Some(d));
            assert!(!d.settings().is_empty());
        }
    }
}
