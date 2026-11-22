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

use crate::ground::Ground;
use geoforge_core::grid::Cell;
use geoforge_core::seed::Seed;
use geoforge_core::surface::GeologicProvince as P;

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

const OROGENS: &[P] = &[P::CollisionOrogen, P::PaleoOrogen, P::OrogenicPlateau];
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
    Setting { name: "oolitic ironstone", rate: 0.1, applies: |g, c| on(land(g, c) && matches!(g.rock[c], "carbonate" | "clastic_sedimentary") && g.province_is(c, BASINS) && (50.0..500.0).contains(&age(g, c))) },
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
        applies: |g, c| on(land(g, c) && matches!(g.rock[c], "rift_fill" | "clastic_sedimentary") && g.province_is(c, &[P::ContinentalRift, P::BackarcBasin, P::IntracratonicBasin]) && age(g, c) > 500.0),
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
    Setting { name: "foreland basin coal (high grade)", rate: 0.8, applies: |g, c| on(land(g, c) && g.rock[c] == "clastic_sedimentary" && g.province_is(c, &[P::ForelandBasin]) && (2.0..360.0).contains(&age(g, c))) },
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
