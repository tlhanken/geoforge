//! Surface lithology and rock age from geologic provinces.
//!
//! Each province kind has a characteristic mix of rock (an orogen exposes metamorphic cores,
//! granite plutons and deformed sediments; a platform is sediment over basement; an arc is
//! volcanic with plutonic roots). A smooth noise field picks within that mix, so rock types
//! form map-scale patches rather than speckle, and the same noise value means the same rock
//! across a province. Ages follow each rock's [`AgeRule`]: basement as old as the craton,
//! cover sediments younger, arc volcanics youngest.
//!
//! Rocks come from a [`Catalog`]: Earth's are built in with their province mixes, and
//! `[params.geology] extra_rocks` adds more (fantasy or otherwise), each saying where it
//! occurs and how old it is. The lithology partition's table is the catalog, so a painted
//! edit can use any of its rocks.
//!
//! Lithology feeds erosion (hard basement resists, soft sediment wears down) and, later,
//! resources (porphyry copper in arcs, gold and iron in greenstone, salt and hydrocarbons in
//! basins).

use crate::noise::{uniform_ranks, TerrainNoise};
use geoforge_core::grid::{Cell, Grid};
use geoforge_core::key::Key;
use geoforge_core::shapes::{Field, Partition, Table};
use geoforge_core::surface::{CrustType, GeologicProvince as G, Province, Rock};
use geoforge_core::Error;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

/// How old the exposures of a rock are.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum AgeRule {
    /// As old as the crust beneath (ocean floor and its cover).
    Crust,
    /// As old as the basement.
    Basement,
    /// Archean: 2.5 Ga or older.
    Archean,
    /// Plutons and metamorphic cores: the belt's age in orogens and arcs, part of the
    /// basement's elsewhere.
    Plutonic,
    /// Up to the province's age.
    Volcanic,
    /// The eruption's age on land; part of the crust's age at sea.
    Flood,
    /// Obducted ocean floor, somewhat older than its orogen.
    Ophiolite,
    /// Basin and wedge fill: most of the province's age.
    Infill,
    /// Sedimentary cover: older than the belt in orogens, young deposits widespread elsewhere.
    Cover,
}

/// Where a catalog rock occurs: a weight added to the mix of each matching setting. Built-in
/// mixes sum to 1, so 0.05 is about one part in twenty of such a province.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Occurrence {
    /// Province kind; absent = any, including ground outside named provinces.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub province: Option<G>,
    /// Crust type; absent = both.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub crust: Option<CrustType>,
    pub weight: f64,
}

/// A rock the geology layer can place.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RockDef {
    #[serde(flatten)]
    pub rock: Rock,
    pub age: AgeRule,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub occurs: Vec<Occurrence>,
}

/// The built-in rocks: Earth's, plus the impact, anorthosite and ice crusts of airless and
/// rigid-lid worlds. Keys are the variant number + 1, so stored lithology stays stable.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Earth {
    OceanicBasalt,
    PelagicSediment,
    CratonicBasement,
    Greenstone,
    Granitoid,
    Metamorphic,
    ArcVolcanic,
    FloodBasalt,
    ClasticSedimentary,
    Carbonate,
    Evaporite,
    Ophiolite,
    AccretionaryMelange,
    RiftFill,
    ImpactBreccia,
    Anorthosite,
    WaterIce,
}

impl Earth {
    const ALL: [Self; 17] = [
        Self::OceanicBasalt,
        Self::PelagicSediment,
        Self::CratonicBasement,
        Self::Greenstone,
        Self::Granitoid,
        Self::Metamorphic,
        Self::ArcVolcanic,
        Self::FloodBasalt,
        Self::ClasticSedimentary,
        Self::Carbonate,
        Self::Evaporite,
        Self::Ophiolite,
        Self::AccretionaryMelange,
        Self::RiftFill,
        Self::ImpactBreccia,
        Self::Anorthosite,
        Self::WaterIce,
    ];

    fn key(self) -> Key {
        Key(self as u64 + 1)
    }

    fn def(self) -> RockDef {
        use AgeRule as A;
        let (slug, name, color, erodibility, age) = match self {
            Self::OceanicBasalt => ("oceanic_basalt", "Oceanic basalt", [70, 80, 110], 0.35, A::Crust),
            Self::PelagicSediment => ("pelagic_sediment", "Pelagic sediment", [120, 140, 170], 1.0, A::Crust),
            Self::CratonicBasement => ("cratonic_basement", "Cratonic basement", [205, 120, 120], 0.25, A::Basement),
            Self::Greenstone => ("greenstone", "Greenstone belt", [70, 140, 90], 0.25, A::Archean),
            Self::Granitoid => ("granitoid", "Granitoid", [235, 150, 180], 0.25, A::Plutonic),
            Self::Metamorphic => ("metamorphic", "Metamorphic", [170, 110, 170], 0.35, A::Plutonic),
            Self::ArcVolcanic => ("arc_volcanic", "Arc volcanic", [220, 70, 60], 0.5, A::Volcanic),
            Self::FloodBasalt => ("flood_basalt", "Flood basalt", [110, 60, 60], 0.35, A::Flood),
            Self::ClasticSedimentary => ("clastic_sedimentary", "Clastic sedimentary", [225, 200, 130], 0.8, A::Cover),
            Self::Carbonate => ("carbonate", "Carbonate", [150, 200, 230], 0.5, A::Cover),
            Self::Evaporite => ("evaporite", "Evaporite", [245, 240, 225], 1.0, A::Cover),
            Self::Ophiolite => ("ophiolite", "Ophiolite", [40, 110, 100], 0.35, A::Ophiolite),
            Self::AccretionaryMelange => ("accretionary_melange", "Accretionary melange", [150, 150, 110], 0.8, A::Infill),
            Self::RiftFill => ("rift_fill", "Rift fill", [235, 165, 80], 1.0, A::Infill),
            Self::ImpactBreccia => ("impact_breccia", "Impact breccia", [140, 125, 115], 0.6, A::Crust),
            Self::Anorthosite => ("anorthosite", "Anorthosite", [225, 225, 215], 0.3, A::Basement),
            Self::WaterIce => ("water_ice", "Water ice", [215, 230, 240], 0.4, A::Crust),
        };
        RockDef { rock: Rock { slug: slug.into(), name: name.into(), color, erodibility, tags: Vec::new() }, age, occurs: Vec::new() }
    }

    /// The built-in rock mix of a setting (province kind, continental or oceanic crust).
    fn mix(kind: Option<G>, continental: bool) -> &'static [(Self, f64)] {
        use Earth as L;
        match (kind, continental) {
            (Some(G::Craton), true) => &[(L::CratonicBasement, 0.62), (L::Greenstone, 0.16), (L::Granitoid, 0.12), (L::ClasticSedimentary, 0.10)],
            (Some(G::Platform), true) | (None, true) => &[(L::ClasticSedimentary, 0.48), (L::Carbonate, 0.32), (L::Evaporite, 0.04), (L::CratonicBasement, 0.10), (L::Granitoid, 0.06)],
            (Some(G::IntracratonicBasin), true) => &[(L::ClasticSedimentary, 0.55), (L::Carbonate, 0.25), (L::Evaporite, 0.15), (L::FloodBasalt, 0.05)],
            (Some(G::ContinentalRift), _) => &[(L::RiftFill, 0.65), (L::FloodBasalt, 0.20), (L::Granitoid, 0.05), (L::Evaporite, 0.10)],
            (Some(G::ExtendedCrust), true) => &[(L::ClasticSedimentary, 0.35), (L::Metamorphic, 0.30), (L::RiftFill, 0.25), (L::Granitoid, 0.10)],
            (Some(G::CollisionOrogen), true) => &[(L::Metamorphic, 0.38), (L::Granitoid, 0.20), (L::ClasticSedimentary, 0.24), (L::Carbonate, 0.12), (L::Ophiolite, 0.06)],
            (Some(G::OrogenicPlateau), true) => &[(L::ClasticSedimentary, 0.30), (L::Metamorphic, 0.24), (L::Granitoid, 0.20), (L::Carbonate, 0.14), (L::ArcVolcanic, 0.12)],
            (Some(G::ForelandBasin), true) => &[(L::ClasticSedimentary, 0.78), (L::Carbonate, 0.12), (L::Evaporite, 0.10)],
            (Some(G::PaleoOrogen), true) => &[(L::Metamorphic, 0.42), (L::Granitoid, 0.25), (L::ClasticSedimentary, 0.27), (L::Ophiolite, 0.06)],
            // Old basement thrust up through its sedimentary cover, which fills the basins between.
            (Some(G::BasementUplift), true) => &[(L::CratonicBasement, 0.30), (L::Granitoid, 0.20), (L::Metamorphic, 0.15), (L::ClasticSedimentary, 0.30), (L::Carbonate, 0.05)],
            (Some(G::VolcanicArc), true) => &[(L::ArcVolcanic, 0.58), (L::Granitoid, 0.30), (L::ClasticSedimentary, 0.12)],
            (Some(G::VolcanicArc), false) => &[(L::ArcVolcanic, 0.85), (L::Granitoid, 0.10), (L::PelagicSediment, 0.05)],
            (Some(G::AccretionaryWedge), _) => &[(L::AccretionaryMelange, 0.70), (L::Ophiolite, 0.15), (L::ClasticSedimentary, 0.15)],
            (Some(G::ForearcBasin), true) => &[(L::ClasticSedimentary, 0.80), (L::ArcVolcanic, 0.20)],
            (Some(G::BackarcBasin), true) => &[(L::ClasticSedimentary, 0.50), (L::ArcVolcanic, 0.30), (L::RiftFill, 0.20)],
            (Some(G::ContinentalFloodBasalt), true) => &[(L::FloodBasalt, 0.85), (L::ClasticSedimentary, 0.15)],
            (Some(G::ContinentalHotspotTrack), true) => &[(L::FloodBasalt, 0.55), (L::ArcVolcanic, 0.20), (L::Granitoid, 0.10), (L::ClasticSedimentary, 0.15)],
            (Some(G::OceanicPlateau | G::OceanicHotspotTrack | G::HotspotTrack), false) => &[(L::FloodBasalt, 0.80), (L::PelagicSediment, 0.20)],
            (Some(G::MidOceanRidge | G::OceanicFractureZone), false) => &[(L::OceanicBasalt, 1.0)],
            (Some(G::OceanTrench | G::ForearcBasin | G::BackarcBasin), false) => &[(L::PelagicSediment, 0.85), (L::OceanicBasalt, 0.15)],
            // Rigid-lid worlds: battered old crust, lava plains, volcanic constructs.
            (Some(G::CrateredHighlands), _) => &[(L::ImpactBreccia, 0.6), (L::FloodBasalt, 0.25), (L::Anorthosite, 0.15)],
            (Some(G::CrateredIce), _) => &[(L::WaterIce, 0.8), (L::ImpactBreccia, 0.2)],
            (Some(G::VolcanicPlains | G::VolcanicRise), _) => &[(L::FloodBasalt, 0.85), (L::ImpactBreccia, 0.15)],
            (Some(G::ImpactBasin), _) => &[(L::ImpactBreccia, 0.6), (L::FloodBasalt, 0.3), (L::Anorthosite, 0.1)],
            (Some(G::ShieldVolcano), _) => &[(L::FloodBasalt, 1.0)],
            (Some(G::Corona), _) => &[(L::FloodBasalt, 0.7), (L::Metamorphic, 0.3)],
            (Some(G::Tessera), _) => &[(L::Metamorphic, 0.5), (L::FloodBasalt, 0.3), (L::Granitoid, 0.2)],
            // Remaining oceanic settings; young abyssal floor is bare basalt (see `build`).
            (_, false) => &[(L::PelagicSediment, 1.0)],
            // Remaining continental settings on unusual provinces.
            (Some(_), true) => &[(L::ClasticSedimentary, 0.5), (L::Granitoid, 0.2), (L::Metamorphic, 0.3)],
        }
    }
}

/// Keys of extra rocks start here, clear of the built-ins.
const EXTRA_KEY_BASE: u64 = 101;

/// One entry of a setting's mix.
#[derive(Debug, Clone, Copy)]
struct Part {
    key: Key,
    age: AgeRule,
    weight: f64,
}

/// Earth's rocks plus a world's extras, keyed as the lithology partition stores them.
#[derive(Debug, Clone)]
pub struct Catalog {
    rocks: Vec<(Key, RockDef)>,
}

impl Catalog {
    /// # Errors
    /// If a slug repeats, or an occurrence weight is not positive.
    pub fn new(extra: &[RockDef]) -> Result<Self, Error> {
        let mut rocks: Vec<(Key, RockDef)> = Earth::ALL.iter().map(|e| (e.key(), e.def())).collect();
        rocks.extend(extra.iter().enumerate().map(|(i, d)| (Key(EXTRA_KEY_BASE + i as u64), d.clone())));
        let mut slugs = BTreeSet::new();
        for (_, d) in &rocks {
            if !slugs.insert(d.rock.slug.as_str()) {
                return Err(Error::Invalid(format!("rock '{}' is defined twice", d.rock.slug)));
            }
            if d.occurs.iter().any(|o| o.weight.is_nan() || o.weight <= 0.0) {
                return Err(Error::Invalid(format!("rock '{}': occurrence weights must be positive", d.rock.slug)));
            }
        }
        Ok(Self { rocks })
    }

    /// The lithology partition's table.
    #[must_use]
    pub fn table(&self) -> Table<Rock> {
        Table { rows: self.rocks.iter().map(|(k, d)| (*k, d.rock.clone())).collect() }
    }

    /// The mix of a setting: Earth's, plus every extra rock that occurs there.
    fn mix(&self, kind: Option<G>, continental: bool) -> Vec<Part> {
        let mut m: Vec<Part> = Earth::mix(kind, continental).iter().map(|&(e, weight)| Part { key: e.key(), age: e.def().age, weight }).collect();
        let crust = if continental { CrustType::Continental } else { CrustType::Oceanic };
        for (key, d) in &self.rocks {
            let weight: f64 = d.occurs.iter().filter(|o| o.province.is_none_or(|p| Some(p) == kind) && o.crust.is_none_or(|c| c == crust)).map(|o| o.weight).sum();
            if weight > 0.0 {
                m.push(Part { key: *key, age: d.age, weight });
            }
        }
        m
    }
}

fn pick(m: &[Part], u: f64) -> Part {
    let total: f64 = m.iter().map(|x| x.weight).sum();
    let mut acc = 0.0;
    for &p in m {
        acc += p.weight / total;
        if u < acc {
            return p;
        }
    }
    m[m.len() - 1]
}

/// What a rock's age depends on at one cell.
struct Setting {
    kind: Option<G>,
    continental: bool,
    province_age: f64,
    crust_age: f64,
    basement: f64,
}

fn age_ma(rule: AgeRule, s: &Setting, jitter: f64) -> f64 {
    let belt = matches!(s.kind, Some(G::CollisionOrogen | G::OrogenicPlateau | G::PaleoOrogen | G::VolcanicArc | G::ExtendedCrust));
    let p = s.province_age;
    match rule {
        AgeRule::Crust => s.crust_age,
        AgeRule::Basement => s.basement,
        AgeRule::Archean => s.basement.max(2500.0) + 300.0 * jitter,
        AgeRule::Plutonic if belt => p * (1.0 + 0.4 * jitter),
        AgeRule::Plutonic => s.basement * (0.6 + 0.4 * jitter),
        AgeRule::Volcanic => p * jitter,
        AgeRule::Flood if s.continental => p.max(1.0),
        AgeRule::Flood => s.crust_age * (0.3 + 0.6 * jitter),
        AgeRule::Ophiolite => p + 20.0 + 150.0 * jitter,
        AgeRule::Infill => p * (0.2 + 0.8 * jitter),
        AgeRule::Cover => match s.kind {
            Some(G::CollisionOrogen | G::PaleoOrogen) => p + 50.0 + 250.0 * jitter,
            Some(G::ForearcBasin | G::BackarcBasin | G::ContinentalRift | G::ForelandBasin) => p * jitter,
            // Young deposits are the most widespread.
            _ => (2.0 + 538.0 * jitter.powf(1.8)).min(s.basement),
        },
    }
}

/// Lithology and rock age (Ma) per cell.
#[must_use]
pub fn build(grid: &dyn Grid, provinces: &Partition<Province>, crust: &[CrustType], crust_age: Option<&[f32]>, seed: u64, catalog: &Catalog) -> (Partition<Rock>, Vec<f32>) {
    let n = grid.len();
    let pick_noise = TerrainNoise::new(seed);
    let age_noise = TerrainNoise::new(seed ^ 0x5EED_A6E5);
    let raw: Vec<(f64, f64)> = (0..n as Cell)
        .into_par_iter()
        .map(|c| {
            let p = grid.center(c);
            (pick_noise.fbm(p, 7.0, 4) + 0.35 * pick_noise.fbm(p, 23.0, 3), age_noise.fbm(p, 5.0, 3))
        })
        .collect();
    let u = uniform_ranks(&raw.iter().map(|r| r.0).collect::<Vec<_>>());
    let v = uniform_ranks(&raw.iter().map(|r| r.1).collect::<Vec<_>>());
    let mut mixes: BTreeMap<(Option<G>, bool), Vec<Part>> = BTreeMap::new();
    let young_floor = Part { key: Earth::OceanicBasalt.key(), age: AgeRule::Crust, weight: 1.0 };
    let (owner, age): (Vec<u32>, Vec<f32>) = (0..n)
        .map(|i| {
            let continental = crust[i] == CrustType::Continental;
            let prov = provinces.row_at(i as Cell);
            let kind = prov.map(|p| p.kind);
            let crust_age = crust_age.map_or(if continental { 2000.0 } else { 80.0 }, |a| f64::from(a[i]));
            let province_age = prov.map_or(0.0, |p| p.age_ma);
            let basement = if kind == Some(G::Craton) { province_age.max(crust_age) } else { crust_age };
            let s = Setting { kind, continental, province_age, crust_age, basement };
            let mut rock = pick(mixes.entry((kind, continental)).or_insert_with(|| catalog.mix(kind, continental)), u[i]);
            if !continental && matches!(kind, None | Some(G::AbyssalPlain)) && crust_age < 12.0 && rock.key == Earth::PelagicSediment.key() {
                rock = young_floor;
            }
            (rock.key.0 as u32, age_ma(rock.age, &s, v[i]).clamp(0.0, 4400.0) as f32)
        })
        .unzip();
    (Partition { owner: Field::from_values(grid.spec(), owner), table: catalog.table() }, age)
}

/// Erodibility per cell from a lithology partition (0.5 where a cell has no rock).
#[must_use]
pub fn erodibility(lithology: &Partition<Rock>) -> Vec<f64> {
    let by_key: BTreeMap<u32, f64> = lithology.table.rows.iter().map(|(k, r)| (k.0 as u32, r.erodibility)).collect();
    lithology.owner.values.iter().map(|o| by_key.get(o).copied().unwrap_or(0.5)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn earth_mixes_are_normalisable_and_pick_covers_them() {
        let cat = Catalog::new(&[]).unwrap();
        for kind in <G as geoforge_core::surface::Category>::ALL.iter().map(|k| Some(*k)).chain([None]) {
            for cont in [true, false] {
                let m = cat.mix(kind, cont);
                assert!(m.iter().all(|x| x.weight > 0.0));
                assert_eq!(pick(&m, 0.0).key, m[0].key);
                assert_eq!(pick(&m, 0.999_999).key, m[m.len() - 1].key);
            }
        }
        // Built-in keys are stable.
        assert_eq!(cat.table().rows[&Key(1)].slug, "oceanic_basalt");
        assert_eq!(cat.table().rows[&Key(14)].slug, "rift_fill");
    }

    #[test]
    fn extra_rocks_join_their_settings() {
        let star = RockDef {
            rock: Rock { slug: "starmetal".into(), name: "Starmetal".into(), color: [200, 210, 255], erodibility: 0.1, tags: vec!["ore".into()] },
            age: AgeRule::Volcanic,
            occurs: vec![Occurrence { province: Some(G::VolcanicArc), crust: Some(CrustType::Continental), weight: 0.1 }],
        };
        let cat = Catalog::new(std::slice::from_ref(&star)).unwrap();
        let arc = cat.mix(Some(G::VolcanicArc), true);
        assert!(arc.iter().any(|p| p.key == Key(EXTRA_KEY_BASE) && (p.weight - 0.1).abs() < 1e-12 && p.age == AgeRule::Volcanic));
        assert!(!cat.mix(Some(G::VolcanicArc), false).iter().any(|p| p.key == Key(EXTRA_KEY_BASE)));
        assert!(!cat.mix(Some(G::Craton), true).iter().any(|p| p.key == Key(EXTRA_KEY_BASE)));
        assert!(Catalog::new(&[star.clone(), star]).is_err(), "duplicate slug");
    }

    #[test]
    fn extra_rocks_parse_from_toml() {
        let t = r##"
            slug = "dreamstone"
            name = "Dreamstone"
            color = "#9a7fd1"
            erodibility = 0.3
            tags = ["magical"]
            age = "archean"
            occurs = [{ province = "craton", weight = 0.02 }]
        "##;
        let d: RockDef = toml::from_str(t).unwrap();
        assert_eq!(d.rock.color, [0x9a, 0x7f, 0xd1]);
        assert_eq!(d.occurs[0].province, Some(G::Craton));
    }
}
