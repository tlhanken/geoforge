//! Surface lithology and rock age from geologic provinces.
//!
//! Each province kind has a characteristic mix of rock (an orogen exposes metamorphic cores,
//! granite plutons and deformed sediments; a platform is sediment over basement; an arc is
//! volcanic with plutonic roots). A smooth noise field picks within that mix, so rock types
//! form map-scale patches rather than speckle, and the same noise value means the same rock
//! across a province. Ages follow the province's history: basement as old as the craton,
//! cover sediments younger, arc volcanics youngest.
//!
//! Lithology feeds erosion (hard basement resists, soft sediment wears down) and, later,
//! resources (porphyry copper in arcs, gold and iron in greenstone, salt and hydrocarbons in
//! basins).

use crate::noise::TerrainNoise;
use geoforge_core::grid::{Cell, Grid};
use geoforge_core::shapes::Partition;
use geoforge_core::surface::{CrustType, GeologicProvince as G, Lithology as L, Province};
use rayon::prelude::*;

type Mix = &'static [(L, f64)];

/// The rock mix of a province (continental or oceanic setting).
#[must_use]
pub fn mix(kind: Option<G>, continental: bool) -> Mix {
    match (kind, continental) {
        (Some(G::Craton), true) => &[(L::CratonicBasement, 0.62), (L::Greenstone, 0.16), (L::Granitoid, 0.12), (L::ClasticSedimentary, 0.10)],
        (Some(G::Platform), true) | (None, true) => {
            &[(L::ClasticSedimentary, 0.48), (L::Carbonate, 0.32), (L::Evaporite, 0.04), (L::CratonicBasement, 0.10), (L::Granitoid, 0.06)]
        }
        (Some(G::IntracratonicBasin), true) => &[(L::ClasticSedimentary, 0.55), (L::Carbonate, 0.25), (L::Evaporite, 0.15), (L::FloodBasalt, 0.05)],
        (Some(G::ContinentalRift), true) => &[(L::RiftFill, 0.65), (L::FloodBasalt, 0.20), (L::Granitoid, 0.05), (L::Evaporite, 0.10)],
        (Some(G::ExtendedCrust), true) => &[(L::ClasticSedimentary, 0.35), (L::Metamorphic, 0.30), (L::RiftFill, 0.25), (L::Granitoid, 0.10)],
        (Some(G::CollisionOrogen), true) => {
            &[(L::Metamorphic, 0.38), (L::Granitoid, 0.20), (L::ClasticSedimentary, 0.24), (L::Carbonate, 0.12), (L::Ophiolite, 0.06)]
        }
        (Some(G::PaleoOrogen), true) => &[(L::Metamorphic, 0.42), (L::Granitoid, 0.25), (L::ClasticSedimentary, 0.27), (L::Ophiolite, 0.06)],
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
        // Remaining oceanic settings; abyssal plains get basalt or sediment by crust age below.
        (_, false) => &[(L::PelagicSediment, 1.0)],
        // Remaining continental settings on unusual provinces.
        (Some(_), true) => &[(L::ClasticSedimentary, 0.5), (L::Granitoid, 0.2), (L::Metamorphic, 0.3)],
    }
}

fn pick(m: Mix, u: f64) -> L {
    let total: f64 = m.iter().map(|x| x.1).sum();
    let mut acc = 0.0;
    for &(l, w) in m {
        acc += w / total;
        if u < acc {
            return l;
        }
    }
    m[m.len() - 1].0
}

/// Values mapped to their rank in [0, 1): a uniform field with the noise's spatial structure.
fn uniform(v: &[f64]) -> Vec<f64> {
    let mut idx: Vec<usize> = (0..v.len()).collect();
    idx.sort_by(|a, b| v[*a].total_cmp(&v[*b]).then(a.cmp(b)));
    let mut out = vec![0.0; v.len()];
    let n = v.len() as f64;
    for (r, i) in idx.into_iter().enumerate() {
        out[i] = (r as f64 + 0.5) / n;
    }
    out
}

/// Lithology and rock age (Ma) per cell.
#[must_use]
pub fn build(grid: &dyn Grid, provinces: &Partition<Province>, crust: &[CrustType], crust_age: Option<&[f32]>, seed: u64) -> (Vec<L>, Vec<f32>) {
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
    let u = uniform(&raw.iter().map(|r| r.0).collect::<Vec<_>>());
    let v = uniform(&raw.iter().map(|r| r.1).collect::<Vec<_>>());
    (0..n)
        .map(|i| {
            let cont = crust[i] == CrustType::Continental;
            let prov = provinces.owner_of(i as Cell).and_then(|k| provinces.table.get(k));
            let kind = prov.map(|p| p.kind);
            let p_age = prov.map_or(0.0, |p| p.age_ma);
            let c_age = crust_age.map_or(if cont { 2000.0 } else { 80.0 }, |a| f64::from(a[i]));
            let mut lith = pick(mix(kind, cont), u[i]);
            if !cont && matches!(kind, None | Some(G::AbyssalPlain)) && c_age < 12.0 {
                lith = L::OceanicBasalt;
            }
            let basement = if kind == Some(G::Craton) { p_age.max(c_age) } else { c_age };
            let jitter = v[i];
            let age = match lith {
                L::OceanicBasalt | L::PelagicSediment => c_age,
                L::CratonicBasement => basement,
                L::Greenstone => basement.max(2500.0) + 300.0 * jitter,
                L::Granitoid | L::Metamorphic => {
                    if matches!(kind, Some(G::CollisionOrogen | G::PaleoOrogen | G::VolcanicArc | G::ExtendedCrust)) {
                        p_age * (1.0 + 0.4 * jitter)
                    } else {
                        basement * (0.6 + 0.4 * jitter)
                    }
                }
                L::ArcVolcanic => p_age * jitter,
                L::FloodBasalt => {
                    if cont {
                        p_age.max(1.0)
                    } else {
                        c_age * (0.3 + 0.6 * jitter)
                    }
                }
                L::Ophiolite => p_age + 20.0 + 150.0 * jitter,
                L::AccretionaryMelange | L::RiftFill => p_age * (0.2 + 0.8 * jitter),
                L::ClasticSedimentary | L::Carbonate | L::Evaporite => {
                    if matches!(kind, Some(G::CollisionOrogen | G::PaleoOrogen)) {
                        p_age + 50.0 + 250.0 * jitter
                    } else if matches!(kind, Some(G::ForearcBasin | G::BackarcBasin | G::ContinentalRift)) {
                        p_age * jitter
                    } else {
                        // Cover sediments: young deposits are the most widespread.
                        (2.0 + 538.0 * jitter.powf(1.8)).min(basement)
                    }
                }
            };
            (lith, age.clamp(0.0, 4400.0) as f32)
        })
        .unzip()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mixes_are_normalisable_and_pick_covers_them() {
        for kind in <G as geoforge_core::surface::Category>::ALL.iter().map(|k| Some(*k)).chain([None]) {
            for cont in [true, false] {
                let m = mix(kind, cont);
                assert!(m.iter().all(|x| x.1 > 0.0));
                assert_eq!(pick(m, 0.0), m[0].0);
                assert_eq!(pick(m, 0.999_999), m[m.len() - 1].0);
            }
        }
    }
}
