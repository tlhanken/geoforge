//! Integration-style tests for cross-module invariants.

use crate::{
    geology::{GeologicProvince, GeologicLayerData, ProvinceInfo},
    pipeline::PipelineLayerId,
    planetary::PlanetaryParams,
    seed::Seed,
    tectonics::{CrustType, PlateType, TectonicLayerData, default_crust_age_ma},
};

#[test]
fn plate_type_default_crust_matches_character() {
    assert_eq!(
        PlateType::Oceanic.default_crust_type(),
        CrustType::Oceanic
    );
    assert_eq!(
        PlateType::Continental.default_crust_type(),
        CrustType::Continental
    );
}

#[test]
fn all_provinces_have_characteristics() {
    for &p in GeologicProvince::all() {
        let c = p.default_characteristics();
        assert!(c.width_km > 0.0);
        let _ = p.default_elevation_class();
    }
}

#[test]
fn serde_roundtrip_planetary_params() {
    let earth = PlanetaryParams::earth();
    let json = serde_json::to_string(&earth).unwrap();
    let back: PlanetaryParams = serde_json::from_str(&json).unwrap();
    assert_eq!(earth, back);
}

#[test]
fn pipeline_layers_parse_roundtrip() {
    for &layer in PipelineLayerId::pipeline_order() {
        let s = layer.as_str();
        assert_eq!(PipelineLayerId::parse_str(s), Some(layer));
    }
}

#[test]
fn seed_hierarchy_is_stable() {
    let root = Seed::new(1);
    let t = root.child(b"tectonics");
    let g = root.child(b"geology");
    assert_ne!(t, g);
    assert_eq!(t.child(b"ridges"), t.child(b"ridges"));
}

#[test]
fn tectonic_metadata_assigns_plate_types() {
    use crate::tectonics::{PlateId, PlateMotion, PlateSeed, PlateStats};
    use crate::coordinates::{LatLon, PixelCoord};

    let seed = PlateSeed {
        id: 1,
        pixel: PixelCoord { x: 0, y: 0 },
        latlon: LatLon::new(0.0, 0.0),
        motion: PlateMotion {
            direction_deg: 0.0,
            speed_cm_per_year: 5.0,
        },
    };
    let mut data = TectonicLayerData {
        plate_seeds: vec![seed.clone()],
        plate_stats: [
            (
                1_u16 as PlateId,
                PlateStats {
                    pixels: 100,
                    percentage: 10.0,
                    area_km2: 1e7,
                    seed: seed.clone(),
                    plate_type: PlateType::Oceanic,
                },
            ),
            (
                2,
                PlateStats {
                    pixels: 900,
                    percentage: 90.0,
                    area_km2: 1e8,
                    seed,
                    plate_type: PlateType::Oceanic,
                },
            ),
        ]
        .into_iter()
        .collect(),
        boundaries: vec![],
        boundary_statistics: Default::default(),
    };
    data.assign_plate_types_from_size();
    assert_eq!(data.plate_stats[&2].plate_type, PlateType::Continental);
    assert_eq!(data.plate_stats[&1].plate_type, PlateType::Oceanic);
}

#[test]
fn geology_register_roundtrip() {
    use crate::geology::ElevationClass;
    use crate::tectonics::CrustType;

    let mut data = GeologicLayerData::default();
    let info = ProvinceInfo {
        id: 1,
        province_type: GeologicProvince::Craton,
        plate_ids: vec![1],
        primary_crust_type: CrustType::Continental,
        tectonic_context: None,
        characteristics: GeologicProvince::Craton.default_characteristics(),
        pixel_count: 0,
        elevation_class: ElevationClass::Lowland,
    };
    assert_eq!(data.register(info), 1);
    assert!(data.provinces.contains_key(&1));
}

#[test]
fn default_crust_ages_are_positive() {
    for crust in [CrustType::Oceanic, CrustType::Continental, CrustType::Transitional] {
        assert!(default_crust_age_ma(crust) > 0.0);
    }
}
