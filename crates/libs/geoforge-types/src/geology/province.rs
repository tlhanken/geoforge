//! Geologic province taxonomy and terrain characteristics.

#![allow(missing_docs)]


use serde::{Deserialize, Serialize};

/// Major geologic province types generated from tectonic data.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[repr(u16)]
pub enum GeologicProvince {
    // --- Orogenic ---
    /// Continent–continent collision (Himalayas, Alps).
    CollisionOrogen,
    /// Ancient, eroded collision belt (Appalachians, Urals).
    PaleoOrogen,
    /// Subduction accretion prism.
    AccretionaryWedge,

    // --- Large igneous provinces ---
    ContinentalFloodBasalt,
    OceanicPlateau,
    HotspotTrack,

    // --- Arc systems ---
    VolcanicArc,
    ForearcBasin,
    BackarcBasin,
    OceanTrench,

    // --- Stable continental ---
    Craton,
    Platform,
    ExtendedCrust,

    // --- Extensional ---
    ContinentalRift,

    // --- Oceanic ---
    AbyssalPlain,
    MidOceanRidge,
    FractureZone,
}

/// Province id on the raster (`0` = unset / default basin).
pub type ProvinceId = u16;

/// Parameters that modulate elevation and roughness in Stage 3.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ProvinceCharacteristics {
    /// Elevation multiplier (1.0 = orogen-scale uplift potential).
    pub elevation_intensity: f64,
    /// Terrain roughness 0–1.
    pub roughness: f64,
    /// Typical width (km).
    pub width_km: f64,
    /// Associated tectonic rate (cm/yr); magnitude only.
    pub tectonic_rate_cm_per_year: f64,
}

impl GeologicProvince {
    /// Default characteristics for this province type.
    #[must_use]
    pub fn default_characteristics(self) -> ProvinceCharacteristics {
        use GeologicProvince::*;
        match self {
            CollisionOrogen => ProvinceCharacteristics {
                elevation_intensity: 1.0,
                roughness: 0.9,
                width_km: 1500.0,
                tectonic_rate_cm_per_year: 5.0,
            },
            PaleoOrogen => ProvinceCharacteristics {
                elevation_intensity: 0.4,
                roughness: 0.5,
                width_km: 800.0,
                tectonic_rate_cm_per_year: 0.0,
            },
            AccretionaryWedge => ProvinceCharacteristics {
                elevation_intensity: 0.2,
                roughness: 0.7,
                width_km: 100.0,
                tectonic_rate_cm_per_year: 5.0,
            },
            ContinentalFloodBasalt => ProvinceCharacteristics {
                elevation_intensity: 0.3,
                roughness: 0.3,
                width_km: 500.0,
                tectonic_rate_cm_per_year: 0.0,
            },
            OceanicPlateau => ProvinceCharacteristics {
                elevation_intensity: -0.3,
                roughness: 0.2,
                width_km: 400.0,
                tectonic_rate_cm_per_year: 0.0,
            },
            HotspotTrack => ProvinceCharacteristics {
                elevation_intensity: 0.6,
                roughness: 0.4,
                width_km: 200.0,
                tectonic_rate_cm_per_year: 0.0,
            },
            VolcanicArc => ProvinceCharacteristics {
                elevation_intensity: 0.8,
                roughness: 0.7,
                width_km: 150.0,
                tectonic_rate_cm_per_year: 5.0,
            },
            ForearcBasin => ProvinceCharacteristics {
                elevation_intensity: -0.2,
                roughness: 0.2,
                width_km: 80.0,
                tectonic_rate_cm_per_year: 5.0,
            },
            BackarcBasin => ProvinceCharacteristics {
                elevation_intensity: -0.4,
                roughness: 0.2,
                width_km: 200.0,
                tectonic_rate_cm_per_year: 3.0,
            },
            OceanTrench => ProvinceCharacteristics {
                elevation_intensity: -1.2,
                roughness: 0.3,
                width_km: 75.0,
                tectonic_rate_cm_per_year: 8.0,
            },
            Craton => ProvinceCharacteristics {
                elevation_intensity: 0.2,
                roughness: 0.1,
                width_km: 2000.0,
                tectonic_rate_cm_per_year: 0.0,
            },
            Platform => ProvinceCharacteristics {
                elevation_intensity: 0.15,
                roughness: 0.1,
                width_km: 1500.0,
                tectonic_rate_cm_per_year: 0.0,
            },
            ExtendedCrust => ProvinceCharacteristics {
                elevation_intensity: 0.1,
                roughness: 0.2,
                width_km: 300.0,
                tectonic_rate_cm_per_year: 0.0,
            },
            ContinentalRift => ProvinceCharacteristics {
                elevation_intensity: 0.3,
                roughness: 0.5,
                width_km: 100.0,
                tectonic_rate_cm_per_year: 2.0,
            },
            AbyssalPlain => ProvinceCharacteristics {
                elevation_intensity: -1.0,
                roughness: 0.05,
                width_km: 5000.0,
                tectonic_rate_cm_per_year: 0.0,
            },
            MidOceanRidge => ProvinceCharacteristics {
                elevation_intensity: -0.5,
                roughness: 0.3,
                width_km: 120.0,
                tectonic_rate_cm_per_year: -5.0,
            },
            FractureZone => ProvinceCharacteristics {
                elevation_intensity: -0.8,
                roughness: 0.4,
                width_km: 50.0,
                tectonic_rate_cm_per_year: 0.0,
            },
        }
    }

    /// All variants for testing and validation.
    #[must_use]
    pub fn all() -> &'static [Self] {
        use GeologicProvince::*;
        &[
            CollisionOrogen,
            PaleoOrogen,
            AccretionaryWedge,
            ContinentalFloodBasalt,
            OceanicPlateau,
            HotspotTrack,
            VolcanicArc,
            ForearcBasin,
            BackarcBasin,
            OceanTrench,
            Craton,
            Platform,
            ExtendedCrust,
            ContinentalRift,
            AbyssalPlain,
            MidOceanRidge,
            FractureZone,
        ]
    }
}
