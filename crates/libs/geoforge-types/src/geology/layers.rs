#![allow(missing_docs)]

//! Stage 2 layer descriptions.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

use super::province::ProvinceId;
use super::{GeologicProvince, ProvinceCharacteristics, TectonicContext};
use crate::tectonics::CrustType;

/// Per-pixel geology raster fields (storage in `geoforge-grid`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum GeologyRaster {
    /// `TerrainMap<ProvinceId>`
    ProvinceId,
    /// `TerrainMap<f32>` — activity / distance weight 0–1.
    Intensity,
}

/// Expected elevation category for a province (Stage 3 hint).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ElevationClass {
    DeepMarine,
    ShallowMarine,
    Lowland,
    Moderate,
    High,
    Extreme,
}

/// Lookup metadata for one province id.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ProvinceInfo {
    /// Raster id.
    pub id: ProvinceId,
    /// Province type.
    pub province_type: GeologicProvince,
    /// Contributing plate ids.
    pub plate_ids: Vec<crate::tectonics::PlateId>,
    /// Dominant crust under this province.
    pub primary_crust_type: CrustType,
    /// Subduction/collision context if applicable.
    pub tectonic_context: Option<TectonicContext>,
    /// Terrain parameters.
    pub characteristics: ProvinceCharacteristics,
    /// Pixel count when rasterized.
    pub pixel_count: usize,
    /// Elevation category hint.
    pub elevation_class: ElevationClass,
}

/// Non-raster Stage 2 output.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, Default)]
pub struct GeologicLayerData {
    /// Province catalog keyed by id.
    pub provinces: HashMap<ProvinceId, ProvinceInfo>,
}

impl GeologicLayerData {
    /// Register a province; returns assigned id.
    pub fn register(&mut self, info: ProvinceInfo) -> ProvinceId {
        let id = info.id;
        self.provinces.insert(id, info);
        id
    }
}

impl GeologicProvince {
    /// Default elevation class for this province.
    #[must_use]
    pub fn default_elevation_class(self) -> ElevationClass {
        use ElevationClass::*;
        use GeologicProvince::*;
        match self {
            CollisionOrogen => Extreme,
            PaleoOrogen => Moderate,
            AccretionaryWedge | ForearcBasin | BackarcBasin => ShallowMarine,
            ContinentalFloodBasalt | Craton | Platform => Lowland,
            OceanicPlateau | HotspotTrack | VolcanicArc => Moderate,
            OceanTrench | AbyssalPlain | FractureZone => DeepMarine,
            ExtendedCrust | ContinentalRift => Lowland,
            MidOceanRidge => ShallowMarine,
        }
    }
}
