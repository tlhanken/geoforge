//! Stage 1 layer *descriptions* (raster data lives in `geoforge-grid`).

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

use super::boundaries::{BoundarySegmentMeta, BoundaryStatistics};
use super::{PlateId, PlateSeed, PlateStats, PlateType};
use crate::tectonics::CrustType;

/// Per-pixel tectonic layer field identifiers.
///
/// Actual `TerrainMap<T>` storage is provided by downstream crates; this enum
/// documents which typed rasters Stage 1 produces.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum TectonicRaster {
    /// `TerrainMap<PlateId>`
    PlateId,
    /// `TerrainMap<CrustType>`
    CrustType,
    /// `TerrainMap<f32>` — crust age in Ma.
    CrustAgeMa,
    /// `TerrainMap<BoundaryType>`
    BoundaryType,
    /// `TerrainMap<f32>` — convergence rate cm/yr (negative = divergent).
    ConvergenceRateCmPerYear,
}

/// Non-raster metadata produced by Stage 1.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TectonicLayerData {
    /// Plate seed points.
    pub plate_seeds: Vec<PlateSeed>,
    /// Per-plate statistics keyed by id.
    pub plate_stats: HashMap<PlateId, PlateStats>,
    /// Boundary segments.
    pub boundaries: Vec<BoundarySegmentMeta>,
    /// Aggregate boundary stats.
    pub boundary_statistics: BoundaryStatistics,
}

impl TectonicLayerData {
    /// Assign [`PlateType`] from size percentiles for all plates.
    pub fn assign_plate_types_from_size(&mut self) {
        let total = self.plate_stats.len();
        if total == 0 {
            return;
        }

        let mut by_area: Vec<(PlateId, f64)> = self
            .plate_stats
            .iter()
            .map(|(&id, s)| (id, s.area_km2))
            .collect();
        by_area.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

        for (rank, (id, _)) in by_area.iter().enumerate() {
            let percentile = rank as f64 / total as f64;
            if let Some(stats) = self.plate_stats.get_mut(id) {
                stats.plate_type = PlateType::from_size_percentile(percentile);
            }
        }
    }
}

/// Typical oceanic crust age range (Ma) for initialization.
pub const OCEANIC_CRUST_AGE_MA_RANGE: (f32, f32) = (0.0, 180.0);

/// Typical continental crust age range (Ma) for initialization.
pub const CONTINENTAL_CRUST_AGE_MA_RANGE: (f32, f32) = (1000.0, 3000.0);

/// Pick a default crust age (Ma) from crust type (mid-range).
#[must_use]
pub fn default_crust_age_ma(crust: CrustType) -> f32 {
    match crust {
        CrustType::Oceanic => {
            (OCEANIC_CRUST_AGE_MA_RANGE.0 + OCEANIC_CRUST_AGE_MA_RANGE.1) / 2.0
        }
        CrustType::Continental => {
            (CONTINENTAL_CRUST_AGE_MA_RANGE.0 + CONTINENTAL_CRUST_AGE_MA_RANGE.1) / 2.0
        }
        CrustType::Transitional => 500.0,
    }
}
