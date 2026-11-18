//! Entity and category types of the surface products (data only; algorithms live in
//! `geoforge-surface`).

use crate::geom::UnitVec3;
use crate::grid::Cell;
use crate::key::Key;
use serde::{Deserialize, Serialize};

/// Categorical cell values stored as `u8` in fields.
pub trait Category: Copy + Sized + 'static {
    const ALL: &'static [Self];
    fn to_u8(self) -> u8;
    fn from_u8(v: u8) -> Option<Self>;
    fn name(self) -> &'static str;
    fn color(self) -> [u8; 3];
}

macro_rules! category {
    ($ty:ident { $($var:ident = $n:literal, $label:literal, [$r:literal, $g:literal, $b:literal]),+ $(,)? }) => {
        #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
        #[serde(rename_all = "snake_case")]
        #[repr(u8)]
        pub enum $ty { $($var = $n),+ }
        impl Category for $ty {
            const ALL: &'static [Self] = &[$($ty::$var),+];
            fn to_u8(self) -> u8 { self as u8 }
            fn from_u8(v: u8) -> Option<Self> { match v { $($n => Some($ty::$var),)+ _ => None } }
            fn name(self) -> &'static str { match self { $($ty::$var => $label),+ } }
            fn color(self) -> [u8; 3] { match self { $($ty::$var => [$r, $g, $b]),+ } }
        }
    };
}

category!(CrustType {
    Oceanic = 0, "Oceanic", [52, 101, 164],
    Continental = 1, "Continental", [196, 160, 0],
});

/// A tectonic plate. Its motion is a rigid rotation about `euler_pole`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Plate {
    pub name: String,
    pub euler_pole: UnitVec3,
    /// Rotation rate, degrees per million years (right-hand rule about the pole).
    pub rate_deg_per_myr: f64,
    /// Where the plate was seeded (its "home" point).
    pub seed_point: UnitVec3,
}

impl Plate {
    /// Surface velocity at `p`, km per million years (= mm/yr), as a tangent vector.
    #[must_use]
    pub fn velocity_km_per_myr(&self, p: UnitVec3, radius_km: f64) -> crate::geom::Vec3 {
        let omega = self.euler_pole.v() * self.rate_deg_per_myr.to_radians();
        omega.cross(p.v()) * radius_km
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum BoundaryKind {
    Convergent,
    Divergent,
    Transform,
}

impl BoundaryKind {
    #[must_use]
    pub fn color(self) -> [u8; 3] {
        match self {
            BoundaryKind::Convergent => [220, 40, 40],
            BoundaryKind::Divergent => [40, 110, 230],
            BoundaryKind::Transform => [40, 200, 70],
        }
    }
}

/// A stretch of boundary between two plates with one kind of motion. Derived from the
/// plates partition and the plates' Euler poles; never stored as truth.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BoundarySegment {
    pub a: Key,
    pub b: Key,
    pub kind: BoundaryKind,
    /// Cells of both plates that touch the other plate along this stretch, ascending.
    pub cells: Vec<Cell>,
    /// Mean normal convergence (+) or divergence rate for convergent/divergent, mean
    /// relative speed for transform; cm/yr.
    pub rate_cm_yr: f64,
    pub length_km: f64,
}

category!(GeologicProvince {
    CollisionOrogen = 0, "Collision Orogen", [139, 69, 19],
    PaleoOrogen = 1, "Paleo-Orogen", [188, 143, 143],
    AccretionaryWedge = 2, "Accretionary Wedge", [222, 184, 135],
    ContinentalFloodBasalt = 3, "Continental Flood Basalt", [128, 0, 128],
    OceanicPlateau = 4, "Oceanic Plateau", [148, 0, 211],
    HotspotTrack = 5, "Hotspot Track", [255, 0, 255],
    VolcanicArc = 6, "Volcanic Arc", [255, 69, 0],
    ForearcBasin = 7, "Forearc Basin", [176, 196, 222],
    BackarcBasin = 8, "Backarc Basin", [119, 136, 153],
    Craton = 9, "Shield (Craton)", [255, 140, 0],
    Platform = 10, "Platform", [255, 228, 196],
    IntracratonicBasin = 11, "Intracratonic Basin", [210, 180, 140],
    ContinentalRift = 12, "Continental Rift", [255, 215, 0],
    ExtendedCrust = 13, "Extended Crust (Passive Margin)", [240, 230, 140],
    MidOceanRidge = 14, "Mid-Ocean Ridge", [0, 255, 255],
    AbyssalPlain = 15, "Abyssal Plain", [65, 105, 225],
    OceanTrench = 16, "Ocean Trench", [0, 0, 139],
    OceanicFractureZone = 17, "Oceanic Fracture Zone", [72, 209, 204],
    OceanicHotspotTrack = 18, "Oceanic Hotspot Track", [138, 43, 226],
    ContinentalHotspotTrack = 19, "Continental Hotspot Track", [153, 50, 204],
});

/// A geologic province: one painted region of the provinces partition.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Province {
    pub kind: GeologicProvince,
    /// Relative contribution to elevation (positive up, negative down), 0..1 scale.
    pub elevation_intensity: f64,
    pub roughness: f64,
    pub width_km: f64,
    pub intensity: f64,
    pub convergence_rate_cm_yr: f64,
    /// Index of the boundary segment this province was built from, if any.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub boundary: Option<usize>,
    /// Age of the province's formation, Ma (for the timeline and elevation decay).
    #[serde(default)]
    pub age_ma: f64,
    /// A defining point: a hotspot's active end, a LIP's or shield's centre.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub anchor: Option<UnitVec3>,
}

category!(Biome {
    Ocean = 0, "Ocean", [28, 62, 120],
    Lake = 1, "Lake", [60, 110, 190],
    IceSheet = 2, "Ice Sheet", [240, 248, 255],
    Tundra = 3, "Tundra", [176, 184, 160],
    BorealForest = 4, "Boreal Forest", [58, 98, 72],
    TemperateRainforest = 5, "Temperate Rainforest", [38, 120, 80],
    TemperateForest = 6, "Temperate Deciduous Forest", [92, 148, 70],
    TemperateGrassland = 7, "Temperate Grassland", [186, 196, 110],
    Mediterranean = 8, "Mediterranean Woodland", [160, 150, 70],
    ColdDesert = 9, "Cold Desert", [196, 182, 150],
    HotDesert = 10, "Hot Desert", [234, 208, 140],
    Savanna = 11, "Savanna", [200, 180, 80],
    TropicalSeasonalForest = 12, "Tropical Seasonal Forest", [120, 160, 50],
    TropicalRainforest = 13, "Tropical Rainforest", [20, 110, 30],
    Alpine = 14, "Alpine", [150, 140, 130],
    Wetland = 15, "Wetland", [70, 130, 120],
});

/// A river node: a cell the network passes through.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RiverNode {
    pub cell: Cell,
    /// Mean discharge, m³/s.
    pub discharge_m3_s: f64,
}

/// A reach between two river nodes, flowing `from` → `to`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RiverReach {
    pub strahler_order: u8,
    pub discharge_m3_s: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Lake {
    pub surface_elevation_m: f64,
    /// No outflow to the sea (salt lake / closed basin).
    pub endorheic: bool,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Watershed {
    /// Where the basin drains: an ocean cell, or the lowest cell of a closed basin.
    pub outlet: Cell,
    pub reaches_sea: bool,
}
