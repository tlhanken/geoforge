//! The seas: ice, reefs and how much life they hold, worked out from the stored climate,
//! currents and seabed (derived views, never stored).
//!
//! - **Sea ice**: seasonal where the coldest month's air is under −4 °C, perennial where the
//!   warmest month's stays under −6 °C. Sea water freezes at −1.8 °C, but the modelled polar
//!   seasons are muted (the Arctic's warmest month comes out −8 °C, observed ~0), so the
//!   thresholds are set on the ETOPO-forced Earth run to NSIDC's extents.
//! - **Coral reefs**: shallow sea (above −60 m, or along a shore, where fringing reefs grow
//!   on cells too coarse to show the shallows) whose coldest month is at least 18 °C and
//!   warmest at most 31 °C, without a cold upwelling or a big river's muddy plume (150 km
//!   across at 2000 m³/s, as the square root of the flow beyond). **Atolls** on hotspot tracks in
//!   such water more than 300 km behind the active end: islands that have sunk while their
//!   reefs grew (Darwin 1842).
//! - **Productivity**, 0–1: plankton where nutrients reach sunlit water: coastal and
//!   equatorial upwelling (a cold surface anomaly), high-latitude seas mixed by winter storms,
//!   shallow shelves and river plumes; the subtropical gyres are deserts.

use geoforge_core::grid::{cells_within, Cell, Grid};
use geoforge_core::product::Products;
use geoforge_core::surface::GeologicProvince as P;

/// Sea ice over a year.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum SeaIce {
    Land,
    Open,
    Seasonal,
    Perennial,
}

impl SeaIce {
    pub const ALL: [Self; 4] = [Self::Land, Self::Open, Self::Seasonal, Self::Perennial];

    #[must_use]
    pub fn info(self) -> (&'static str, &'static str, [u8; 3]) {
        match self {
            Self::Land => ("land", "Land", [205, 200, 190]),
            Self::Open => ("open", "Open water all year", [60, 110, 170]),
            Self::Seasonal => ("seasonal", "Ice in winter", [160, 200, 230]),
            Self::Perennial => ("perennial", "Ice all year", [250, 250, 255]),
        }
    }
}

/// Coral reefs and atolls.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Reef {
    Land,
    None,
    Reef,
    Atoll,
}

impl Reef {
    pub const ALL: [Self; 4] = [Self::Land, Self::None, Self::Reef, Self::Atoll];

    #[must_use]
    pub fn info(self) -> (&'static str, &'static str, [u8; 3]) {
        match self {
            Self::Land => ("land", "Land", [205, 200, 190]),
            Self::None => ("none", "No reef", [60, 110, 170]),
            Self::Reef => ("reef", "Coral reef", [255, 120, 110]),
            Self::Atoll => ("atoll", "Atoll", [255, 200, 60]),
        }
    }
}

/// Air temperature of the coldest month below which sea ice forms in winter, °C.
const ICE_WINTER_C: f64 = -4.0;
/// Air temperature of the warmest month below which the ice lasts all year, °C.
const ICE_SUMMER_C: f64 = -6.0;
/// River flow whose plume keeps reefs away, m³/s.
const MUDDY_FLOW_M3_S: f64 = 2000.0;
/// Deepest sea a reef grows in, m (sampled over a cell, so a little deeper than the 50 m
/// reefs build to).
const REEF_DEPTH_M: f64 = -60.0;

/// What the sea views read.
pub struct Sea<'a> {
    pub grid: &'a dyn Grid,
    pub radius_km: f64,
    pub elevation: &'a [f32],
    pub temperature: &'a [f32],
    pub range: &'a [f32],
    /// The sea surface's departure from its latitude's mean, °C (0 where not modelled).
    pub sst_anomaly: Option<&'a [f32]>,
    pub products: &'a Products,
    /// Below sea level on a world with seas.
    pub sea: Vec<bool>,
}

impl<'a> Sea<'a> {
    /// `None` until the climate has run.
    #[must_use]
    pub fn read(p: &'a Products, grid: &'a dyn Grid, radius_km: f64, sea_world: bool) -> Option<Self> {
        let elevation = &p.elevation().or(p.rough_elevation())?.values;
        Some(Self {
            grid,
            radius_km,
            elevation,
            temperature: &p.temperature()?.values,
            range: &p.temperature_range()?.values,
            sst_anomaly: p.sst_anomaly().map(|f| f.values.as_slice()),
            products: p,
            sea: elevation.iter().map(|&e| sea_world && e < 0.0).collect(),
        })
    }

    fn coldest(&self, c: usize) -> f64 {
        f64::from(self.temperature[c]) - 0.5 * f64::from(self.range[c])
    }

    fn warmest(&self, c: usize) -> f64 {
        f64::from(self.temperature[c]) + 0.5 * f64::from(self.range[c])
    }

    fn anomaly(&self, c: usize) -> f64 {
        self.sst_anomaly.map_or(0.0, |a| f64::from(a[c]))
    }

    /// Sea ice at every cell.
    #[must_use]
    pub fn ice(&self) -> Vec<SeaIce> {
        (0..self.sea.len())
            .map(|c| match () {
                () if !self.sea[c] => SeaIce::Land,
                () if self.warmest(c) < ICE_SUMMER_C => SeaIce::Perennial,
                () if self.coldest(c) < ICE_WINTER_C => SeaIce::Seasonal,
                () => SeaIce::Open,
            })
            .collect()
    }

    /// Reefs and atolls at every cell.
    #[must_use]
    pub fn reefs(&self) -> Vec<Reef> {
        let n = self.sea.len();
        // Muddy, fresh water off big rivers' mouths: a plume 150 km across at 2000 m³/s,
        // growing as the square root of the flow (the Amazon's reaches ~1000 km).
        let mut muddy = vec![false; n];
        if let Some(r) = self.products.rivers() {
            for node in r.nodes.rows.values().filter(|nd| nd.discharge_m3_s >= MUDDY_FLOW_M3_S && self.near_sea(nd.cell as usize)) {
                let reach_km = (150.0 * (node.discharge_m3_s / MUDDY_FLOW_M3_S).sqrt()).min(1000.0);
                for c in cells_within(self.grid, &[node.cell], reach_km / self.radius_km, |_| true) {
                    muddy[c as usize] = true;
                }
            }
        }
        let provinces = self.products.provinces();
        // Old, sunken parts of a hotspot track: more than 300 km behind its active end.
        let atoll_track = |c: usize| {
            provinces.and_then(|p| p.row_at(c as Cell)).is_some_and(|p| {
                matches!(p.kind, P::OceanicHotspotTrack | P::HotspotTrack) && p.anchor.is_none_or(|a| self.grid.center(c as Cell).angle_to(a) * self.radius_km > 300.0)
            })
        };
        (0..n)
            .map(|c| {
                if !self.sea[c] {
                    return Reef::Land;
                }
                let warm_water = self.coldest(c) >= 18.0 && self.warmest(c) <= 31.0 && self.anomaly(c) > -1.5 && !muddy[c];
                if !warm_water {
                    Reef::None
                } else if atoll_track(c) {
                    Reef::Atoll
                } else if f64::from(self.elevation[c]) >= REEF_DEPTH_M || self.touches_land(c) {
                    Reef::Reef
                } else {
                    Reef::None
                }
            })
            .collect()
    }

    fn touches_land(&self, c: usize) -> bool {
        self.grid.neighbors(c as Cell).iter().any(|n| !self.sea[n as usize])
    }

    fn near_sea(&self, c: usize) -> bool {
        self.sea[c] || self.grid.neighbors(c as Cell).iter().any(|n| self.sea[n as usize])
    }

    /// Productivity at every sea cell, 0–1; `None` on land and under perennial ice.
    #[must_use]
    pub fn productivity(&self) -> Vec<Option<f64>> {
        let ice = self.ice();
        let n = self.sea.len();
        let coast: Vec<Cell> = (0..n).filter(|&c| !self.sea[c]).map(|c| c as Cell).collect();
        let shore_km: Vec<f64> = geoforge_core::grid::distance_field(self.grid, coast, 400.0 / self.radius_km, |c| self.sea[c as usize]).into_iter().map(|d| d * self.radius_km).collect();
        (0..n)
            .map(|c| {
                if !self.sea[c] || ice[c] == SeaIce::Perennial {
                    return None;
                }
                let lat = self.grid.center(c as Cell).to_lat_lon().lat_deg.abs();
                // Upwelling: each degree of cold anomaly brings up nutrients.
                let upwelling = (-self.anomaly(c) / 3.0).clamp(0.0, 1.0);
                // Winter storms mix deep water up poleward of ~40°.
                let mixing = ((lat - 35.0) / 25.0).clamp(0.0, 1.0) * 0.6;
                // Shelves and coasts: shallow water and runoff.
                let shelf = if f64::from(self.elevation[c]) > -200.0 { 0.35 } else { 0.0 };
                let coastal = if shore_km[c].is_finite() { 0.25 * (1.0 - shore_km[c] / 400.0) } else { 0.0 };
                Some((0.05 + upwelling + mixing + shelf + coastal).min(1.0))
            })
            .collect()
    }

    /// Area of `ice` classes in each hemisphere, km²: (north, south).
    #[must_use]
    pub fn ice_area_km2(&self, ice: &[SeaIce], classes: &[SeaIce]) -> (f64, f64) {
        let mut out = (0.0, 0.0);
        for (c, i) in ice.iter().enumerate() {
            if classes.contains(i) {
                let a = self.grid.area_sr(c as Cell) * self.radius_km * self.radius_km;
                if self.grid.center(c as Cell).to_lat_lon().lat_deg >= 0.0 {
                    out.0 += a;
                } else {
                    out.1 += a;
                }
            }
        }
        out
    }
}
