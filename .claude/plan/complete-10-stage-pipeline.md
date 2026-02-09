# Complete 10-Stage Pipeline Architecture

## Design Principles

1. **Strictly Linear**: No stage modifies previous stages' data (enables manual injection)
2. **Per-Pixel Architecture**: All data as `TerrainMap<T>` for queryability
3. **Physics-Based**: Climate driven by planetary parameters, not arbitrary rules
4. **Cached Computations**: Slope computed once from elevation, stored for performance
5. **Geological Realism**: Deposit types, erosion modeling, age-based properties

---

## Stage 0: Planetary Parameters

**Output:** Single `PlanetaryParams` struct (not per-pixel)

```rust
pub struct PlanetaryParams {
    // Physical
    pub radius_km: f64,              // Default: 6371 (Earth)
    pub mass_earth: f64,             // Default: 1.0
    pub gravity_ms2: f64,            // Default: 9.81
    pub rotation_period_hours: f64,  // Default: 24.0
    pub axial_tilt_degrees: f64,     // Default: 23.5

    // Stellar/Orbital
    pub insolation_wm2: f64,         // Default: 1361 (Earth)
    pub orbital_eccentricity: f64,   // Default: 0.017
    pub orbital_period_days: f64,    // Default: 365.25

    // Internal Heat
    pub core_temperature_k: f64,     // Default: 5200
    pub heat_flow_mwm2: f64,         // Default: 87
    pub magnetic_field_ut: f64,      // Default: 50 (microtesla)

    // Composition
    pub ocean_coverage: f64,         // Default: 0.71
    pub atmosphere_pressure_kpa: f64,// Default: 101.325
    pub greenhouse_warming_k: f64,   // Default: 33 (atmosphere effect)

    // Fantasy
    pub fantasy_params: Option<FantasyParams>,
}
```

**Presets:** Earth, Mars (0.43 insolation, thin atmo), Venus (2600 W/m², thick atmo), SuperEarth (1.5 mass, 1.2 radius)

---

## Stage 1: Tectonic Foundation

**Input:** PlanetaryParams
**Output:** `TectonicLayers`

```rust
pub struct TectonicLayers {
    pub plate_id: TerrainMap<u16>,
    pub crust_type: TerrainMap<CrustType>,           // Per-pixel!
    pub crust_age_ma: TerrainMap<f32>,               // Million years
    pub lithosphere_thickness_km: TerrainMap<f32>,
    pub boundary_type: TerrainMap<BoundaryType>,     // Convergent/Divergent/Transform/None
    pub convergence_rate_cmyr: TerrainMap<f32>,      // Negative = divergent
}

pub enum CrustType {
    Oceanic,        // Mafic, dense (3.0 g/cm³), thin (7 km)
    Continental,    // Felsic, light (2.7 g/cm³), thick (35 km)
    Transitional,   // Island arcs, obducted terranes
}
```

**Key Insight:** Plate character (oceanic/continental) ≠ crust composition. Need per-pixel crust type.

**Sub-Stages:**
- 1.1: Electrostatic plate generation
- 1.2: Boundary refinement (roughening)
- 1.3: Island removal
- 1.4: Motion assignment & boundary classification
- 1.5: **Crust type assignment** (new!) based on plate type + boundary proximity
- 1.6: **Age gradient from ridges** (new!) oceanic crust ages away from spreading centers

---

## Stage 2: Geologic Provinces

**Input:** TectonicLayers, PlanetaryParams
**Output:** `GeologyLayers`

```rust
pub struct GeologyLayers {
    pub province_id: TerrainMap<u16>,
    pub province_type: TerrainMap<GeologicProvince>,
    pub intensity: TerrainMap<f32>,              // 0.0-1.0 strength
    pub volcanic_activity: TerrainMap<f32>,      // 0.0-1.0 current activity
    pub formation_age_ma: TerrainMap<f32>,       // When feature formed
}

pub enum GeologicProvince {
    // Active Convergent (19 types total now)
    CollisionOrogen,        // Himalayas (active, continent-continent)
    PaleoOrogen,           // Appalachians (ancient, eroded, in plate interior)
    VolcanicArc,           // Andes (oceanic-continental), Aleutians (oceanic-oceanic)
    AccretionaryWedge,
    ForearcBasin,
    BackarcBasin,
    OceanTrench,

    // Stable Continental
    CratonShield,
    CratonPlatform,
    ExtendedCrust,

    // Oceanic
    MidOceanRidge,
    AbyssalPlain,
    FractureZone,
    OceanicPlateau,

    // Large Igneous Provinces
    ContinentalFloodBasalt,
    OceanicPlateau,         // Duplicate? Check
    HotspotTrack,

    // Rifts (future)
    ContinentalRift,
}
```

**PaleoOrogen Implementation:**
```rust
pub struct ProvinceCharacteristics {
    pub province_type: GeologicProvince,
    pub age_ma: f32,
    pub erosion_factor: f32,  // Exponential decay: e^(-age/300)
}

// Generation: Scan plate interiors for old collision zones
// Placement: >500 km from active boundaries
// Age: 200-600 Ma
// Erosion: Reduces elevation by 60-80%
```

**Sub-Stages:**
- 2.1: Collision orogens (active only)
- 2.2: Large igneous provinces
- 2.3: Arc systems (now handles BOTH oceanic-oceanic AND oceanic-continental!)
- 2.4: Stable continental
- 2.5: Continental rifts (deferred)
- 2.6: Oceanic domains
- 2.7: **Paleo-orogens** (new!) ancient mountains in plate interiors

---

## Stage 3: Elevation & Terrain

**Input:** TectonicLayers, GeologyLayers, PlanetaryParams
**Output:** `ElevationLayers`

```rust
pub struct ElevationLayers {
    pub elevation_m: TerrainMap<f32>,        // -11000 to +8848
    pub roughness: TerrainMap<f32>,          // 0.0-1.0 PRIMARY DATA (not computed!)
    pub is_ocean: TerrainMap<bool>,          // elevation < sea_level

    // Cached computations (computed ONCE, stored)
    pub slope_degrees: TerrainMap<f32>,      // From elevation
    pub aspect_degrees: TerrainMap<f32>,     // From elevation
    pub curvature: TerrainMap<f32>,          // From elevation
}
```

**Elevation Formula:**
```rust
elevation = isostatic_base + tectonic_contribution + erosion_adjustment

// Isostatic base from crust type/thickness
isostatic_base = match crust_type {
    Oceanic => -4000.0 + (lithosphere_thickness - 100) * 10,
    Continental => 0.0 + (lithosphere_thickness - 35) * 20,
}

// Tectonic contribution from provinces
tectonic = match province_type {
    CollisionOrogen => 4000 to 8000 (intensity-based),
    PaleoOrogen => (4000 to 8000) * erosion_factor,  // Eroded!
    VolcanicArc => 2000 to 4000,
    OceanTrench => -7000 to -11000,
    MidOceanRidge => -2500 (elevated above abyssal),
    AbyssalPlain => -5000,
    // ...
}
```

**Roughness Purpose:** Multi-scale terrain generation. High roughness (0.9) = jagged peaks, low roughness (0.2) = smooth highlands.

**Why Cache Slope?**
- Computation: 195 seconds for 1800x900 map
- Storage: 6.5 MB
- Used by: Stage 4 (orographic rain), Stage 5 (watersheds), Stage 6 (biomes), Stage 8 (landslides)
- **Decision: CACHE IT**

---

## Stage 4: Climate Modeling

**Input:** ElevationLayers, PlanetaryParams, TectonicLayers
**Output:** `ClimateLayers`

```rust
pub struct ClimateLayers {
    pub temperature_c: TerrainMap<f32>,
    pub precipitation_mm: TerrainMap<f32>,
    pub wind_direction_deg: TerrainMap<f32>,
    pub wind_speed_ms: TerrainMap<f32>,
    pub ocean_current_velocity: TerrainMap<Vec2>,  // For ocean pixels
    pub koppen_zone: TerrainMap<KoppenClimate>,
}
```

**Physics-Based Temperature:**
```rust
// Base from solar insolation
solar_temp = (insolation / STEFAN_BOLTZMANN)^0.25 - 273.15
solar_temp += planetary_params.greenhouse_warming_k

// Latitude adjustment
lat_factor = cos(latitude * π/180)
temp_base = solar_temp * lat_factor

// Elevation lapse rate
temp = temp_base - (elevation_m * 0.0065)  // -6.5°C per km

// Maritime effect (if near ocean)
if distance_to_ocean < 500_km {
    temp += ocean_moderation_factor
}
```

**Ocean Currents (NEW!):**
- Driven by wind patterns + Coriolis effect
- Warm currents: Increase coastal temps by +5-10°C
- Cold currents: Decrease by -5-10°C (Atacama, Namib deserts)

**Sub-Stages:**
- 4.1: Solar insolation patterns
- 4.2: Atmospheric circulation (Hadley cells, westerlies)
- 4.3: **Ocean currents** (Coriolis + wind-driven)
- 4.4: Precipitation (orographic + circulation)
- 4.5: Seasonal variation

---

## Stage 5: Hydrology

**Input:** ElevationLayers, ClimateLayers
**Output:** `HydrologyLayers`

```rust
pub struct HydrologyLayers {
    pub watershed_id: TerrainMap<u32>,
    pub river_flow_m3s: TerrainMap<f32>,      // 0 = no river
    pub is_lake: TerrainMap<bool>,            // SEPARATE from Stage 3!
    pub water_table_depth_m: TerrainMap<f32>,
    pub soil_moisture: TerrainMap<f32>,       // 0.0-1.0
}
```

**Water Classification:**
```rust
// Stage 3 output
is_ocean = elevation < sea_level

// Stage 5 output
is_lake = endorheic_basin && precipitation > evaporation

// Combined query helper
fn is_water(x, y) -> bool {
    elevation_layers.is_ocean[x][y] || hydrology_layers.is_lake[x][y]
}
```

**Key Insight:** Lakes can be above sea level (Lake Titicaca at 3812m!). Need separate flag.

**Sub-Stages:**
- 5.1: Watershed delineation (from elevation)
- 5.2: River network generation
- 5.3: Lake placement (endorheic basins)
- 5.4: Groundwater modeling

---

## Stage 6: Biomes

**Input:** ClimateLayers, ElevationLayers, HydrologyLayers
**Output:** `BiomeLayers`

```rust
pub struct BiomeLayers {
    pub biome_type: TerrainMap<BiomeType>,
    pub vegetation_density: TerrainMap<f32>,
    pub biodiversity_index: TerrainMap<f32>,
}

pub enum BiomeType {
    // Whittaker classification
    TropicalRainforest,
    TropicalSeasonalForest,
    Desert,
    Grassland,
    TemperateForest,
    BorealForest,
    Tundra,
    Alpine,
    // ... 15 total
}
```

**Whittaker Diagram:** Based on temperature + precipitation

---

## Stage 7: Resources & Deposits

**Input:** GeologyLayers, TectonicLayers, HydrologyLayers
**Output:** `ResourceLayers`

```rust
pub struct ResourceLayers {
    pub mineral_flags: TerrainMap<MineralFlags>,  // Quick bitflags
    pub deposits: HashMap<(usize, usize), Vec<DepositInfo>>,  // Detailed data
}

pub struct DepositInfo {
    pub mineral: MineralType,
    pub deposit_type: DepositType,  // GEOLOGICAL CONTEXT!
    pub ore_grade_percent: f64,
    pub tonnage_mt: f64,            // Million tonnes
    pub depth_m: f32,
    pub accessibility: f32,         // 0.0-1.0
}

pub enum DepositType {
    // Igneous-hosted
    Porphyry,              // Cu-Au-Mo at volcanic arcs (Andes)
    Kimberlite,            // Diamonds in cratons (South Africa)
    Layered,               // PGE-Cr in mafic intrusions

    // Hydrothermal
    VMS,                   // Cu-Zn at mid-ocean ridges
    SEDEX,                 // Pb-Zn in sedimentary basins
    OrogenicGold,          // Au in collision zones

    // Sedimentary
    Placer,                // Au-Sn in rivers (from erosion)
    BIF,                   // Fe in ancient ocean floors
    Evaporite,             // Salt, gypsum in closed basins
    Coal,                  // Swamp burial + compression

    // Metamorphic
    Skarn,                 // Cu-Fe-W at intrusion contacts

    // ... 15 total deposit types
}

pub enum MineralType {
    Gold, Silver, Copper, Iron, Coal, Oil, Uranium, Diamonds, Salt, Phosphate, ...
}
```

**Geological Context Examples:**
- **Porphyry Copper**: Only at volcanic arcs where oceanic plate subducts
- **VMS (Volcanogenic Massive Sulfide)**: Only at mid-ocean ridges (black smokers)
- **Placer Gold**: Only in rivers downstream of orogenic gold sources
- **Kimberlite Diamonds**: Only in cratonic shields >2.5 Ga old
- **Coal**: Only in low-lying swamps (elevation < 200m, high precipitation)

**Sub-Stages:**
- 7.1: Igneous-hosted deposits (from geology)
- 7.2: Hydrothermal deposits (from boundaries + heat flow)
- 7.3: Sedimentary deposits (from climate + hydrology)
- 7.4: Placer enrichment (from erosion + rivers)

---

## Stage 8: Natural Hazards

**Input:** TectonicLayers, GeologyLayers, ClimateLayers, HydrologyLayers
**Output:** `HazardLayers`

```rust
pub struct HazardLayers {
    pub earthquake_risk: TerrainMap<f32>,      // 0.0-1.0
    pub volcanic_risk: TerrainMap<f32>,
    pub flood_risk: TerrainMap<f32>,
    pub landslide_risk: TerrainMap<f32>,
    pub tsunami_risk: TerrainMap<f32>,         // Coastal only
}
```

**Risk Calculations:**
- **Earthquake**: High at convergent boundaries, moderate at transform, low at divergent
- **Volcanic**: High at volcanic arcs + hotspots, zero elsewhere
- **Flood**: High in low-lying areas with high precipitation + river flow
- **Landslide**: High slope + high precipitation + weak geology
- **Tsunami**: Coastal pixels within 100 km of ocean trenches

---

## Stage 9: Settlement Suitability

**Input:** ElevationLayers, ClimateLayers, HydrologyLayers, ResourceLayers, HazardLayers
**Output:** `SettlementLayers`

```rust
pub struct SettlementLayers {
    pub habitability: TerrainMap<f32>,         // 0.0-1.0 overall score
    pub agriculture_yield: TerrainMap<f32>,
    pub fresh_water_access: TerrainMap<f32>,
    pub resource_value: TerrainMap<f32>,
    pub defensibility: TerrainMap<f32>,
    pub settlement_clusters: Vec<SettlementCluster>,
}

pub struct SettlementCluster {
    pub center: (usize, usize),
    pub population_capacity: u64,
    pub settlement_type: SettlementType,  // Village, Town, City, Capitol
}
```

**Habitability Formula:**
```rust
habitability =
    0.3 * climate_score +        // Temp 10-30°C, precip 500-2000mm
    0.2 * water_access +          // Rivers, lakes, groundwater
    0.2 * agriculture_potential + // Flat land, good soil, water
    0.15 * resource_proximity +   // Minerals, lumber, game
    0.1 * defensibility -         // Hills, rivers as barriers
    0.05 * hazard_risk            // Subtract earthquake/flood risk
```

---

## Stage 10: Kingdoms & Cultures

**Input:** SettlementLayers, BiomeLayers, GeologyLayers
**Output:** `CulturalLayers`

```rust
pub struct CulturalLayers {
    pub kingdom_id: TerrainMap<u16>,
    pub culture_id: TerrainMap<u16>,
    pub language_family: TerrainMap<u16>,
    pub trade_routes: Vec<TradeRoute>,
}

pub struct Kingdom {
    pub id: u16,
    pub capital: (usize, usize),
    pub territory: Vec<(usize, usize)>,
    pub population: u64,
    pub government_type: GovernmentType,
    pub resources: Vec<MineralType>,
}
```

**Generation:**
- Start with settlement clusters from Stage 9
- Grow territories using Voronoi-like expansion weighted by:
  - Habitability (prefer good land)
  - Natural barriers (mountains, rivers block expansion)
  - Resource value (fight over rich areas)
- Trade routes connect kingdoms via navigable rivers + passes

---

## Pipeline Linearity Verification

**Strict Rule:** No stage modifies previous stages' output

| Stage | Reads From | Writes To | Violations? |
|-------|------------|-----------|-------------|
| 0 | User input | PlanetaryParams | ✅ None |
| 1 | Stage 0 | TectonicLayers | ✅ None |
| 2 | Stages 0, 1 | GeologyLayers | ✅ None |
| 3 | Stages 0, 1, 2 | ElevationLayers | ✅ None |
| 4 | Stages 0, 1, 3 | ClimateLayers | ✅ None |
| 5 | Stages 3, 4 | HydrologyLayers | ✅ None (is_lake separate!) |
| 6 | Stages 3, 4, 5 | BiomeLayers | ✅ None |
| 7 | Stages 1, 2, 5 | ResourceLayers | ✅ None |
| 8 | Stages 1, 2, 4, 5 | HazardLayers | ✅ None |
| 9 | Stages 3, 4, 5, 7, 8 | SettlementLayers | ✅ None |
| 10 | Stages 2, 6, 9 | CulturalLayers | ✅ None |

**Result:** Fully linear! User can inject custom data at any stage boundary.

---

## Memory Estimates (1800x900 map)

| Stage | Primary Data | Size | Notes |
|-------|--------------|------|-------|
| 1 | 6 TerrainMaps | ~40 MB | u16 + 5× f32 |
| 2 | 5 TerrainMaps | ~33 MB | u16 + 4× f32 |
| 3 | 5 TerrainMaps | ~33 MB | Includes cached slope |
| 4 | 6 TerrainMaps | ~40 MB | 5× f32 + enum |
| 5 | 5 TerrainMaps | ~33 MB | u32 + 3× f32 + bool |
| 6 | 3 TerrainMaps | ~20 MB | enum + 2× f32 |
| 7 | 1 TerrainMap + HashMap | ~10 MB | Sparse deposits |
| 8 | 5 TerrainMaps | ~33 MB | 5× f32 |
| 9 | 5 TerrainMaps + clusters | ~33 MB | 5× f32 + Vec |
| 10 | 3 TerrainMaps + kingdoms | ~20 MB | 3× u16 + Vec |

**Total:** ~295 MB for full 10-stage pipeline (1800×900 map)

**1 TerrainMap<f32>:** 1800 × 900 × 4 bytes = 6.48 MB

---

## Implementation Roadmap

### Phase 1: Refactor Stage 1 (Tectonics)
- Add `crust_type: TerrainMap<CrustType>` per-pixel assignment
- Add `crust_age_ma: TerrainMap<f32>` gradient from ridges
- Add `lithosphere_thickness_km: TerrainMap<f32>`

### Phase 2: Refactor Stage 2 (Geology)
- Fix oceanic-continental arc systems (extend line 220)
- Add `PaleoOrogen` province type
- Implement paleo-orogen generation (plate interior scan)

### Phase 3: Refactor Stage 3 (Elevation)
- Add `roughness: TerrainMap<f32>` as PRIMARY input to terrain generation
- Compute and cache `slope_degrees` from elevation
- Rename `is_water` → `is_ocean` (from `elevation < 0`)

### Phase 4: New Stage 4 (Climate)
- Implement physics-based temperature from `PlanetaryParams.insolation`
- Add ocean current simulation
- Implement orographic precipitation using cached slope

### Phase 5: New Stage 5 (Hydrology)
- Watershed delineation from elevation
- River network generation
- **Add `is_lake` flag** (separate from Stage 3's `is_ocean`)

### Phase 6-10: Remaining Stages
- Implement biomes (Stage 6)
- Implement deposit types with geological context (Stage 7)
- Implement hazards (Stage 8)
- Implement settlements (Stage 9)
- Implement kingdoms (Stage 10)

---

## Key Takeaways

1. **Per-Pixel Crust Type**: Foundation for correct geology (oceanic-continental convergence)
2. **Roughness is Primary**: Not computed from elevation, used to GENERATE elevation detail
3. **Slope is Cached**: Computed once, stored (saves 195 seconds across stages)
4. **Split Water Flags**: `is_ocean` (Stage 3) + `is_lake` (Stage 5) = strictly linear
5. **Physics-Based Climate**: Driven by planetary params (insolation, rotation, tilt)
6. **Deposit Types**: Geological realism (Porphyry at arcs, VMS at ridges, placer in rivers)
7. **Paleo-Orogens**: Ancient eroded mountains in plate interiors (Appalachians)
8. **Strictly Linear**: No backwards modification, enables manual injection at any stage

---

## Outstanding Questions

1. **Stage 2 Province Count**: Do we have duplicate `OceanicPlateau` in LIP and Oceanic categories?
2. **Roughness Generation**: How is roughness initially assigned in Stage 3? From province type?
3. **Ocean Current Model**: Simple wind-driven or full thermohaline circulation?
4. **Settlement Clustering**: Voronoi or hierarchical (village → town → city)?

---

## Next Steps

1. **Review this architecture** for completeness
2. **Decide on outstanding questions**
3. **Begin Phase 1 implementation** (refactor Stage 1 for per-pixel crust)
