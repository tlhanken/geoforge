# Geoforge Pipeline Architecture - Complete Redesign

## Executive Summary

**Problem Identified:** The current Stage 2 (geology) architecture has fundamental flaws:
1. **Missing per-pixel crust type** - Critical for elevation, climate, and biomes
2. **Province regions stored as overlays** - Not a proper per-pixel layer
3. **Plate type ≠ Crust type** - Confusion between plate character and crust composition
4. **Context-dependent provinces** - Same province type means different things (island arc vs continental arc)

**Proposed Solution:** Redesign the entire pipeline with proper data flow architecture from Stage 1 → Stage 6.

---

## Current Architecture (Flawed)

### Data Structures

```rust
// Per-pixel layers (TerrainMap<T>)
tectonics: TerrainMap<u16>        // Plate IDs only
elevation: TerrainMap<f32>        // Future
temperature: TerrainMap<f32>      // Future
precipitation: TerrainMap<f32>    // Future
biomes: TerrainMap<u8>            // Future

// NOT per-pixel! (Vector of regions)
geology: Vec<ProvinceRegion>      // PROBLEM: Overlays, not fundamental data

// Metadata (per-plate)
TectonicMetadata {
    plate_seeds: Vec<PlateSeed>,
    plate_stats: HashMap<u16, PlateStats>,  // Contains PlateType (Oceanic/Continental)
    plate_boundaries: Vec<BoundarySegment>,
}
```

### Critical Gaps

| What's Missing | Why It Matters | Blocks Which Stages |
|----------------|----------------|---------------------|
| **Per-pixel crust type** | Continental crust is buoyant (high elevation), oceanic crust is dense (low elevation) | Stage 3 (Elevation) |
| **Per-pixel lithosphere age** | Older oceanic crust is colder, denser, deeper | Stage 3 (Elevation) |
| **Per-pixel geology ID** | Need to query "what province is at pixel (x,y)?" for elevation/climate | Stage 3, 4, 5 |
| **Marine vs terrestrial flag** | Climate over ocean vs land is completely different | Stage 4 (Climate) |
| **Sediment thickness** | Affects elevation, resource distribution | Stage 3, 7 |

---

## Proposed Architecture: Layered Per-Pixel Design

### Philosophy: **Everything is per-pixel, with metadata for groupings**

```
┌─────────────────────────────────────────────────────────┐
│  Fundamental Layers (Per-Pixel TerrainMaps)            │
├─────────────────────────────────────────────────────────┤
│  Stage 1: Tectonic Foundation                           │
│    - plate_id: TerrainMap<u16>         [EXISTING]      │
│    - crust_type: TerrainMap<CrustType> [NEW]           │
│    - crust_age: TerrainMap<f32>        [NEW]           │
│    - lithosphere_thickness: TerrainMap<f32> [FUTURE]   │
├─────────────────────────────────────────────────────────┤
│  Stage 2: Geologic Provinces                            │
│    - province_id: TerrainMap<u16>      [NEW]           │
│    - province_intensity: TerrainMap<f32> [NEW]         │
├─────────────────────────────────────────────────────────┤
│  Stage 3: Elevation                                     │
│    - elevation: TerrainMap<f32>        [EXISTING]      │
│    - bathymetry: TerrainMap<f32>       [FUTURE]        │
├─────────────────────────────────────────────────────────┤
│  Stage 4: Climate                                       │
│    - temperature: TerrainMap<f32>      [EXISTING]      │
│    - precipitation: TerrainMap<f32>    [EXISTING]      │
│    - wind_u/v: TerrainMap<f32>         [FUTURE]        │
├─────────────────────────────────────────────────────────┤
│  Stage 5: Biomes                                        │
│    - biome_id: TerrainMap<u8>          [EXISTING]      │
├─────────────────────────────────────────────────────────┤
│  Stage 6: Hydrology                                     │
│    - watershed_id: TerrainMap<u32>     [FUTURE]        │
│    - river_flow: TerrainMap<f32>       [FUTURE]        │
└─────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────┐
│  Metadata Layers (Lookups and Groupings)               │
├─────────────────────────────────────────────────────────┤
│  TectonicMetadata                                       │
│    - plate_stats: HashMap<u16, PlateStats>             │
│    - plate_boundaries: Vec<BoundarySegment>            │
├─────────────────────────────────────────────────────────┤
│  GeologyMetadata [NEW]                                  │
│    - provinces: HashMap<u16, ProvinceInfo>             │
│    - province_stats: HashMap<u16, ProvinceStats>       │
├─────────────────────────────────────────────────────────┤
│  ElevationMetadata [FUTURE]                             │
│    - sea_level: f32                                     │
│    - elevation_stats: ElevationStatistics              │
└─────────────────────────────────────────────────────────┘
```

---

## Stage-by-Stage Design

### Stage 1: Tectonic Foundation

#### Outputs (Per-Pixel)

```rust
/// Stage 1 per-pixel outputs
pub struct TectonicLayers {
    /// Which tectonic plate (1-255, 0 = unassigned)
    pub plate_id: TerrainMap<u16>,

    /// Type of crustal material at this pixel
    pub crust_type: TerrainMap<CrustType>,

    /// Age of crust in millions of years (Ma)
    /// - Oceanic: 0-200 Ma (gets subducted)
    /// - Continental: 500-4000 Ma (ancient, preserved)
    pub crust_age: TerrainMap<f32>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CrustType {
    /// Oceanic crust - Basaltic, dense (~7 km thick), young (<200 Ma)
    Oceanic,

    /// Continental crust - Granitic, buoyant (~35 km thick), ancient
    Continental,

    /// Transitional crust - Island arcs, volcanic buildups, rifted margins
    /// Intermediate density and composition
    Transitional,
}
```

#### Initial Assignment Rules

**At plate generation time:**
1. All pixels start with crust type matching plate character:
   - "Oceanic" plates → `CrustType::Oceanic`
   - "Continental" plates → `CrustType::Continental`

2. Crust age initialization:
   - Oceanic: Random 0-180 Ma (simulates spreading history)
   - Continental: Random 1000-3000 Ma (ancient cratons)

**During Stage 2 (geology):**
3. Volcanic arcs CREATE new crust:
   - Oceanic-oceanic arcs → `CrustType::Transitional` (building toward continental)
   - Oceanic-continental arcs → Stays `CrustType::Continental` (adds to existing)

4. Mid-ocean ridges set crust age = 0 Ma (brand new crust)

5. Continental rifts → Transition from Continental → Transitional → Oceanic (over time)

---

### Stage 2: Geologic Provinces

#### Purpose
Identify **what tectonic process is active** at each pixel, which determines how elevation forms.

#### Outputs (Per-Pixel)

```rust
/// Stage 2 per-pixel outputs
pub struct GeologyLayers {
    /// Province ID (0-65535, 0 = unassigned/default)
    /// Lookup in GeologyMetadata.provinces for details
    pub province_id: TerrainMap<u16>,

    /// Intensity/activity level at this pixel (0.0-1.0)
    /// - Distance from plate boundary
    /// - Convergence rate
    /// - Used to modulate elevation in Stage 3
    pub intensity: TerrainMap<f32>,
}

/// Metadata for each province
pub struct GeologyMetadata {
    /// Province information by ID
    pub provinces: HashMap<u16, ProvinceInfo>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProvinceInfo {
    /// Unique ID for this province
    pub id: u16,

    /// Type of province
    pub province_type: GeologicProvince,

    /// Which plate(s) this province belongs to
    pub plate_ids: Vec<u16>,

    /// Crust type this province sits on (can be mixed)
    pub primary_crust_type: CrustType,

    /// Tectonic context (if subduction-related)
    pub tectonic_context: Option<TectonicContext>,

    /// Typical characteristics
    pub characteristics: ProvinceCharacteristics,

    /// Number of pixels
    pub pixel_count: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum TectonicContext {
    /// Oceanic plate subducting under oceanic plate (island arcs)
    OceanicOceanic { subducting: u16, overriding: u16 },

    /// Oceanic plate subducting under continental plate (Andes-style)
    OceanicContinental { subducting: u16, overriding: u16 },

    /// Continental plates colliding (Himalayan-style)
    ContinentalContinental { plate_a: u16, plate_b: u16 },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProvinceCharacteristics {
    /// Base elevation intensity (0.0-1.5+)
    pub elevation_intensity: f64,

    /// Terrain roughness (0.0-1.0)
    pub roughness: f64,

    /// Width in km
    pub width_km: f64,

    /// Convergence/spreading rate (cm/year)
    pub tectonic_rate: f64,

    /// Expected elevation class
    pub elevation_class: ElevationClass,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ElevationClass {
    /// Deep ocean (< -3000m) - Trenches, abyssal plains
    DeepMarine,

    /// Shallow ocean (0 to -3000m) - Continental shelf, shallow seas
    ShallowMarine,

    /// Lowland (0-500m) - Coastal plains, sedimentary basins
    Lowland,

    /// Moderate elevation (500-2000m) - Hills, plateaus
    Moderate,

    /// High elevation (2000-5000m) - Mountains
    High,

    /// Extreme elevation (>5000m) - Collision orogens (Himalayas, Andes)
    Extreme,
}
```

#### Province Assignment Algorithm

**Input:** Plate map, crust type map, boundaries, plate stats
**Output:** Province ID map, intensity map, metadata

**Algorithm:**
1. **Initialize base layers** (foundation):
   ```rust
   for each pixel (x, y):
       crust = crust_type_map[x, y]
       if crust == Continental:
           province_id[x, y] = assign_continental_base(plate_id, distance_to_boundary)
           // Craton (far from boundary) or Platform (medium distance)
       else if crust == Oceanic:
           province_id[x, y] = assign_oceanic_base(plate_id, distance_to_boundary)
           // AbyssalPlain (far from boundary)
   ```

2. **Overlay boundary features** (active tectonics):
   ```rust
   for each boundary:
       if boundary.type == Convergent:
           if oceanic_subducts_under_oceanic:
               create_island_arc_system()
               // Trench, AccretionaryWedge, ForearcBasin, VolcanicArc, BackarcBasin
               // All on Transitional crust
           else if oceanic_subducts_under_continental:
               create_continental_arc_system()
               // Trench (oceanic side), Wedge/Forearc/Arc/Backarc (continental side)
               // Arc on Continental crust
           else if continental_collides_continental:
               create_collision_orogen()
               // CollisionOrogen (both sides, Continental crust)

       else if boundary.type == Divergent:
           if oceanic_oceanic:
               create_mid_ocean_ridge()
               // Set crust_age = 0 Ma (new crust forming)
           else if continental_continental:
               create_continental_rift()
               // May transition to Transitional crust (proto-oceanic)
   ```

3. **Overlay volcanic features** (hotspots, flood basalts):
   ```rust
   // These sit ON TOP of base layers
   generate_hotspot_tracks()
   generate_flood_basalts()
   ```

4. **Calculate intensity field**:
   ```rust
   for each pixel with province_id:
       intensity[x, y] = calculate_intensity(
           distance_to_boundary,
           convergence_rate,
           province_type
       )
   ```

---

### Stage 3: Elevation Generation

#### Inputs
- `plate_id: TerrainMap<u16>`
- `crust_type: TerrainMap<CrustType>` ⭐ **KEY INPUT**
- `crust_age: TerrainMap<f32>` ⭐ **KEY INPUT**
- `province_id: TerrainMap<u16>`
- `intensity: TerrainMap<f32>`
- `GeologyMetadata` (province info)
- `TectonicMetadata` (plate info)

#### Outputs
- `elevation: TerrainMap<f32>` (meters relative to sea level)

#### Algorithm

```rust
for each pixel (x, y):
    // Step 1: Base elevation from crust type
    base_elevation = match crust_type[x, y] {
        Oceanic => {
            age = crust_age[x, y]
            // Older oceanic crust is colder, denser, deeper
            // New crust at ridges: -2500m
            // Old crust at trenches: -5500m
            -2500.0 - (age * 15.0)  // Subsides ~15m per million years
        },
        Continental => {
            // Continental crust is buoyant
            // Default elevation ~200m (low plains)
            200.0
        },
        Transitional => {
            // Island arcs, volcanic buildups
            // Above sea level but not high mountains
            50.0
        },
    }

    // Step 2: Modify by province type
    province = get_province_info(province_id[x, y])
    if province exists:
        elevation_modifier = province.characteristics.elevation_intensity
        elevation_adjustment = match province.province_type {
            CollisionOrogen => {
                // Extreme uplift
                intensity[x, y] * 6000.0 * elevation_modifier
            },
            VolcanicArc => {
                context = province.tectonic_context
                match context {
                    OceanicOceanic => {
                        // Island arc (low elevation)
                        intensity[x, y] * 2000.0 * elevation_modifier
                    },
                    OceanicContinental => {
                        // Continental arc (high elevation - Andes)
                        intensity[x, y] * 5000.0 * elevation_modifier
                    },
                }
            },
            MidOceanRidge => {
                // Elevated seafloor
                2000.0 * elevation_modifier  // Rises to ~-2500m
            },
            OceanTrench => {
                // Deepest ocean
                -6000.0 * intensity[x, y]  // Down to -11000m
            },
            ForearcBasin | BackarcBasin => {
                // Context-dependent!
                if province.primary_crust_type == Continental:
                    // Terrestrial basin (Great Valley, Altiplano)
                    -500.0 * intensity[x, y]
                else:
                    // Marine basin
                    -2000.0 * intensity[x, y]
            },
            AbyssalPlain => 0.0,  // Base already set
            Craton | Platform => 100.0,  // Slight uplift
            // ... other types
        }

        elevation[x, y] = base_elevation + elevation_adjustment
    else:
        elevation[x, y] = base_elevation

    // Step 3: Add noise for terrain variation
    elevation[x, y] += generate_terrain_noise(
        x, y,
        province.characteristics.roughness,
        seed
    )
```

**Key Insight:** By having per-pixel `crust_type` and `crust_age`, Stage 3 can properly calculate:
- Oceanic depth (older = deeper)
- Continental buoyancy (always elevated)
- Volcanic arc height (depends on tectonic context)

---

### Stage 4: Climate

#### Inputs
- `elevation: TerrainMap<f32>` ⭐ **Determines land vs ocean**
- `crust_type: TerrainMap<CrustType>` (for land/ocean discrimination if elevation ambiguous)
- `latitude` (from projection)
- `planetary_params` (insolation, axial tilt)

#### Outputs
- `temperature: TerrainMap<f32>`
- `precipitation: TerrainMap<f32>`
- `wind_u, wind_v: TerrainMap<f32>` (future)

#### Key Decision
```rust
for each pixel (x, y):
    is_ocean = (elevation[x, y] < 0.0) || (crust_type[x, y] == Oceanic && elevation[x, y] < 100.0)

    if is_ocean:
        temperature[x, y] = calculate_ocean_temperature(lat, elevation, ocean_currents)
    else:
        temperature[x, y] = calculate_land_temperature(lat, elevation, distance_to_coast)
```

---

### Stage 5: Biomes

#### Inputs
- `elevation: TerrainMap<f32>`
- `temperature: TerrainMap<f32>`
- `precipitation: TerrainMap<f32>`
- `crust_type: TerrainMap<CrustType>` (for marine vs terrestrial discrimination)

#### Outputs
- `biome_id: TerrainMap<u8>`

#### Decision Tree
```rust
if elevation[x, y] < 0.0:
    // Marine biomes
    biome = classify_marine_biome(elevation, temperature)
else:
    // Terrestrial biomes (Whittaker diagram)
    biome = classify_terrestrial_biome(temperature, precipitation)
```

---

## Implementation Roadmap

### Phase 1: Add Fundamental Layers (Stage 1 Enhancement)
**Status:** FOUNDATION - Do this first

**Files to modify:**
- `src/map/world.rs` - Add new TerrainMap fields
- `src/tectonics/plates.rs` - Add CrustType enum
- `src/tectonics/generator.rs` - Initialize crust_type and crust_age maps

**Steps:**
1. Define `CrustType` enum
2. Add `crust_type: TerrainMap<CrustType>` to WorldMap
3. Add `crust_age: TerrainMap<f32>` to WorldMap
4. Initialize during plate generation:
   - Match plate character → initial crust type
   - Random age assignment
5. Export functions (PNG, binary) to visualize crust

**Tests:**
- Verify all pixels have crust type assigned
- Verify oceanic plates have oceanic crust (initially)
- Verify continental plates have continental crust
- Verify crust ages are reasonable (0-200 Ma oceanic, 1000-3000 Ma continental)

---

### Phase 2: Refactor Stage 2 to Per-Pixel
**Status:** DEPENDS ON PHASE 1

**Files to modify:**
- `src/map/world.rs` - Replace `geology: Vec<ProvinceRegion>` with TerrainMaps
- `src/geology/provinces.rs` - Add GeologyMetadata structure
- `src/geology/generator.rs` - Output per-pixel maps instead of region vectors

**Changes:**
```rust
// OLD
pub geology: Option<Vec<ProvinceRegion>>,

// NEW
pub geology_province_id: Option<TerrainMap<u16>>,
pub geology_intensity: Option<TerrainMap<f32>>,
pub geology_metadata: Option<GeologyMetadata>,
```

**Steps:**
1. Implement province ID assignment algorithm
2. Calculate intensity field
3. Build GeologyMetadata lookup table
4. Update all export functions
5. Update all tests

**Tests:**
- Verify every pixel has province_id
- Verify province metadata matches pixel counts
- Verify intensity field is smooth (no discontinuities)
- Test province lookup by pixel coordinate

---

### Phase 3: Fix Oceanic-Continental Convergence
**Status:** DEPENDS ON PHASE 2

Now that we have:
- Per-pixel crust type
- Tectonic context in metadata
- Proper province assignment

The oceanic-continental fix becomes simple:

```rust
// Arc systems generator
match (subducting_crust, overriding_crust) {
    (Oceanic, Oceanic) => {
        // Island arc system
        create_volcanic_arc(context: TectonicContext::OceanicOceanic)
    },
    (Oceanic, Continental) => {
        // Continental arc system (ANDES!)
        create_volcanic_arc(context: TectonicContext::OceanicContinental)
    },
    _ => {}
}
```

**Tests:**
- Test oceanic-continental boundaries generate arcs
- Test arc sits on continental crust
- Test trench sits on oceanic crust
- Test tectonic context is recorded correctly

---

### Phase 4: Implement Stage 3 (Elevation)
**Status:** DEPENDS ON PHASE 2 & 3

With proper crust type and province data, elevation generation is straightforward:
- Read crust_type → base elevation
- Read crust_age → oceanic depth adjustment
- Read province_id → lookup characteristics
- Apply elevation formula
- Add terrain noise

**Tests:**
- Verify oceanic crust is below sea level
- Verify continental crust is above sea level
- Verify collision orogens have extreme elevation
- Verify volcanic arcs have context-appropriate elevation

---

## Migration Strategy

### Option A: Clean Break (Recommended)
1. Implement Phase 1 & 2 in a separate branch
2. Update ALL tests to new API
3. Merge when complete
4. Old API deprecated

**Pros:** Clean architecture, no technical debt
**Cons:** Large PR, significant testing effort

### Option B: Gradual Migration
1. Add new fields alongside old ones
2. Implement both old and new code paths
3. Migrate one stage at a time
4. Deprecate old API when all stages migrated

**Pros:** Incremental, less risky
**Cons:** Maintains technical debt temporarily, confusing dual APIs

---

## Open Questions

1. **Should crust type be mutable during Stage 2?**
   - YES: Volcanic arcs create new crust (Transitional)
   - NO: Keep Stage 1 immutable, handle in Stage 3

2. **Should we track lithosphere thickness?**
   - Affects elevation (thicker = more buoyant)
   - Affects subduction angle
   - May be overkill for v1.0

3. **How to handle crust type at boundaries?**
   - Sharp transition or gradual blend?
   - Matters for elevation continuity

4. **Should province_id be u8 (255 max) or u16 (65535 max)?**
   - Large worlds with many boundaries could exceed 255
   - u16 safer but uses 2x memory

5. **Export format implications?**
   - More layers = more export files
   - Need combined visualizations?
   - GeoTIFF multi-band support?

---

## Recommendation

**Proceed with Phase 1 first** - Add crust_type and crust_age to Stage 1. This is:
- Non-breaking (adds new fields)
- Validates the concept
- Provides immediate value (can visualize crust distribution)
- Required foundation for everything else

Once Phase 1 is solid, we can tackle Phase 2 (refactor Stage 2 to per-pixel).

---

## Summary

**Current Problem:** Missing fundamental per-pixel data (crust type) causes downstream issues.

**Solution:** Redesign as layered per-pixel architecture where each stage produces TerrainMap outputs.

**Key Insight:** Crust type is fundamental - it determines elevation, climate, and biomes. It must be per-pixel and set in Stage 1.

**Next Step:** Implement Phase 1 (add crust_type and crust_age to Stage 1) as foundation for all other improvements.
