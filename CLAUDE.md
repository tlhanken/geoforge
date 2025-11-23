# Geoforge Project

**Realistic geological and climate modeling for procedural world generation**

## Project Overview
Geoforge is a Rust library for generating scientifically-inspired geological features, climate patterns, and biomes for procedural world generation. It starts with tectonic plate simulation and builds up through geological domains, elevation, climate, and biomes.

## Key Features
- Electrostatic physics simulation for natural plate boundaries
- Earth-like size variety with 6000x+ ratios from superplates to micro-plates
- Global projection support with longitude wraparound and polar regions
- Deterministic generation with seed-based random generation
- Multiple export formats: Binary, PNG visualization, and GeoTIFF
- Performance optimized with fast physics simulation

## Current Status
- ✅ **Stage 1: Tectonic Foundation** - COMPLETE (production ready)
  - Electrostatic physics simulation for natural plate boundaries
  - Earth-like size variety with 6000x+ ratios
  - Boundary refinement with domain warping
  - Island removal for plate contiguity
  - Plate motion assignment and boundary classification
  - Boundary visualization export
- 🔧 **Stage 2: Geologic Provinces** - IN PROGRESS (polishing phase)
  - 18 province types implemented (2.1-2.6 complete)
  - ~1,950 lines of code across 3 modules
  - Tests passing
  - Ready for polishing and validation
- ✅ Comprehensive planetary parameters system
- ✅ Stellar luminosity and insolation calculations
- ✅ Physics-based orbital mechanics (inverse square law)
- ✅ 82 tests passing (57 unit, 20 integration, 5 doc)

## Development Workflow

**IMPORTANT:** When implementing new features, always follow this workflow:

### 1. API Design First
- **Start with interface design** - Define public APIs and function signatures first
- **Review and proposal** - Present the API design for review before implementation
- **Consider ergonomics** - Ensure APIs are intuitive and follow Rust best practices
- **Document early** - Write doc comments during API design phase

### 2. Test-Driven Development
- **Write tests second** - After API approval, write comprehensive test cases
- **Cover edge cases** - Include boundary conditions, error cases, and integration scenarios
- **Test organization:**
  - Unit tests: `#[cfg(test)] mod tests` in source files
  - Integration tests: `/tests/stageN_*_integration.rs` for pipeline stages
  - Examples: `/examples/` for user-facing demonstrations
- **Test before implementation** - Tests should fail initially (Red-Green-Refactor)

### 3. Implementation
- **Implement to tests** - Write code to satisfy the test requirements
- **Maintain test coverage** - Ensure all tests pass before marking complete
- **Refactor iteratively** - Improve code quality while keeping tests green

### 4. Code Review
- **Ensure all code is sensible** - Make sure code is sensible and not bloated.  Implement what is needed efficiently, nothing more.
- **Meets style guidelines** - Ensure all code uses best practices for rust
- **Remove unnessary code** - Simplify down to just needed code, and clean it up

### 5. Workflow Summary
```
1. Design API → Review → Approve
2. Write Tests → Review
3. Implement → Tests Pass → Done
4. Code Review
```

### 6. Test Coverage Standards
- **Unit tests** - All public functions and methods
- **Integration tests** - Full pipeline workflows for each major stage
- **Edge cases** - Boundary conditions, empty inputs, extreme values
- **Error handling** - All error paths and failure modes
- **Determinism** - Verify the same seed reproduces the same output every time.

## Pipeline Roadmap

### **Stage 0: Stellar System Generation** ⏳ FUTURE
**Foundation:** Generate the stellar environment that planets orbit
- **0.1** Multiple star systems (binary, trinary stars)
- **0.2** Stellar evolution and lifecycle modeling
- **0.3** Stellar classification and spectral types (O, B, A, F, G, K, M)
- **0.4** Variable star luminosity and solar cycles
- **0.5** Stellar habitable zone calculations
- **0.6** Stellar metallicity and planet formation effects
- **Note:** Currently defaults to Sun-like star (G-type, 1.0 solar luminosity)

### **Stage 1: Tectonic Foundation** ✅ COMPLETE
**1.1 Core Plate Generation** ✅ COMPLETE
- Electrostatic physics simulation for plate generation
- Natural plate boundaries with Earth-like size distribution
- Export system (PNG, binary, GeoTIFF)
- Deterministic generation with seeds

**1.2 Boundary Refinement** ✅ COMPLETE
**Foundation:** Post-process Voronoi boundaries to add realistic irregularity
- **Boundary roughening** - Add natural jaggedness to plate edges
- **Noise-based perturbation** - Perlin/simplex noise for organic shapes
- **Smoothing options** - Configurable smoothing to avoid over-jaggedness
- **Boundary-type variation** - Different roughness for convergent vs divergent vs transform boundaries
- **Seed-based consistency** - Maintain determinism in roughening process
- **Note:** Keeps Stage 1.1 generation isolated and unaffected

**1.3 Island Removal / Plate Contiguity** ✅ COMPLETE
**Foundation:** Ensure all plates are contiguous and remove isolated fragments
- **Connected component analysis** - Identify non-contiguous plate regions using flood fill
- **Fragment detection** - Find isolated plate "islands" separated from main body
- **Size comparison** - Determine smaller vs larger portions of fragmented plates
- **Smart reassignment** - Assign smaller fragments to most common neighboring plate
- **Deterministic processing** - Consistent behavior for reproducible results
- **Note:** Applied after boundary refinement to clean up artifacts

**1.4 Plate Motion & Boundary Classification** ✅ COMPLETE
**Foundation:** Assign motion vectors and classify boundary interactions for geology
- **Motion vector assignment** - Each plate gets direction (azimuth) and velocity (cm/year)
  - Random assignment with deterministic seeding (1-10 cm/year Earth-realistic)
  - Physics-based angular velocity option available
  - Automatically assigned during tectonic generation
- **Relative motion calculation** - Compute motion at each boundary pixel pair
- **Boundary type classification** - Convergent, divergent, or transform based on relative motion
  - Convergent: Plates colliding (moving toward each other)
  - Divergent: Plates spreading apart (moving away from each other)
  - Transform: Plates sliding past (parallel motion)
- **Plate character assignment** - Oceanic vs continental designation
  - Size-based heuristic: larger plates = continental, smaller = oceanic
  - Automatically assigned and stored in metadata
- **Boundary segment storage** - Store type and characteristics per boundary segment
- **Visualization export** - Color-coded boundary types (red=convergent, blue=divergent, green=transform)
- **Note:** Foundation for Stage 2 (geologic provinces) and Stage 3 (elevation)

### **Stage 2: Geologic Provinces** 🔧 IN PROGRESS (Polishing)
**Foundation:** Use tectonic plates to determine geological characteristics and provinces

**Implementation Status:**
- ✅ **2.1 Collision Orogens** - Continent-continent mountain building (1 type: CollisionOrogen)
- ✅ **2.2 Large Igneous Provinces** - Continental flood basalts, oceanic plateaus, hotspot tracks (3 types)
- ✅ **2.3 Subduction Zone Systems** - Complete oceanic→continental transect (5 types: OceanTrench, AccretionaryWedge, ForearcBasin, VolcanicArc, BackarcBasin)
- ✅ **2.4 Stable Continental Regions** - Cratons/shields, platforms, extended crust (3 types)
- ⏳ **2.5 Continental Rifts** - DEFERRED (will implement when needed)
- ✅ **2.6 Oceanic Domains** - Mid-ocean ridges, abyssal plains, fracture zones, hotspot tracks (4 types)
- ✅ **Hotspot Tracks Refinement** - Geologically accurate visible portions only (2,400 km oceanic, 400 km continental)
- ✅ **Module Refactoring** - Clean naming: `geology/{provinces, orogenic, generator}`
- ✅ **Determinism Bug Fix** - Fixed HashMap iteration order causing non-deterministic generation

**Code Quality:**
- ~1,950 lines of code across 3 modules
- **18 province types implemented** (SubductionOrogen removed - now modeled as VolcanicArc + AccretionaryWedge; IntracratonicBasin too small-scale)
- Tests passing
- Clean module structure with clear naming
- Comprehensive documentation
- Deterministic generation verified

**Polishing Phase: Remaining Actions**

**Phase 1: Critical Fixes & Simplification** ⭐ HIGH PRIORITY
- [x] **Simplify naming conventions** - Remove vague names like "comprehensive"
  - Renamed `comprehensive.rs` → `generator.rs` (matches `tectonics/generator.rs` pattern)
  - Renamed `ComprehensiveGeologyGenerator` → `GeologyGenerator`
  - Renamed `ComprehensiveGeologyConfig` → `GeologyConfig`
  - Module structure now clean: `geology/{provinces, orogenic, generator}`
- [x] **Extract constants** - Replace magic numbers with named constants
  - Replaced all `71.0` km/pixel hardcoded values → dynamic `km_per_pixel()` calculation
  - Uses `MapProjection::km_per_pixel()` for accurate scale at any map size
  - Added defensive min value check (0.01 km/pixel minimum)
- [x] **Fix potential bugs**
  - Added defensive checks in `km_per_pixel()` for zero/negative values
  - All edge cases handled with early returns and empty checks

**Phase 2: Code Quality & Refactoring** ⭐ HIGH PRIORITY
- [ ] **Eliminate code duplication**
  - Refactor `expand_boundary()` and `expand_boundary_toward_plate()` - extract common flood-fill pattern
  - Consolidate repeated plate pixel collection loops (performance issue)
- [ ] **Optimize performance**
  - Build `PlatePixelIndex` once during generation instead of scanning per-plate
  - Reduces O(plates × width × height) to O(width × height + plates)
  - Profile and optimize HashSet→Vec conversions in expansion
- [ ] **Simplify complex functions**
  - Break down large functions (e.g., `generate_arc_systems` at 106 lines)
  - Extract helper methods for clarity
  - Remove unnecessary complexity

**Phase 3: Testing & Validation** ⭐ HIGH PRIORITY
- [ ] **Add full geology integration test**
  - Test `GeologyGenerator` full pipeline (currently only `OrogenicBeltGenerator` tested)
  - Verify all province categories generate correctly
  - Test province layering (foundation → active features)
- [ ] **Add edge case tests**
  - Empty boundary lists
  - Single-pixel plates
  - Polar region provinces
  - Longitude wraparound handling
- [ ] **Visual validation of all 18 province types** ⭐ CRITICAL
  - Generate test worlds with known configurations
  - Export PNG visualizations for each province type
  - Validate color assignments match geological meaning:
    - Collision orogens: High elevation (browns/whites)
    - Subduction orogens: Moderate-high elevation (browns)
    - Oceanic features: Blues (depth-based gradients)
    - Continental stable: Greens/tans (low elevation)
    - Volcanic: Reds/oranges (active features)
  - Verify behavior: width, roughness, intensity all sensible
  - Check for overlaps, gaps, or rendering issues

**Phase 4: Documentation** ⭐ MEDIUM PRIORITY
- [ ] **Add module-level usage examples**
  - Complete example in `provinces.rs`, `orogenic.rs`, `geology/mod.rs`
  - Show common workflows and configurations
- [ ] **Enhance province documentation**
  - Add elevation ranges to each province type (+4000m to +8000m for collision orogens)
  - Add typical width ranges (500-2000 km for collision zones)
  - Add more real-world examples per province type
- [ ] **Document constants and formulas**
  - Explain km/pixel calculations
  - Document width scaling formulas
  - Add references to geological literature where appropriate

**Phase 5: API Polish** ⭐ LOW PRIORITY
- [ ] **Consider API improvements** (optional, evaluate need)
  - `GeologyConfig` boolean flags → enum-based `enabled_stages`?
  - Builder pattern for `ProvinceRegion`?
  - Ergonomics review with fresh eyes

### **Stage 3: Elevation Generation** ⏳ PLANNED
**Foundation:** Use geologic provinces to generate realistic elevation
- **3.1** Mountain range generation based on orogenic belts
- **3.2** Ocean floor depth modeling (ridges, trenches, abyssal plains)
  - **Mid-ocean ridges**: Morphology depends on spreading rate (stored in `convergence_rate` field, negative value)
    - Fast-spreading (>10 cm/yr, width ~60 km): Smooth gentle rise, NO deep rift valley (East Pacific Rise)
    - Slow-spreading (2-5 cm/yr, width ~120 km): Deep central rift valley 1-2 km deep (Mid-Atlantic Ridge)
    - Ultra-slow (<1 cm/yr, width ~150 km): Highly irregular, deepest rift valleys (Gakkel Ridge)
    - All ridges: Elevated 2-3 km above abyssal plain (~2500m depth vs ~5000m)
  - **Ocean trenches**: Narrow (75 km), deepest ocean features (-7000m to -11000m)
- **3.3** Continental shelf and slope definition
- **3.4** Volcanic elevation features from LIPs and arcs
- **3.5** Erosion and sedimentation effects over geological time

### **Stage 4: Climate Modeling** ⏳ PLANNED
**Foundation:** Use elevation + global patterns for climate
- **4.1** Trade winds and prevailing wind patterns (Hadley cells, westerlies)
- **4.2** Global temperature patterns (latitude, elevation, maritime effects)
- **4.3** Precipitation modeling (rain shadows, orographic precipitation)
- **4.4** Ocean current simulation (driven by wind patterns)
- **4.5** Seasonal variation patterns
- **4.6** Climate zone classification (Köppen system)

### **Stage 5: Biome Generation** ⏳ PLANNED
**Foundation:** Use climate + elevation for biome placement
- **5.1** Whittaker biome classification system
- **5.2** Transition zones between biomes
- **5.3** Special biomes (coastal, alpine, desert oases)
- **5.4** Biodiversity hotspot identification

### **Stage 6: Hydrological Systems** ⏳ PLANNED
**Foundation:** Use elevation + precipitation for water flow
- **6.1** River network generation
- **6.2** Lake and wetland placement
- **6.3** Watershed analysis
- **6.4** Groundwater simulation

### **Stage 7: Advanced Features** ⏳ FUTURE
- **7.1** Resource distribution (minerals, oil, etc.)
- **7.2** Natural disaster modeling
- **7.3** Human settlement suitability
- **7.4** Temporal evolution simulation

## Development Commands
```bash
# Run tests
cargo test

# Run with all export features
cargo run --features export-full

# Run with specific features
cargo run --features export-png
cargo run --features export-tiff

# Run with optimizations
cargo run --release --features export-full

# Check without building
cargo check --features export-full
```

## Current Branch
- Working on: `stage2_geologic_provinces`
- Main branch: `main`
- Recent focus: Polishing Stage 2 geologic provinces implementation

## Dependencies
- `rand = "0.8"` - Random number generation
- `image = "0.24"` - PNG export (optional)
- `gdal = "0.15"` - GeoTIFF export (optional)

## Architecture
- Main library: `src/lib.rs`
- Spherical geometry: `src/spherical/geometry.rs`
- Tectonics module: `src/tectonics/`
- Binary example: `src/main.rs`

## Export Features
- `export-png`: Enables PNG visualization export
- `export-tiff`: Enables GeoTIFF export for GIS applications
- `export-full`: Enables all export formats

## Notes
- Working on polar region handling and longitude wraparound
- Recent commits focused on spherical coordinate improvements
- Manual re-write of tectonics module in progress