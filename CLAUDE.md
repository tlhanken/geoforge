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
- ✅ Electrostatic physics simulation (production ready)
- ✅ Earth-like size variety and natural boundaries
- ✅ Clean, optimized codebase with electrostatic-only approach
- ✅ Nested output directory structure (`outputs/examples/`)
- ✅ Comprehensive planetary parameters system
- ✅ Stellar luminosity and insolation calculations
- ✅ Physics-based orbital mechanics (inverse square law)

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

**1.4 Plate Motion & Boundary Classification** 🔄 NEXT
**Foundation:** Assign motion vectors and classify boundary interactions for geology
- **Motion vector assignment** - Each plate gets direction (azimuth) and velocity (cm/year)
  - Use physics-based angular velocity on sphere or simplified random assignment
  - Earth-realistic velocities: ~1-10 cm/year
- **Relative motion calculation** - Compute motion at each boundary pixel pair
- **Boundary type classification** - Convergent, divergent, or transform based on relative motion
  - Convergent: Plates colliding (angle < 45° from head-on)
  - Divergent: Plates spreading apart (angle < 45° from opposite)
  - Transform: Plates sliding past (angle near 90°)
- **Plate character assignment** - Oceanic vs continental designation
  - Affects subduction zones (oceanic subducts under continental)
  - Could use size heuristic (larger plates = more continental) or random
- **Boundary segment storage** - Store type and characteristics per boundary segment
- **Visualization export** - Color-coded boundary types for validation
- **Note:** Critical foundation for Stage 2 (geologic provinces) and Stage 3 (elevation)

### **Stage 2: Geologic Provinces** ⏳ PLANNED
**Foundation:** Use tectonic plates to determine geological characteristics and provinces

**2.1 Orogenic Belts (Mountain-building zones)**
- **Collision orogeny** - Continental plate convergence (Himalayas-style)
- **Subduction orogeny** - Oceanic-continental convergence (Andes-style)  
- **Accretionary orogeny** - Terrane accretion and exotic block collision
- **Extensional orogeny** - Core complex formation in rifting zones

**2.2 Large Igneous Provinces (LIPs)**
- **Continental flood basalts** - Massive volcanic provinces (Deccan Traps-style)
- **Oceanic plateaus** - Underwater volcanic provinces
- **Hotspot tracks** - Volcanic island chains and seamount trails
- **Dyke swarms** - Radiating intrusive networks

**2.3 Arc and Basin Systems**
- **Volcanic arcs** - Active subduction zone volcanism
- **Forearc basins** - Sedimentary basins between trench and arc
- **Backarc basins** - Extensional basins behind volcanic arcs
- **Backarc ridges** - Spreading centers in backarc regions

**2.4 Stable Continental Regions**
- **Cratons/Shields** - Ancient, stable continental cores (>1.5 Ga)
- **Platforms** - Stable cratonic areas with thin sedimentary cover
- **Intracratonic basins** - Subsided regions within stable cratons

**2.5 Extensional Zones**
- **Continental rifts** - Active extension and normal faulting
- **Extended crust** - Thinned continental crust from extension
- **Transitional crust** - Continent-ocean boundary zones

**2.6 Oceanic Domains**
- **Mid-ocean ridges** - Active seafloor spreading centers
- **Abyssal plains** - Deep oceanic basins with sediment cover
- **Oceanic fracture zones** - Transform fault systems
- **Deep ocean trenches** - Subduction zone depocenters

### **Stage 3: Elevation Generation** ⏳ PLANNED
**Foundation:** Use geologic provinces to generate realistic elevation
- **3.1** Mountain range generation based on orogenic belts
- **3.2** Ocean floor depth modeling (ridges, trenches, abyssal plains)
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
- Working on: `spherical_coords_tectonics`
- Main branch: `main`
- Recent focus: Spherical coordinate system and tectonics improvements

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