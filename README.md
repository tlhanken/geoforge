# Geoforge 🌍

**Realistic tectonic plate generation using electrostatic physics simulation**

> **Geoforge v2 (branch `v2-implementation`):** a solar system and its world, generated layer by
> layer (tectonics → geology → elevation → climate → hydrology → biomes), with every layer
> exportable to PNG/GIS files, editable or hand-drawn, and re-importable; downstream layers
> regenerate from your edits.
>
> ```bash
> cargo build --release -p geoforge-cli
> ./target/release/geoforge-v2 new worlds/demo --seed 42 --habitable
> ./target/release/geoforge-v2 generate worlds/demo --render
> ```
>
> See [docs/V2_PIPELINE.md](docs/V2_PIPELINE.md), [docs/V2_REALISM.md](docs/V2_REALISM.md) and the roadmap, [docs/V2_ROADMAP.md](docs/V2_ROADMAP.md).
> The v1 crate below (`crates/bins/geoforge`) is kept as the reference implementation.

Geoforge is a Rust library for generating scientifically-inspired tectonic plates for procedural world generation. It uses electrostatic physics simulation to create natural plate boundaries with Earth-like size variety.

## Features

- ⚡ **Electrostatic Physics Simulation** - Point charges reach equilibrium for natural plate spacing
- 🌍 **Earth-like Size Variety** - 6000x+ size ratios from superplates to micro-plates
- 🧲 **Natural Boundaries** - Clean, curved boundaries from physics-based positioning
- 🌐 **Global Projection Support** - Proper handling of longitude wraparound and polar regions  
- 🎲 **Deterministic Generation** - Reproducible results with seed-based random generation
- ⚡ **Performance Optimized** - Fast physics simulation and memory-efficient algorithms
- 📁 **Multiple Export Formats** - Binary, PNG visualization, and GeoTIFF for GIS applications
- 🎯 **Organized Output** - Clean file organization in dedicated output directories
- 🪐 **Planetary Parameters** - Configurable planetary size, gravity, orbital mechanics, and stellar luminosity
- ☀️ **Accurate Insolation** - Physics-based solar radiation calculations for climate modeling

## Quick Start

Add this to your `Cargo.toml`:

```toml
[dependencies]
geoforge = { version = "0.1", features = ["export-full"] }
```

Or with specific export features:

```toml
[dependencies]
geoforge = { version = "0.1", features = ["export-png", "export-tiff"] }
```

Generate your first world with the complete Stage 1 pipeline:

```rust
use geoforge::WorldMap;

// Create a new world (1800x900 at 0.2° resolution)
let mut world = WorldMap::new(1800, 900, 42)?;

// Stage 1.1: Generate tectonic plates using electrostatic physics
world.generate_tectonics(20)?;

// Stage 1.2: Refine boundaries for realistic irregular edges
world.refine_boundaries(None)?;

// Stage 1.3: Remove islands to ensure plate contiguity
let island_stats = world.remove_islands(None)?;

// Export the result
#[cfg(feature = "export-png")]
world.export_tectonics_png("outputs", "tectonics.png")?;

// Get statistics about the generated plates
if let Some(stats) = world.get_tectonic_stats() {
    for (plate_id, stat) in stats.iter().take(5) {
        println!("Plate {}: {:.1}% of surface ({} km²)",
                 plate_id, stat.percentage, stat.area_km2);
    }
}
```

## Examples

Run the example application:

```bash
cargo run --features export-png
```

This will execute the complete Stage 1 pipeline:
1. **Generate** 20 tectonic plates using electrostatic physics simulation
2. **Refine** boundaries to add realistic irregularity
3. **Remove islands** to ensure all plates are contiguous
4. **Export** the result as PNG visualization

## Export Formats

Geoforge supports multiple export formats for different use cases:

### 📊 **Binary Format** (.bin)
- Raw u16 data (little-endian) for high-performance applications
- Always available, no feature flags required
- Ideal for further processing or embedding in applications

### 🖼️ **PNG Visualization** (.png) 
- Color-coded bitmap showing plate boundaries
- Each plate gets a distinct color for easy identification
- Perfect for validation and presentation
- Requires `export-png` feature

### 🗺️ **GeoTIFF** (.tiff)
- Industry-standard format with proper georeferencing
- Compatible with GIS software (QGIS, ArcGIS, etc.)
- Includes coordinate system and projection information
- Requires `export-tiff` feature

### Export Examples

```rust
// Export all formats
generator.export_all("outputs", "world_name")?;

// Export specific formats
generator.export_binary("outputs", "plates.bin")?;
generator.export_png("outputs", "visualization.png")?;     // if export-png enabled
generator.export_geotiff("outputs", "georef.tiff")?;       // if export-tiff enabled
```

## How It Works

### Stage 1: Tectonic Foundation (Complete)

The complete tectonic generation pipeline consists of three stages:

**Stage 1.1: Core Plate Generation**
1. **⚡ Charge Placement** - Random point charges distributed on the sphere
2. **🎯 Size Variety** - Power-law charge distribution creates Earth-like size hierarchy
3. **🧲 Physics Simulation** - Charges repel until reaching equilibrium
4. **🗺️ Boundary Generation** - Voronoi diagram from equilibrium positions
5. **✨ Smoothing** - Optional geodesic-aware boundary smoothing

**Stage 1.2: Boundary Refinement**
- Adds realistic irregularity to plate edges using 3D noise-based domain warping
- Creates natural jagged boundaries instead of perfect Voronoi edges
- Configurable noise parameters for different levels of roughness

**Stage 1.3: Island Removal**
- Ensures all plates are contiguous using connected component analysis
- Removes isolated plate fragments created during refinement
- Reassigns island pixels to surrounding plates intelligently

## Pipeline Overview

Geoforge implements a scientifically-grounded generation pipeline:

- **⭐ Stage 0: Stellar Systems** ⏳ - Star generation, luminosity, habitable zones (defaults to Sun-like)
- **🔥 Stage 1: Tectonic Foundation** ✅ - Complete pipeline with plate generation, boundary refinement, and island removal
- **🏔️ Stage 2: Geologic Provinces** ⏳ - Orogenic belts, LIPs, cratons, oceanic domains
- **⛰️ Stage 3: Elevation Generation** ⏳ - Mountains, ocean floors, continental margins  
- **🌤️ Stage 4: Climate Modeling** ⏳ - Trade winds, temperature, precipitation patterns
- **🌿 Stage 5: Biome Generation** ⏳ - Realistic biome distribution from climate data
- **🏞️ Stage 6: Hydrological Systems** ⏳ - Rivers, lakes, watersheds from elevation
- **🔮 Stage 7: Advanced Features** ⏳ - Resources, hazards, settlement patterns

*See CLAUDE.md for detailed roadmap and geological specifications*


## Development

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

## License

Licensed under MIT OR Apache-2.0
