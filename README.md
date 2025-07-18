# Geoforge 🌍

**Realistic geological and climate modeling for procedural world generation**

Geoforge is a Rust library for generating scientifically-inspired geological features, climate patterns, and biomes for procedural world generation. It starts with tectonic plate simulation and builds up through geological domains, elevation, climate, and biomes.

## Features

- 🗺️ **Tectonic Plate Generation** - Realistic plate boundaries using Voronoi or region growing algorithms
- 🌐 **Global Projection Support** - Proper handling of longitude wraparound and polar regions  
- 🎲 **Deterministic Generation** - Reproducible results with seed-based random generation
- ⚡ **Performance Optimized** - Fast distance calculations and memory-efficient algorithms
- 🧪 **Scientifically Inspired** - Based on real geological and climate processes
- 📁 **Multiple Export Formats** - Binary, PNG visualization, and GeoTIFF for GIS applications
- 🎯 **Organized Output** - Clean file organization in dedicated output directories

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

Generate your first world:

```rust
use geoforge::TectonicPlateGenerator;

// Generate a world with 15 tectonic plates at 0.2° resolution
let mut generator = TectonicPlateGenerator::with_seed(1800, 900, 15, 42)?;
let plate_map = generator.generate("region_growing", true)?;

// Export in all available formats to organized directory
generator.export_all("outputs", "my_world")?;

// Get statistics about the generated plates
let stats = generator.get_plate_stats();
for (plate_id, stat) in stats {
    println!("Plate {}: {:.1}% of surface ({} km²)", 
             plate_id, stat.percentage, stat.area_km2);
}
```

## Examples

Run the example application:

```bash
cargo run
```

This will generate a tectonic plate map with 15 plates and export it in multiple formats.

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

## World Building Pipeline

The planned full pipeline:

1. **✅ Tectonic Plates** - Foundation layer defining major crustal boundaries
2. **🚧 Geological Domains** - Rock types, orogenies, basins based on plate interactions  
3. **🚧 Elevation** - Height maps from geological processes
4. **🚧 Climate** - Temperature, precipitation, wind patterns
5. **🚧 Biomes** - Ecosystems based on climate and geography


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
