# Geoforge 🌍

**Realistic geological and climate modeling for procedural world generation**

Geoforge is a Rust library for generating scientifically-inspired geological features, climate patterns, and biomes for procedural world generation. It starts with tectonic plate simulation and builds up through geological domains, elevation, climate, and biomes.

## Features

- 🗺️ **Tectonic Plate Generation** - Realistic plate boundaries using Voronoi or region growing algorithms
- 🌐 **Global Projection Support** - Proper handling of longitude wraparound and polar regions  
- 🎲 **Deterministic Generation** - Reproducible results with seed-based random generation
- ⚡ **Performance Optimized** - Fast distance calculations and memory-efficient algorithms
- 🧪 **Scientifically Inspired** - Based on real geological and climate processes

## Quick Start

Add this to your `Cargo.toml`:

```toml
[dependencies]
geoforge = "0.1"
```

Generate your first world:

```rust
use geoforge::TectonicPlateGenerator;

// Generate a world with 15 tectonic plates at 0.2° resolution
let mut generator = TectonicPlateGenerator::with_seed(1800, 900, 15, 42)?;
let plate_map = generator.generate("region_growing", true)?;

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

This will generate several example worlds and show performance comparisons between different algorithms.

## Algorithms

### Tectonic Plate Generation

- **Voronoi Method**: Fast, clean boundaries based on distance to nearest seed point
- **Region Growing**: More realistic, irregular boundaries with natural variation

Both methods support:
- Proper longitude wraparound at ±180°
- Configurable resolution (default 0.2° = ~22km at equator)
- Boundary smoothing for natural-looking edges
- Motion vectors for each plate (direction and speed)

## World Building Pipeline

The planned full pipeline:

1. **✅ Tectonic Plates** - Foundation layer defining major crustal boundaries
2. **🚧 Geological Domains** - Rock types, orogenies, basins based on plate interactions  
3. **🚧 Elevation** - Height maps from geological processes
4. **🚧 Climate** - Temperature, precipitation, wind patterns
5. **🚧 Biomes** - Ecosystems based on climate and geography

## Technical Details

- **Coordinate System**: WGS84 Geographic (EPSG:4326)
- **Data Format**: 16-bit unsigned integers for plate IDs
- **Performance**: ~3 seconds for 3600×1800 grid (0.1° resolution)
- **Memory Usage**: ~13MB for 3600×1800 grid

## Export Formats

- Raw binary data (little-endian u16)
- GeoTIFF metadata for GIS applications
- Statistics and analysis data

## Development

```bash
# Run tests
cargo test

# Run with optimizations
cargo run --release

# Run specific examples
cargo run --example detailed_analysis
```

## License

Licensed under either of Apache License, Version 2.0 or MIT license at your option.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
