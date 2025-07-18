# Geoforge Project

**Realistic geological and climate modeling for procedural world generation**

## Project Overview
Geoforge is a Rust library for generating scientifically-inspired geological features, climate patterns, and biomes for procedural world generation. It starts with tectonic plate simulation and builds up through geological domains, elevation, climate, and biomes.

## Key Features
- Tectonic plate generation using Voronoi or region growing algorithms
- Global projection support with longitude wraparound and polar regions
- Deterministic generation with seed-based random generation
- Multiple export formats: Binary, PNG visualization, and GeoTIFF
- Performance optimized with fast distance calculations

## Current Status
- ✅ Tectonic plates implementation (foundation layer)
- 🚧 Geological domains, elevation, climate, biomes (planned)
- Working on spherical coordinate system improvements

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