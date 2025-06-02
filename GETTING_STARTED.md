# Getting Started with Geoforge

Welcome to your new Geoforge repository! Here's how to get started:

## 🚀 Quick Start

1. **Navigate to the project directory:**
   ```bash
   cd C:\Users\tlhan\Downloads\geoforge
   ```

2. **Run with full export features:**
   ```bash
   cargo run --features export-full
   ```

3. **Run the simple example:**
   ```bash
   cargo run --example simple
   ```

4. **Run tests:**
   ```bash
   cargo test
   ```

## 📁 Output Organization

All generated files are organized in the `outputs/` directory:

```
outputs/
├── random_world.bin           # Different every run
├── random_world.png
├── random_world.tiff
├── reproducible_world.bin      # Always same (seed 42)
├── reproducible_world.png
├── reproducible_world.tiff
├── voronoi_method.bin          # Algorithm comparison
├── region_growing_method.bin
└── detailed_world.*            # High resolution
```

## 📁 Project Structure

```
geoforge/
├── Cargo.toml          # Package configuration
├── README.md           # Project documentation  
├── src/
│   ├── lib.rs          # Main library code
│   └── main.rs         # Example application
├── examples/
│   └── simple.rs       # Simple usage example
└── .gitignore          # Git ignore file
```

## 🔧 Development Commands

```bash
# Run with optimizations (faster)
cargo run --release

# Run specific example
cargo run --example simple

# Check code without running
cargo check

# Format code
cargo fmt

# Run clippy linter
cargo clippy

# Generate documentation
cargo doc --open
```

## 📊 Expected Output

When you run `cargo run --features export-full`, you should see:
- 4 different examples of tectonic plate generation
- Performance comparisons between algorithms  
- Detailed world statistics
- Export progress for each format
- File organization summary

## 🖼️ Viewing Results

**Visual Validation:**
- **🖼️ PNG files** - Open with any image viewer to see colored plate boundaries
- **🗺️ TIFF files** - Import into QGIS, ArcGIS, or other GIS software
- **📊 Binary files** - For programmatic processing (raw u16 data)

**Recommended viewers:**
- **Windows:** Photos app, Paint, or any image viewer for .png
- **QGIS:** Free GIS software - perfect for .tiff files with coordinates
- **Web browsers:** Can display .png files by dragging and dropping

## 🛠️ Next Steps

1. **Test the basic functionality** - Run the examples to make sure everything works
2. **Explore the API** - Look at `src/lib.rs` to understand the available functions
3. **Add geological domains** - This will be the next layer in the worldbuilding pipeline
4. **Customize parameters** - Try different grid sizes, plate counts, and seeds

## 📖 Usage Examples

```rust
// Basic generation
let mut gen = TectonicPlateGenerator::new(1800, 900, 15)?;
let plates = gen.generate("region_growing", true)?;

// Reproducible generation  
let mut gen = TectonicPlateGenerator::with_seed(1800, 900, 15, 42)?;
let plates = gen.generate("voronoi", false)?;

// Get statistics
let stats = gen.get_plate_stats();
for (id, stat) in stats {
    println!("Plate {}: {:.1}%", id, stat.percentage);
}
```

## 🎯 Performance Notes

- Default resolution: 0.2° per pixel (1800×900 grid)
- Generation time: ~3 seconds for full resolution
- Memory usage: ~13MB for 3600×1800 grid
- Algorithm choice: "region_growing" is more realistic, "voronoi" is faster

Enjoy building worlds! 🌍
