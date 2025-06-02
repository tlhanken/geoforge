# Getting Started with Geoforge

Welcome to your new Geoforge repository! Here's how to get started:

## 🚀 Quick Start

1. **Navigate to the project directory:**
   ```bash
   cd C:\Users\tlhan\Downloads\geoforge
   ```

2. **Run the main example:**
   ```bash
   cargo run
   ```

3. **Run the simple example:**
   ```bash
   cargo run --example simple
   ```

4. **Run tests:**
   ```bash
   cargo test
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

When you run `cargo run`, you should see:
- 4 different examples of tectonic plate generation
- Performance comparisons between algorithms
- Detailed world statistics
- Export of binary data file

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
