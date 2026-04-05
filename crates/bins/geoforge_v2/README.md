# GeoForge CLI Tool

GeoForge is a procedural world generation CLI tool designed for generating world procedurally and reproducibly.

## Project Structure

- **`crates/bins/geoforge_v2`**: This CLI tool.
- **`crates/libs/types`**: Core domain models and coordinate systems.

## Getting Started

To run the GeoForge V2 generator from the project root:

```bash
cargo run -p geoforge_v2 -- generate [OPTIONS]
```

### Examples

- **Basic generation:**
  ```bash
  cargo run -p geoforge_v2 -- generate
  ```
- **Custom seed and output:**
  ```bash
  cargo run -p geoforge_v2 -- generate --seed 42 --output-directory ./my_results
  ```
- **Specify output format:**
  ```bash
  cargo run -p geoforge_v2 -- generate -f sqlite
  ```

## CLI Configuration

The `generate` command supports the following options:

- `-s, --seed <SEED>`: Seed for random number generation.
- `--from-layer <LAYER>`: Starting layer (default: `planetary-body`).
- `--to-layer <LAYER>`: Ending layer (default: `political-entities`).
- `-d, --output-directory <DIR>`: Directory for generated files (default: `outputs`).
- `-f, --output-file-format <FORMAT>`: Output format (`bin`, `sqlite`, `png`, `bin-and-png`, `sqlite-and-png`, `geotiff`).
- `-v`: Increase CLI verbosity.

For more details, run:
```bash
cargo run -p geoforge_v2 -- --help
```
