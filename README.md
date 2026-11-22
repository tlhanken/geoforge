# Geoforge

**Procedural worlds you can edit at any layer.**

Geoforge draws a solar system from a galaxy and generates its world layer by layer —
plate tectonics, geology, elevation, winds and climate, erosion, rivers and lakes, biomes,
ore and soils — each checked against real Earth. Any layer can be exported to images and
GIS files, edited or redrawn, and imported back; the layers below follow. Then you can ask
about any place and day: its climate, deposits and hazards, the weather, the storms
crossing the map, the night sky and its eclipses.

```bash
cargo build --release -p geoforge-cli
./target/release/geoforge new worlds/demo --seed 42
./target/release/geoforge generate worlds/demo --render
./target/release/geoforge weather worlds/demo --lat 45 --lon 10
```

- [Getting started](docs/GETTING_STARTED.md)
- [The pipeline: commands, layers, editing, models](docs/V2_PIPELINE.md)
- [Realism: how each layer compares with Earth](docs/V2_REALISM.md)
- [Roadmap](docs/V2_ROADMAP.md)
- [The crates](crates/libs/README.md)

Licensed under MIT or Apache-2.0.
