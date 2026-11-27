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

With [`just`](https://github.com/casey/just) (`just` lists every task):

```bash
just quick 7                      # a 1024×512 world from seed 7, every layer, with previews
just report worlds/quick-7        # one page: every map with its key, places, weather, the sky
just view worlds/quick-7          # look around its universe in a browser (127.0.0.1:8137)
just place worlds/quick-7 45 10   # what is at a place
just weather worlds/quick-7 45 10 100      # two weeks of weather from day 100
just sky worlds/quick-7 30        # the night sky on day 30
just body mars                    # a real body instead (earth, mars, venus, titan, …)
just edit worlds/quick-7 rough_elevation   # export a layer; edit the PNG; then
just apply worlds/quick-7 rough_elevation  # import it and regenerate what follows
just ab                           # the Earth A/B page (needs the observed data, see below)
```

Worlds go in `worlds/` (ignored by git); each command prints where it wrote its renders,
report or edit files. `just world 42` makes the full 2048×1024 version (a few minutes).

- [Getting started](docs/GETTING_STARTED.md)
- [What is next](docs/TODO.md)
- [The Earth A/B check and its data](scripts/real-earth/README.md)
- [The pipeline: commands, layers, editing, models](docs/V2_PIPELINE.md)
- [Realism: how each layer compares with Earth](docs/V2_REALISM.md)
- [Roadmap](docs/V2_ROADMAP.md)
- [The crates](crates/libs/README.md)

Licensed under MIT or Apache-2.0.
