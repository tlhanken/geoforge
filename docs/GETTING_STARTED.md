# Getting started

geoforge makes a solar system and generates its world layer by layer: tectonics, geology,
elevation, winds and climate, erosion, rivers, biomes, resources. Every layer can be
exported as images, edited or redrawn, and imported back; the layers below follow. Then you
can ask about any place and any day: what is there, the weather, the sky.

## Build

You need Rust (stable, 1.85 or newer; `rustup` installs it). With Nix, `nix develop` gives a
shell with everything. [`just`](https://github.com/casey/just) is optional but handy.

```bash
cargo build --release -p geoforge-cli        # or: just build
./target/release/geoforge --help
```

The command is `target/release/geoforge`; put it on your `PATH` or use `just`, which builds
it first.

## Your first world

```bash
geoforge new worlds/demo --seed 42           # a system drawn from the Milky Way and its most Earth-like body
geoforge generate worlds/demo --render       # every layer, with previews in worlds/demo/renders
geoforge show worlds/demo                    # the star, planets, and the world: air, water, field, tides, calendar
```

Or in one go: `just world 42` (makes `worlds/world-42`).

- `--habitable` forces a Sun-like star with an Earth-like world; `--galaxy elliptical`
  (or `dwarf`, `starburst`) draws from another kind of galaxy; `geoforge census` lists
  seeds by what they give.
- `--archetype earth` (or `mars`, `titan`, …) starts from a real body, in the real Solar
  System at its last northern spring equinox (Earth's: 2026-03-20 14:46 UTC), so `sky`
  shows the real planets and eclipses with dates; `--body my.toml` from your own
  description of one.
- `--width 1024` makes a smaller, faster grid (the default is 2048 × 1024 cells).

The previews in `renders/` show each layer (elevation, biomes, climate, rivers, plates …)
and the derived views (Köppen climates, soils, deposits, hazards, sea ice, pressure).
`geoforge render worlds/demo koppen hazards` redraws chosen ones.

## Ask about a place

```bash
geoforge place worlds/demo --lat 45 --lon 10                 # climate, soil, deposits, hazards within 100 km
geoforge weather worlds/demo --lat 45 --lon 10 --day 100     # two weeks of daily weather from day 100
geoforge almanac worlds/demo --lat 45 --lon 10 --year 1      # a year day by day: sun, moons, weather, soil, droughts
geoforge almanac worlds/demo --lat 45 --lon 10 --calendar lunar   # the same days in its biggest moon's months
geoforge weather-map worlds/demo --day 100                   # isobars, fronts, rain and storms that day
geoforge sky worlds/demo --day 100 --to worlds/demo/renders  # the night sky, eclipses, a star chart
geoforge positions worlds/demo --day 100                     # where every body is
```

Days are Earth days from the present, which is the world's northern spring equinox; any day
works, far in the future or the past, and gives the same answer every time. `--json` prints
machine-readable output.

The almanac dates its days in a calendar: `--calendar twelfths` (twelve equal months),
`lunar` or `lunar:<moon>` (new moon to new moon), or `world` (the default): the world's own,
if `bodies/<name>/calendar.json` gives one, for example

```json
{"name": "Telari reckoning", "present_year": 1042, "era": "AR",
 "week": ["Moonday", "Tideday", "Ashday", "Starday", "Restday"],
 "months": [{"name": "Thaw", "days": 30}, {"name": "Bloom", "days": 31}, {"name": "Frost"}]}
```

Months are whole local days; the last takes what the year leaves, so it gains a day in leap
years. The day, the year and the second never change with the calendar.

## Take its pulse

```bash
geoforge report worlds/demo                                  # worlds/demo/report/index.html
geoforge report worlds/earth --place "London:51.5,-0.1" --place "Cairo:30.0,31.2"
```

One page to read at a glance: every map in themed groups with its key (classes, ramps,
chances), the world's numbers against Earth's, three places (picked by climate unless you
name them) with their hazards and two weeks of weather, the weather map, and the night sky
with a key to star sizes, colours and planets. The folder opens in a browser or publishes
as is.

To set worlds against an observed one (authored products holding real data, like the
real-Earth worlds of `scripts/real-earth/ab.sh`):

```bash
geoforge compare worlds/real worlds/model --probe "Lima:-12.05,-77.05" --probe "Alps:44..48,6..15" \
  --source "temperature=ERA-Interim 2 m, 2005"                # worlds/real/compare/index.html
```

Each observed product gets the observed map beside each model's, with one key, and the
numbers between them (bias, RMSE, pattern correlation; wind direction agreement; class
agreement), plus the values at the probes.

## Edit a layer

```bash
geoforge export worlds/demo rough_elevation --caps   # worlds/demo/edit/rough_elevation.png (+ polar caps)
#   … paint in GIMP, Krita, QGIS or Blender …
geoforge import worlds/demo rough_elevation          # only what you changed, as a patch
geoforge generate worlds/demo --render               # climate, rivers and biomes follow
```

Heights are 16-bit greyscale; biomes, rocks and plates are exact colours listed in the
sidecar files. `--scale 2` exports at twice the grid's size for detail; any 2:1 image comes
back. `--caps` adds north and south polar views without the map's stretching. `geoforge
status` shows what is edited; `geoforge reset worlds/demo rough_elevation` drops the edit.

To move mountains or coasts edit `rough_elevation` (the climate reads it); edit `biome`,
`lithology` or `plates` to repaint those. `--authored` replaces a layer with your own map
outright.

## Where next

- [V2_PIPELINE.md](V2_PIPELINE.md): every command, the layers and products, editing in
  depth (Blender round trip), and the models in short.
- [V2_REALISM.md](V2_REALISM.md): how each layer compares with Earth and other worlds.
- [V2_ROADMAP.md](V2_ROADMAP.md): what comes next.
- `geoforge help <command>`: every option.
