# Getting started

geoforge makes a solar system and generates its world layer by layer: tectonics, geology,
elevation, winds and climate, erosion, rivers, biomes, resources. Every layer can be
exported as images, edited or redrawn, and imported back; the layers below follow. Then you
can ask about any place and any day: what is there, the weather, the sky.

## Build

You need Rust (stable, 1.85 or newer; `rustup` installs it). With Nix, `nix develop` gives a
shell with everything (the flake is a [Bloomery](https://github.com/actionable-work/bloomery)
workspace: `nix build` makes `geoforge`, `nix build .#geoforge-viewer` the viewer, and
`bloomery check` runs the checks). [`just`](https://github.com/casey/just) is optional but handy.

```bash
cargo build --release -p geoforge-cli        # or: just build
./target/release/geoforge --help
```

The command is `target/release/geoforge`; put it on your `PATH` or use `just`, which builds
it first.

## Your first world

```bash
geoforge new worlds/demo --seed 42           # universe 42 and its default body: a Sun-like star's wet world
geoforge generate worlds/demo --render       # every layer, with previews in worlds/demo/renders
geoforge show worlds/demo                    # the star, planets, and the world: air, water, field, tides, calendar
```

Or in one go: `just world 42` (makes `worlds/world-42`).

- `worlds/demo` is a universe directory: universe 42, and each body made in it at its
  address (`7/h6:0,0,-1:17.0/f12:3,-4,0:5/4` reads universe / galaxy / star / body). The
  default body is in the universe's galaxy most like the Milky Way: the nearest Sun-like
  star with a wet world to a point on that galaxy's Sun-like circle.
- `--at <address>` makes any other galaxy's, system's or body's (a galaxy is taken to its
  default star, a system to its most Earth-like body) in the same directory, and makes it
  the focus every command works on; `--focus <name>` picks another body of the system.
  `geoforge focus worlds/demo` lists the bodies made and switches between them;
  `geoforge visit worlds/demo <star>` goes to a star of the sky. `--at real` is the real
  universe: the Earth, the Solar System, the Milky Way and the Local Group.
- `--archetype earth` (or `mars`, `titan`, …) makes a world on its own from a real body,
  in the real Solar System at its last northern spring equinox (Earth's: 2026-03-20 14:46 UTC), so `sky`
  shows the real planets and eclipses with dates; `--body my.toml` from your own
  description of one.
- `--width 1024` makes a smaller, faster grid (the default is 2048 × 1024 cells).

The previews in `renders/` show each layer (elevation, biomes, climate, rivers, plates …)
and the derived views (Köppen climates, soils, deposits, hazards, sea ice, pressure).
`geoforge render worlds/demo koppen hazards` redraws chosen ones.

## Look around

```bash
geoforge-viewer worlds/demo                  # then open http://127.0.0.1:8137/
geoforge-viewer worlds/demo --host 0.0.0.0   # reachable from a phone or tablet on your network
```

The viewer zooms from the universe's galaxies into a galaxy (its shape, the systems made in
it), the stars around a system, the system's orbits and a body: its true colour on a globe
or a flat map, and every layer and view once it is generated. Search finds names among what
is made or nearby; type an address to go anywhere. A body you find can be made, and its
layers generated, from its panel. Round the home star the galaxy view glows with nebulae
(emission, reflection, planetary, supernova remnants, dark clouds) and shows the galaxy's
globular clusters and centre; each tells about itself and leads to its star.

## Your setting's own materials

A universe directory holds `exotics.toml`, a commented list of fictional materials. Each
names where it is found — in asteroids by type, on worlds by kind, in nebulae by kind,
optionally only round certain stars — with a chance per candidate:

```toml
[[material]]
name = "Xylarium"
note = "a metastable superheavy alloy"
colour = "#7fffd4"

[[material.found]]
in = "asteroid"     # asteroid | world | nebula
kinds = ["M"]       # S C M V P D icy | rocky icy ocean giant small | emission reflection dark planetary supernova_remnant
chance = 0.002
```

They appear only as deposits and markers, always flagged fictional: pinned on a body's
globe, listed in its system and in `geoforge show`, marked in the star view and found by
search. Each material draws from its own stream, so the list changes nothing else.

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

## The stars around a world

Its stars are its galaxy's, generated where they are looked at, the same from every world
of the galaxy. Positions are given in the world's own sky (RA and Dec in degrees, distance in
light years); a star by its identity or name:

```bash
geoforge stars near worlds/demo --radius-ly 20                 # the nearest systems, faint ones too
geoforge star add worlds/demo --name Beacon --ra 83 --dec -5 --distance-ly 40 --mass 3
geoforge star move worlds/demo Beacon --distance-ly 60
geoforge star cluster worlds/demo --name Hearth --ra 150 --dec -10 --distance-ly 400 --members 800 --age-gyr 0.05
geoforge constellation add worlds/demo Kite --line Beacon-Vaepemor --line Vaepemor-Baene
geoforge stars export worlds/demo --to stars.csv --limit 5      # edit in a spreadsheet, then
geoforge stars import worlds/demo stars.csv
geoforge visit worlds/demo Baene --to worlds/baene              # a new world at that star
```

Edits live in `galaxy/stars.json`, constellations in `constellations.json`; `sky` draws both.
A world without constellations of its own is shown figures drawn from its brightest stars;
`geoforge constellation auto worlds/demo` saves them, to rename (`remove`, `add`) or keep.

Beyond the galaxy, a universe of other galaxies is generated around it; the sky and the
report list those brighter than magnitude 8 and how they are seen:

```bash
geoforge universe worlds/demo --radius-mpc 5     # the home galaxy, its satellites and neighbours
```
A galaxy's model (components, arms, dust; see `geoforge_cosmos::galaxy`) can be edited:
put it as `galaxies/<galaxy>/galaxy.json` in the universe directory, and every body made
there after draws its stars from it.

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
geoforge report worlds/demo          # worlds/demo/report/index.html: the lite report, seconds
geoforge report worlds/demo --full   # every map, by layer, and the system's bodies
```

One page to read at a glance: the world whole on a globe in true colour, its key figures,
13 main maps with their keys (true colour, plates, geologic provinces, surface rock,
elevation, January and July winds, currents, temperature, precipitation, Köppen, rivers,
biomes) and its numbers against Earth's. `--full` adds every other map in themed groups
(hazards, pressure, seasons, sea ice, soils, resources, auroras) and the system. Places,
weather and the sky have their own commands (`place`, `almanac`, `weather`, `sky`). The
folder opens in a browser or publishes as is.

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
