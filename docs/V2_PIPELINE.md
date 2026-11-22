# Geoforge v2

A world is a directory. Seeds generate every layer; any layer can be exported to files that
mapping or paint software opens, edited or drawn from scratch, and imported back. Everything
downstream then regenerates from your version.

```text
cosmology (thin: one solar system)  →  BodyContext  →  surface layers  →  (refinement, later)
                  tectonics → geology → elevation → circulation → climate → erosion → hydrology → biomes
```

Any solid body can be a world: Earth-like planets, Mars-, Venus- and Mercury-like worlds,
moons, frozen ice bodies, asteroids. Each layer is a *slot*: the body's tectonic regime and
what it has (air, surface water) decide which generator fills it, or that it does not apply.

## Quick start

```bash
cargo build --release -p geoforge-cli          # binary: target/release/geoforge-v2

geoforge-v2 new worlds/demo --seed 42 --habitable   # a Sun-like star and its Earth-like world
geoforge-v2 generate worlds/demo --render      # ~40 s at 2048×1024; previews in worlds/demo/renders
geoforge-v2 stats worlds/demo                  # summary numbers next to Earth's
```

Other starting points:

- `new --seed 42` without `--habitable` draws the system from the Milky Way's stellar
  population: three in four are red dwarfs, one in sixteen a white dwarf, a few in a thousand
  a giant, a blue star or a neutron star. The focus is the system's most Earth-like solid
  body (often a tidally locked planet of a red dwarf); about one seed in six has nothing
  solid (a lone black hole, say) and `new` says so. `--galaxy elliptical|dwarf|starburst`
  draws from another population. `census` counts thousands of seeds and lists example
  seeds of each kind of star.

- `new --archetype mars` starts from a real body instead of a generated system: `earth`,
  `mars`, `venus`, `mercury`, `moon`, `titan`, `europa`, `gas_giant`, `ice_giant` (real
  sizes, orbits, shells and air).
- `new --seed 42 --focus "Kibedo VII b"` generates the system and focuses on another of its
  bodies: a planet or moon by name (`show` lists them with what each is), or
  `"<belt name> #<n>"` for the n-th largest member of a belt (generated on demand).
- `new --body my_world.toml` uses a body you wrote (a `BodyContext` in TOML or JSON;
  `bodies/<name>/body.json` of any world is a template).

## Commands

| Command | What it does |
|---|---|
| `new <dir> --seed N [--habitable \| --galaxy NAME] [--width 2048] [--archetype NAME \| --focus BODY \| --body file]` | Create a world |
| `generate <dir> [--render] [--from LAYER] [--until LAYER] [--force]` | Build what is missing or stale |
| `plan <dir>` | Which generator fills each layer, what `generate` would run, and why |
| `status <dir>` | Source (generate / patched / authored) and state of every product |
| `export <dir> [products\|layers\|views\|all] [--to DIR]` | Write edit views (default `<dir>/edit`); derived views as PNG + legend + world file |
| `import <dir> <product> [--from DIR] [--authored]` | Bring an edit back as a patch, or as authored data |
| `reset <dir> <product>` | Drop a product's edits |
| `render <dir> [products\|layers\|views] [--to DIR]` | Preview PNGs (all by default) |
| `show <dir>` | The solar system (each body with its regime, air and water) and the focus: circulation, magnetic field, breathability, tides, calendar |
| `positions <dir> [--day T] [--to FILE.json]` | Where every star, planet and moon is T Earth days from the present (Kepler orbits), the focus's local date, and how lit its moons are |
| `seed <dir> <product> [--owner KEY]` | The seed a product is generated from |
| `stats <dir>` | Realism summary with Earth reference values |
| `census [--galaxy NAME] [--count N] [--from SEED]` | Draw many systems: kinds of star, planets, wet worlds, example seeds |

## The world directory

```text
world.toml                       seed, grid, layer parameters, product sources, provenance
system.json                      the solar system (absent for a hand-written body)
bodies/<body>/body.json          BodyContext: the only thing the surface layers see
bodies/<body>/products/          generated products (cache: always rebuildable)
bodies/<body>/edits/             your authored data and patches (never overwritten)
cache/baselines/                 what each edit view showed when exported (for diffs)
cache/bases/                     generated product under each patch
edit/                            default edit-view folder
renders/                         preview images
```

Parameters of the generators this body uses are written out in `world.toml`
(`[params.<generator>]`) with their defaults, so they can be edited. Changing one makes
exactly that layer and everything after it stale.

## Layers and generators

| Layer | Plate-tectonic worlds | Rigid-lid worlds (stagnant, episodic, inactive) | Needs |
|---|---|---|---|
| tectonics | `tectonics` | `lid_tectonics` | rocky crust |
| geology | `geology` | `lid_geology` | rocky crust |
| elevation | `elevation` | `lid_elevation` | rocky crust |
| circulation | `circulation` | `circulation` | solid surface, atmosphere |
| climate | `climate` | `climate` | solid surface, atmosphere |
| erosion | `erosion` | `erosion` (copies the rough heightmap: no water to cut with) | rocky crust |
| hydrology | `hydrology` | `hydrology` | liquid water on the surface |
| biomes | `biomes` | `biomes` | liquid water on the surface |
| bands | — | — | no solid surface, an atmosphere: `giant_bands` |

So Earth runs everything; Mars and Venus get a surface, winds and a climate but no rivers
or biomes; the Moon, Mercury, frozen moons and asteroids get a surface only. Giants get
only the `bands` layer. Ice-shell worlds (Europa) have no generator yet: every layer
reports why. A generator for a
new kind of body is a new `Generator` registered in `geoforge-pipeline/src/registry.rs`.

## Products

| Layer | Products | Shape |
|---|---|---|
| tectonics | `plates`, `crust_type`, `crust_age`, `tectonic_events` | partition, field, field, timeline |
| geology | `provinces`, `geology_events`, `lithology` (rock catalog), `rock_age` (Ma) | partition, timeline, catalog partition, field |
| elevation | `rough_elevation`: the heightmap from tectonics, before rivers cut it (m above sea level; above the mean radius on dry worlds) | field |
| circulation | `wind_jan`, `wind_jul` (m/s), `ocean_currents` (m/s), `sst_anomaly` (°C) | vector fields, field |
| climate | `temperature` (°C), `temperature_range` (warmest − coldest month, °C), `precipitation` (mm/yr), `precipitation_warm_share` (share in the warm half-year) | fields |
| erosion | `elevation`: the final heightmap, valleys cut where the rain runs off (m) | field |
| hydrology | `rivers`, `lakes`, `watersheds` | network, partition, partition |
| biomes | `biome` (biome catalog) | catalog partition |
| bands | `zonal_wind` (m/s, eastward), `cloud_bands` (cloud-top brightness, 0 darkest belt – 1 brightest zone or oval) | fields |

Derived values (plate boundaries, slope, sea level, flow directions) are recomputed, not stored.

**Derived views** are maps that are pure functions of stored products. They are worked out
when rendered or exported, never stored, so they cannot go stale and are not edited (edit
the products they read instead):

| View | Reads | Export |
|---|---|---|
| `koppen` | `temperature`, `temperature_range`, `precipitation`, `precipitation_warm_share` | `koppen.png` (exact class colours, Beck et al. 2018's palette), `koppen.legend.json` (code, name, colour, share of land), `.pgw`, `.prj` |

The Köppen–Geiger view has 30 classes by Peel et al. (2007)'s rules: arid (B) first, then
polar, tropical, and temperate or continental at a 0 °C coldest month
(`KoppenRules::cold_boundary_c`; Köppen's −3 °C is the alternative). Monthly figures are
rebuilt from the stored ones, with temperature and rain as sinusoids in phase. The rain
curve's amplitude gives the stored warm-half share.

## Editing a layer

```bash
geoforge-v2 export worlds/demo rough_elevation   # edit/rough_elevation.png (+ .pgw, .prj, .edit.json)
#   … paint in GIMP, Krita, QGIS, Photoshop …
geoforge-v2 import worlds/demo rough_elevation   # a patch: only the pixels you changed
geoforge-v2 generate worlds/demo --render        # climate, valleys, rivers and biomes follow
```

**Which heightmap to edit.** The climate sees `rough_elevation`; `erosion` then cuts valleys
into it where the rain runs off and writes `elevation`. To move mountains or coasts, edit
`rough_elevation`. A finished heightmap drawn elsewhere goes in as `rough_elevation
--authored` with `[params.erosion] k = 0`, so it passes through unchanged and the climate
uses it. Edits to `elevation` itself (a lake basin, a gorge) reach rivers and biomes but not
the climate.

Edit views:

| Product kind | View |
|---|---|
| Continuous field | 16-bit greyscale PNG; `value = offset + scale × pixel` (range in the sidecar; elevation −12 000 … 12 000 m) |
| Vector field (winds, currents) | 16-bit RGB PNG: red = eastward, green = northward, each `offset + scale × channel` |
| Category field (crust type) | RGB PNG in the category palette; off-palette pixels (anti-aliasing) snap to the nearest colour |
| Partition (plates, provinces, lakes, watersheds) | RGB PNG, one colour per key, black = none, plus `<p>.table.json` (key, colour, row). Add a region by adding a row with a new key and colour |
| Catalog (lithology, biome) | As a partition, but each row is a kind (rock, biome) painted in its own colour, and every row is kept even if unused. Add a kind by adding a row (key, colour, slug, name, …) and painting with its colour |
| Timeline, network | The JSON document |

Every raster has a `.pgw` world file and a `.prj` sphere definition (the body's radius), so GIS
tools place it on the globe. The grid is equirectangular with pixel centres at
`lon = −180 + (x + ½)·360/W`, `lat = 90 − (y + ½)·180/H`.

**Patch or authored.** By default an import is diffed against the view as it was exported and
stored as a patch: deltas for continuous fields, new values for categories and owners. If the
layer is later regenerated underneath (a parameter changed, an upstream edit), the patch is
re-applied on top and `status` shows it as *rebased*; edits that no longer fit (a cell outside
the grid, a key with no row) are counted as *orphaned*. `--authored` replaces the product
outright: draw a whole map and the generator for that product is skipped. A new plate key
without a row gets a seed point at its centroid and a seeded motion; a new lake gets its
surface at the highest cell.

## Rocks and biomes are catalogs

Rock types and biomes are data, not code: Earth's are built in, and a world adds its own
(fantasy or otherwise) in `world.toml`:

```toml
[[params.geology.extra_rocks]]
slug = "starmetal"
name = "Starmetal"
color = "#c8d2ff"
erodibility = 0.1            # 1 = soft sediment, ~0.1 = very hard
tags = ["ore", "magical"]    # free-form, for later layers
age = "volcanic"             # crust | basement | archean | plutonic | volcanic | flood | ophiolite | infill | cover
occurs = [{ province = "volcanic_arc", crust = "continental", weight = 0.08 }]

[[params.biomes.extra_biomes]]
slug = "glowcap_forest"
name = "Glowcap Forest"
color = "#5a3a8a"
within = ["boreal_forest", "temperate_rainforest"]   # Earth biomes it may replace; empty = any land
temperature_c = [-5.0, 12.0]                          # optional bounds: temperature_c, precipitation_mm, elevation_m
cover = 0.3                                           # share of that ground it takes, in patches
```

A rock's `weight` is added to the mix of each matching province (built-in mixes sum to 1).
Extra biomes are checked before Earth's classification, in order. On a rigid-lid world use
`[[params.lid_geology.extra_rocks]]`. Painted edits can use any row of the catalog.

## Provenance

Each stored product records a hash of everything it was built from: the generator and its
version, the root seed, the body, the grid, the generator's parameters, its products'
sources and edit files, and the provenance of every product the layer reads. `generate`
reuses a product whose recorded hash matches and rebuilds the rest, so an edit to elevation
reruns elevation, climate, hydrology and biomes, and nothing else. Bump a generator's
`VERSION` when its algorithm changes output; until then use `generate --from <layer>` while
iterating on code.

## Seeds

A product's seed is `root_seed → product name → owner key`, never an ancestor, an index or a
name. Adding a layer, renaming a body or editing an unrelated product never changes another
product's randomness. `seed <dir> <product>` prints it.

## Models (short)

- **Solar system, drawn from a galaxy** (the default): the galaxy is a population model
  (Milky Way: Kroupa 2001 initial mass function; stars forming at a constant rate for
  11 Gyr; [Fe/H] solar for young stars at the Sun's radius, −0.02 dex per Gyr of age,
  −0.06 dex/kpc outward, 0.2 dex scatter; an exponential disk of 2.6 kpc). A draw gives a
  site, an age, a metallicity and stars at birth: singles, binaries and triples by primary
  mass (25% of red dwarfs, 46% of Sun-like stars and 90% of O stars have companions;
  log-normal periods, flat mass ratios). Each star is evolved to the present: main
  sequence brightening (the Sun 0.72 → 1 → 2 L☉), subgiant, red giant or supergiant, then
  a white dwarf (Cummings et al. 2018 masses, Mestel cooling), a neutron star (8–20 M☉) or
  a black hole (over 20 M☉); a supernova that loses half a pair's mass unbinds it. Planets
  form on a ladder from a few days out, at Kepler-era rates by host mass and metallicity
  (giants ∝ M·10^(1.2[Fe/H]), Johnson et al. 2010; many small planets around red dwarfs,
  Dressing & Charbonneau 2015), stunted by close companions (Moe & Kratter 2021) and kept
  to stable orbits (Holman & Wiegert 1999); giants swallow what is within three of their
  radii, white dwarfs keep only what formed beyond ~2 AU, on widened orbits, neutron stars
  only rare second-generation planets (PSR B1257+12), black holes nothing. Planets within
  0.35 AU·M^(1/3) of their star are tidally locked. Neighbours are kept 8 mutual Hill radii
  apart; closer than 5 they scatter (the lighter thrown out, the survivor eccentric). The
  asteroid belt lies between the first giant's 4:1 and 2:1 resonances, the Kuiper belt
  between the outermost giant's 3:2 and 2:1; moons between the Roche limit and 0.4 of the
  Hill radius. The focus is the most Earth-like solid
  body (surface water first, then plates and air, then starlight and mass near Earth's).
  Against the Milky Way near the Sun: 75% red dwarfs, 6% white dwarfs, 0.5% neutron stars,
  0.2% black holes, 25% of red dwarfs and 46% of Sun-like stars multiple.
- **Solar system, habitable recipe** (`--habitable`): an F8–K5 star; planets on a geometric
  ladder through the habitable zone with the focus world in it, settled by the same spacing rules (the focus stays put); belts at the giants' resonances; 2–5 large moons per gas giant, 0–2 per ice
  giant, sometimes a large moon for the focus. What each body is follows from its mass,
  starlight and where it formed: airless rock (Mercury, the Moon), runaway greenhouse (Venus:
  ~90 bar CO₂, episodic lid), ocean and plate tectonics (Earth), thin CO₂ (Mars); icy bodies
  with an ice shell over an ocean when tidally heated (Europa), sometimes N₂–CH₄ air (Titan),
  else frozen and cratered. On Earth-like worlds the carbonate–silicate cycle sets CO₂
  where the world weathers as fast as Earth (Walker, Hays & Kasting 1981): dim worlds
  build up CO₂ and run a few degrees cooler, bright ones draw it down.
- **Positions and calendar** (`geoforge_core::ephemeris`): each node on a fixed Kepler
  orbit about its parent from its elements (no mutual perturbations); a pair's stars turn
  about their barycentre, opposite each other. Solar day 1/(1/spin − 1/year) (+ for a
  retrograde spinner: Venus 116.75 days), none for a world locked to its star; months are
  synodic, 1/(1/orbit − 1/year).
- **Magnetic field** (`geoforge_core::magnetism`, printed by `show`): a rocky world's
  dynamo lasts 10 Gyr·√M with plates, 1.5 Gyr·√M under a stagnant lid, 1 Gyr·√M under an
  episodic one (Mars's stopped ~4.1 Gyr ago); giants always have one. Surface field 30
  µT·M^(1/3) (Jupiter 420 µT by mass, ice giants ~20 µT multipolar and tilted 40–60°); a day
  over 200 h makes it multipolar and weaker. The star's wind (2 nPa at Earth, ∝ 1/d², stronger
  around young stars) sets the magnetopause by pressure balance (Earth 9.5 radii) and the
  auroral ovals (Earth 71° magnetic latitude). Derived, not stored; the hazards layer will
  read it.
- **Tectonics**: v1's electrostatic plate seeding (power-law charges), now with each cell
  assigned by distance ÷ charge^0.22 so a seed's charge sets its plate's size
  (`size_weighting = 0` gives v1's plates exactly); domain warp and island removal;
  Euler-pole motion, per-cell boundary classification, crust age from ridge distance.
  Continents are embedded in plates with ocean floor to their passive margins; the largest
  plate and about 30% of the rest are oceanic.
- **Geology**: v1's 20 province types (collision belts now swell and pinch along their
  length), plus lithology: each province kind has a rock mix (craton: basement, greenstone,
  granite; platform: sandstone, limestone, salt; orogen: metamorphic core, granite,
  deformed sediment, ophiolite; arc: volcanics and plutons…) picked by smooth noise, with
  rock ages from the province's history.
- **Elevation**: GDH1 ocean depth from crust age, continental freeboard with plateaus,
  province relief (continental collisions: foreland basin, front range, plateau on the
  upper plate; heights × 1/g), shelves and slopes planed to sea level, sea level from the
  ocean's volume. This is the rough heightmap.
- **Erosion** (after the climate): stream-power erosion (Braun & Willett implicit scheme)
  with erodibility from lithology and uplift in active belts, each cell's drainage area
  weighted by its runoff (rain minus evaporation, Budyko) to the power 0.5. Wet ranges are
  cut deep, deserts barely, so dry basins stay closed. Valley floors cut below the rough
  surface are warmer by the lapse rate when hydrology and biomes read the temperature.
- **Rigid lids** (`lid_*`): one plate; crust split into old highlands and younger lowlands
  (one hemisphere on stagnant lids like Mars, scattered on Venus, none on dead bodies);
  provinces for cratered highlands or ice, volcanic plains, tesserae, impact basins (the
  large end of a Neukum crater population, denser on older ground), volcanic rises with
  shield volcanoes (only on bodies big enough for long-lived plumes), rifts and coronae;
  elevation from the step, each landform and smaller craters, with relief scaled by
  gravity (Olympus Mons is 2.5× Mauna Loa at 0.38 g).
- **Circulation**: January and July winds. Cells by rotation and size: a Hadley cell of
  width ∝ Ω^(−2/3) (30° on Earth), then eddy-driven jets spaced by the Rhines scale,
  round((a·Ω / a⊕·Ω⊕)^0.75) per hemisphere (Chemke & Kaspi 2015; GFDL: about three at 4×
  Earth's rotation), each with westerlies, a storm track and easterlies poleward: 1 + 2·jets
  cells, so Earth and Mars 3, a 12-hour Earth 5, a 6-hour one 7, Venus and Titan a Hadley
  cell to a polar cap (`show` prints the count). Tidally locked worlds (day = year) instead
  circulate around the substellar point (0°N 0°E): air rises under the star and returns along
  the ground from the night side, with an equatorial jet below ~20-day periods (Haqq-Misra
  et al. 2018). Then the seasonal ITCZ, monsoon, subtropical highs over each ocean basin); wind-driven ocean gyres
  (Sverdrup interior, Stommel western boundary currents, Ekman drift, equatorial currents,
  circumpolar current); SST anomaly from currents crossing the temperature gradient, carried
  downstream, plus coastal upwelling.
- **Bands** (giants): jets per hemisphere from the Rhines scale with ~60 m/s winds on gas
  giants and ~250 on ice giants (Jupiter 6, Saturn 5, Uranus and Neptune 1); gas giants
  superrotate at the equator, ice giants run retrograde there with a jet near 70°; zones
  bright where the shear is anticyclonic, belts dark where cyclonic, contrast fading toward
  the poles; wavy band edges, white ovals in the zones and, on half of giants, one great
  storm in the widest subtropical zone (a dark red spot on gas giants, a dark spot on ice
  giants; `great_storm_chance`).
- **Climate**: zonal energy balance on the body's own atmosphere (on a locked world, by
  angle from the substellar point, with no seasons, a convective cloud deck under the star
  on wet worlds (Yang et al. 2013) and 1.6× the heat transport: Earth locked would be ~30 °C
  under the star and ~−40 °C opposite): a gray greenhouse
  (OLR = σT⁴/(1 + ¾τ), τ from CO₂, collision-induced absorption and water vapour; fitted to
  Earth, Mars, Venus and Titan), heat transport and seasonal damping from the air column's
  mass, frost albedo from water snow or CO₂ frost; seasons from obliquity, eccentricity and
  the year's length; continentality and current-warmed air
  carried by the winds; water vapour carried by each season's winds, raining out where air
  rises or cools (storm tracks follow the seasons), recycled over land; global
  precipitation from a water budget; lapse rate and relief on the fine grid.
- **Hydrology**: Budyko (Fu) runoff; priority-flood drainage with jittered routing across
  flats; lakes from each depression's water balance (open where overflowing and kept open by
  rifting, depth or glaciation; endorheic and shrunk to balance where not); arid channel
  losses; Strahler-ordered rivers; watersheds by terminal outlet.
- **Biomes**: Köppen thresholds (10 °C tree line, 18 °C tropics, aridity threshold that
  depends on temperature and rain seasonality, dry-summer Mediterranean) and Whittaker's
  rainfall steps, alpine by altitude, wetlands on flat floodplains.

## Realism status

See `docs/V2_REALISM.md` for the numbers across seeds and the known gaps.
