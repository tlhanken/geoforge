# World types, tectonic regimes and catalogs

- Proposal #: 0003
- Start Date: 2026-11-19
- Status: accepted
- Origin: emerged in discussion (v2 cutover)

## Summary

Make every kind of body a first-class world: Earth-like, Mars-like, Venus-like, airless
rocks, icy moons, giants and asteroids. A body carries a **tectonic regime**; each pipeline
layer becomes a **slot** that several generators can fill, and the first whose
preconditions hold for the body runs. Climate takes its greenhouse from the atmosphere's
composition and pressure. Rock types and biomes become **catalogs** (data tables) so fantasy
entries are rows, not code.

## Motivation

Any rocky surface currently gets plate tectonics, the solar system only makes Earth-like or
bare-rock worlds, climate is fitted to Earth, and rocks and biomes are closed enums.

## Design

### Body (`geoforge-core::body`)

Split into modules, each one concern:

| Module | Contents |
|---|---|
| `system` | `SolarSystem`, `OrbitalNode`, `NodeKind`, `Star`, `BodyExterior`, `BodyClass`, `Elements`, `Population`, `Epoch`, `SystemRoot` |
| `matter` | `Species` (with molar mass), `Composition`, `Phase` |
| `interior` | `ShellKind`, `Shell`, `BodyInterior`, `TectonicRegime`; queries: solid surface, lithosphere, atmosphere and partial pressures, surface liquid water |
| `context` | `BodyContext`, `Irradiance`, `TidalPartner` |
| `presets` | `Archetype` (Earth, Mars, Venus, Mercury, Moon, Titan, Europa, gas giant, ice giant) → `BodyContext` with real values |

`TectonicRegime`: `Plates` (Earth), `StagnantLid` (Mars, Mercury, Moon), `Episodic` (Venus),
`IceShell` (Europa, Titan), `Inactive` (giants, rubble). It lives on `BodyInterior`; the
solar system step infers it from mass, water and heat; a hand-written body states it.

### Generator slots (`geoforge-core::product`, `geoforge-pipeline::layers`)

`Generator` gains `const NAME`. A layer is a slot; the pipeline's `Implementation` enum lists
every generator with its layer. For a body, `select(layer)` returns the first implementation
whose `applies(body)` holds. Parameters (`[params.<name>]`) and provenance are keyed by
implementation name, so switching regime rebuilds exactly what changes.

| Layer | Plates | Stagnant lid / episodic |
|---|---|---|
| tectonics | `tectonics` | `lid_tectonics` (one plate, crustal dichotomy, surface age) |
| geology | `geology` | `lid_geology` (impact basins, volcanic rises, shields, rifts, plains, coronae, tesserae) |
| elevation | `elevation` | `lid_elevation` (craters, shields scaled by gravity, dichotomy) |
| circulation, climate, hydrology, biomes | shared | shared (preconditions decide) |

All implementations of a slot produce the same products, so downstream layers do not care.

### Climate (`geoforge-surface::climate`)

Split into `radiation` (energy balance), `grid` (atmosphere grid and transport), `cells`
(Hadley/Ferrel/polar winds), `moisture` (water vapour and rain) and the generator.
Outgoing longwave is a gray atmosphere, `σT⁴ / (1 + ¾τ)`, with τ from partial pressures
(CO₂, CH₄, N₂ collision-induced) plus water vapour feedback where there is surface water;
calibrated so Earth matches the previous fit, Mars ≈ −60 °C and Venus ≈ 460 °C. Base albedo,
ice-albedo strength, heat transport and seasonal damping follow the atmosphere. The
carbonate–silicate thermostat applies only to plate-tectonic worlds with liquid water.
Hydrology and biomes require surface liquid water.

### Catalogs

`lithology` and `biome` become partitions whose tables are catalogs:

- `Rock { id, name, color, erodibility, age_rule, occurs: {province kind → weight}, tags }`.
  Built-in Earth rocks; extra rocks in `[params.geology] extra_rocks`. Geology picks rocks by
  each province kind's weights, so a fantasy rock with weights appears on its own.
- `BiomeDef { id, name, color, priority, rules: [Envelope] }`, where an envelope bounds mean,
  warmest and coldest month temperature, range, precipitation, aridity, warm-season share,
  elevation and setting (land, sea, lake, floodplain). The Köppen/Whittaker classification
  is expressed as built-in rows; fantasy biomes add rows with a priority, or are painted.

Both reuse the partition edit view (palette from the table; add a row to add a kind).

### Solar system (`geoforge-cosmos`)

Atmosphere, hydrosphere and regime are inferred from mass, starlight and water: runaway
greenhouse worlds become Venus-like, small or dry ones stagnant-lid, cold ones keep thin
CO₂ or none, big cold moons can hold N₂/CH₄ (Titan), icy moons get ice shells over oceans.
Giants get moon systems. Belts stay populations; individual asteroids are generated on demand
by index. The CLI can start from an archetype or focus any body of the system.

## Alternatives Considered

- A boolean "tectonically active" flag: too coarse; Venus and Mars differ in kind.
- One generator per layer branching on regime inside: mixes unrelated models in one file and
  hides which model ran.
- Fantasy rocks as code: every addition needs a rebuild and a release.

## Open Questions

- Ice-shell surfaces (Europa) get a slot entry later; until then their surface layers report
  "no generator for this regime".
- Non-water hydrology (Titan's methane) is out of scope.
