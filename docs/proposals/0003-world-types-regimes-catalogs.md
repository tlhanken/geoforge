# World types, tectonic regimes and catalogs

- Proposal #: 0003
- Start Date: 2026-11-19
- Status: implemented
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

### Climate

The column physics is the body's: `geoforge-core::radiation` holds the gray atmosphere,
`σT⁴ / (1 + ¾τ)`, with τ from CO₂ (pressure-broadened), collision-induced absorption
(∝ P²/g, strongest for cold bodies) and water vapour where there is surface water, fitted to
Earth (OLR 231 W/m² and slope 2.09 W/m²/K at 15 °C), Mars, Venus and Titan; and the open
albedo (clouds, thick cloud decks). The spatial model is `geoforge-surface::atmosphere`
(`grid`, `transport`, `cells`, `radiation`), shared by circulation and climate; `climate`
keeps the generator and `moisture`. Heat transport and seasonal damping follow the air
column's mass and the year's length; frost albedo comes from water snow or, under thin CO₂
air, CO₂ frost. The carbonate–silicate thermostat moved to the solar system step, which sets
the CO₂ of plate-tectonic worlds with liquid water. Dry worlds have no sea and no rain;
hydrology and biomes require surface liquid water.

### Catalogs

`lithology` and `biome` become partitions whose tables are catalogs (`ProductId::is_catalog`):
rows are kinds, each with its own colour, kept even when unused.

- Product rows: `Rock { slug, name, color, erodibility, tags }` and
  `BiomeClass { slug, name, color, tags }` (core).
- Rocks: built-in rocks keep stable keys with their per-province mixes; `RockDef` (a `Rock`
  plus an `AgeRule` and `occurs` weights per province kind and crust) adds rocks through
  `extra_rocks`. Weights add to each matching mix.
- Biomes: the validated Köppen/Whittaker classifier stays code over built-in `Earth`
  biomes; `BiomeDef` (a `BiomeClass` plus an `Envelope`: temperature, precipitation and
  elevation bounds, `within` built-in biomes, patchy `cover`) adds biomes through
  `extra_biomes`, checked before the classifier, in order.

Both reuse the partition edit view; the palette is each row's colour, and a new row is a new
kind.

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
