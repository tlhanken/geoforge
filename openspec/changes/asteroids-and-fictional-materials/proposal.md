# asteroids-and-fictional-materials

## Why

A belt's members were rock or ice of one density and nothing more: no spectral types, no
elements, no families, so a belt could not be told from another or mined in a story. There
was no way to put a setting's fictional materials into a universe without editing the
generators.

## What Changes

- **Asteroids made of something** (`geoforge_cosmos::asteroids`): each belt member draws a
  spectral type (S, C, M, V, P, D; icy in an icy belt) from the Solar System belt's mass
  shares by distance against the frost line (DeMeo & Carry 2014); its elements are its
  meteorite analog's (CI and L chondrites, irons, eucrites, Tagish Lake; 67P's mix of ices,
  organics and dust for icy bodies) with a ±15% spread; its porosity follows its size
  (rubble piles 20–50%, comets 40–75%, large bodies and monoliths little) and its density
  the grain density less the voids; a third of members belong to one of 5–30 named families
  that share a type and, within 0.01 AU and a few hundredths, an orbit. Free metal and
  platinum-group mass are told.
- **Belts' mixes** in the viewer's system panel and the CLI's `show`; belt members resolve
  by address (`<system>/<belt key << 32 | rank>`).
- **Fictional materials** (`geoforge_core::exotics`, `geoforge_cosmos::exotics`): a
  universe's `exotics.toml` (a commented template is written with every new universe) lists
  materials and where they are found — in asteroids by type, on worlds by kind (rocky, icy,
  ocean, giant, small), in nebulae by kind — optionally only round given stages of star,
  each with a chance per candidate. Each material draws from its own stream keyed by its
  name, so the list moves nothing else. They appear only as deposits and markers, flagged
  fictional: pinned on a body's globe and map with their tonnage, listed in a system (belt
  bodies expected, with examples), diamonded in the star view, found by search, printed by
  `geoforge show`. A faulty list is named, not fatal.

## Impact

- `geoforge-core`: `exotics` (new: the list's types and checks).
- `geoforge-store`: `UniverseDir::exotics`, the template on `create`.
- `geoforge-cosmos`: `asteroids` (new), `exotics` (new), `context` (belt members made by
  type; `belt_member_made`, `frost_line_of`), `resolve` (belt members by address).
- `geoforge-pipeline`: the report shows an asteroid focus's type and makeup.
- `geoforge-viewer`: body facts and deposit pins, system finds and belt mixes,
  `/api/exotics`, search.
- Belt members' densities, masses and (for family members) orbits change.
