# Design: universe-addresses

## Addresses
`Address { universe: UniverseSite, galaxy: Option<GalaxyRef>, star: Option<StarRef>, body:
Option<Key> }`, written `/`-separated: the universe's seed (or `real`), the galaxy
(`GalaxyId`, or `milky-way`), the star (`StarId`, or `sun`), the body's key. A star needs a
galaxy and a body a star; the Milky Way and the Sun only in the real universe, generated
galaxies only in seeded ones. `level`, `parent`, `up_to(level)`, `segments()` (`:` → `_`
for directory names). Serialized as its text.

## Resolving (`geoforge_cosmos::resolve`)
- `galaxy(address) → GalaxyAt { generated, model, stars_seed }`. A generated galaxy's
  model: its type's preset (the Milky Way's for a barred Sbc) scaled by
  `Galaxy::scaled(mass, length)` — masses ×mass; every length (shapes, heights, arms, dust,
  reference and edge radii) ×length, length = its half-light radius over the half-light
  radius a galaxy of the preset's type and mass has (Shen et al. 2003); metallicity shifted
  0.35 dex per dex of mass (Gallazzi et al. 2005, −0.6 to +0.15); the dust's density by
  mass/length³ (×¼ to ×4). Stars seeded from the galaxy's seed. The Milky Way: its preset,
  stars seeded from `REAL_SEED`.
- `system(address, field)`: the system born with the star (`galactic::system_at`), or the
  real Solar System at the 2026 equinox; its `placement` is the address and its site names
  the galaxy and universe.
- `body(address, system)`: the body's context; the real Solar System's measured ones.
- `complete(address)`: a universe → its galaxy most like the Milky Way (exact type, within
  0.25 dex of stellar mass, nearest the origin; then the same class, then any); the real
  one → the Sun. A galaxy → from a point on its reference circle (azimuth from its stars'
  seed), searching 30, 60, 120, 240 pc, the nearest single main-sequence star of
  5,200–6,200 K whose focus body has surface water; else the nearest star with a wet focus;
  else any with a solid body. A system → its focus. Universe 7 resolves in about a second.

## Universe directory (`geoforge_store::UniverseDir`)
```text
universe.json                  { universe, focus }
galaxies/<galaxy>/             galaxy.json (if edited), stars.json (star edits)
bodies/<galaxy>/<star>/<body>/ a world: world.toml (with its address), system.json, bodies/…
```
`World::open` on a universe directory opens its focus body's world; a world with an
address finds its universe above it, and its galaxy documents are the universe's. Only what
is made or edited is written. A body's layers are seeded from its system's seed (63 bits).

## Commands
`new DIR [--seed N] [--at ADDRESS] [--focus NAME] [--width W]`; in an existing universe
directory the universe part of the address may be left out. `focus DIR [ADDRESS]`.
`visit DIR STAR` makes the star's system's focus in the same universe (`--to` for worlds on
their own).

## Later in this change
- Appearance: a true-colour view for every body (the focus world from its products; other
  bodies from their bulk properties; stars by temperature), low resolution, cached.
- Viewer: `geoforge view DIR` — a local server (Rust: resolving, search, layers, images)
  and a thin web page (three.js) — zoom universe → galaxy → system → body, search by name
  among what is generated nearby or made, and by address anywhere; layers on a sphere or a
  map.
