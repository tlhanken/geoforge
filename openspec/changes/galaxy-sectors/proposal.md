# galaxy-sectors

## Why

The night sky is drawn afresh for each system from one smooth thin disk, in a random frame:
two neighbouring systems do not see each other or the same stars, no star in the sky can be
visited, the bright young stars come out evenly spread (about 6,500 naked-eye stars where a
clumpy real disk gives 9,110), and there is no bulge, bar, arm, thick disk, halo or dust
cloud. Longer term the cosmos must not be built around one planet: a universe should be
made from a seed and zoomed into, its details generated as they are looked at, while a
single system can still be made on its own. And worldbuilders must be able to move, add and
remove stars, for their own constellations and clusters.

## What Changes

- A galaxy made of components (nucleus, bulge and bar, thin disk by age, spiral arms, thick
  disk, halo, dust), each with its own shape, ages and metals, and morphologies built from
  them (spiral Sa–Sd with or without a bar, lenticular, elliptical, dwarf spheroidal, dwarf
  irregular, starburst). The Milky Way is one of them.
- Stars that live in the galaxy: space is cut into nested cubes, a cube's size set by how
  far its stars can be seen, so the bright stars of a region come from a few large cubes and
  its faint ones only from small cubes nearby; every star is generated from its cube's seed
  and has a stable identity, so any observer sees the same star, and zooming in generates
  only the cubes it reaches.
- Young stars born in groups (clusters and associations that spread with age), and dust in a
  smooth layer plus clouds.
- The night sky is read from these stars around the observer, and the galaxy's glow from
  its components through the dust.
- A system has a place in its galaxy: a standalone system draws one from its world's seed;
  a system made at a star takes that star (its mass, age, metals and companions).
- Edits: stars added, moved, changed or removed, clusters added, and an observer's own
  constellations, kept beside the generated stars and applied over them.

Not in this change: other galaxies in the night sky (Andromeda, satellite galaxies), the
universe, superclusters and galaxy rosters above one galaxy (the v2-architecture layers),
and automatic constellations. The nested cubes are meant to serve the universe's galaxies
the same way later.

## Capabilities

- `galaxy-model` — components, morphologies, density, ages and metals at a point
- `star-sectors` — nested cubes, stars generated per cube, identities, queries
- `night-sky` — the observer's stars and glow from the sectors
- `star-edits` — stars and clusters added, moved, changed and removed; constellations

## Impact

`geoforge-cosmos` (`galaxy`, a new `sectors`, `sky` rewritten on sectors, `population` and
`galactic` placing a system), `geoforge-core` (a system's galactic position, star
identities, edit types), `geoforge-store` (`galaxy/` in the world: the galaxy model and its
edits), `geoforge-almanac` (the sky reads the observer's position), `geoforge-cli` (star
and constellation commands). Skies change for every existing world.
