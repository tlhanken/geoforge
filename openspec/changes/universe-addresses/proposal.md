# universe-addresses

## Why

Every generated world sat in a galaxy of its own, called "Milky Way" whatever the seed, with
a universe seeded from that world alone: two worlds could never share a galaxy, and nothing
but the one focus body could be generated. The v2 architecture meant the opposite — a
universe from a seed, its galaxies, systems and bodies nested in it with the folder tree as
the address — and the viewer to come needs exactly that: zoom from the universe to any item
and generate it as it is looked at.

## What Changes

- **One address for everything** (`geoforge-core::address`): universe / galaxy / star's
  system / body, e.g. `7/h6:0,0,-1:17.0/f12:3,-4,0:5/4`; `real/milky-way/sun/5` in the real
  universe. Written, parsed, cut to a level, turned into directory names.
- **Generation from seed and address alone** (`geoforge-cosmos::resolve`): the galaxy at an
  address (a generated galaxy and a model built from it: the preset of its type scaled to
  its mass and size), the system born with a star, a body of that system. Nothing there
  knows about directories or defaults.
- **Defaults are an address chosen by a rule** (`resolve::complete`): a universe to its
  galaxy most like the Milky Way nearest its origin; a galaxy to the nearest Sun-like star
  with a wet world near its reference circle; a system to its most Earth-like body.
- **A universe directory** (`geoforge-store::UniverseDir`): `universe.json`, shared galaxy
  documents (`galaxies/<galaxy>/` with the model if edited and the star edits) and a world
  per body made (`bodies/<galaxy>/<star>/<body>/`). Commands given the directory work on
  its focus body.
- **CLI**: `new --seed N` makes universe N's default body; `new --at <address>` any other
  (added to an existing universe directory); `focus` lists or switches the bodies made;
  `visit` goes to a star's system in the same universe. `--habitable`, `--galaxy` and
  `--galaxy-file` are gone from `new` (a galaxy's model is edited in the universe instead);
  `--archetype` and `--body` still make a world on its own.
- The sky and `universe` find the home galaxy by the address.

Then, in this change's later parts: a true-colour appearance view for every body and the
viewer (zoom universe → galaxy → system → body, search, layers, sphere or map).

## Capabilities

### New Capabilities
- `viewer`: browse a universe directory, generating what is looked at.

### Modified Capabilities
- `universe`: galaxies are reached by address and have models of their own.

## Impact

- `geoforge-core`: `Address`, `GalaxyRef`, `StarRef`, `Level` (replacing the unused
  `key::Address`); `SystemRoot.placement` is the system's address.
- `geoforge-cosmos`: `resolve`; `UniverseGalaxy::model`, `stars_seed`; `Galaxy::scaled`,
  `Shape::scaled`; `Universe::home` matches the type exactly first, `home_at`.
- `geoforge-store`: `UniverseDir`; `Manifest.address`; `World.universe`; galaxy documents
  through `galaxy_dir` (`load_galaxy_model`, `save_galaxy_model`).
- `geoforge-pipeline`: `create_at`. `create_world` stays for worlds on their own (and tests).
- Worlds made before keep working as worlds on their own; their galaxy is still named for
  the preset they were drawn from.
