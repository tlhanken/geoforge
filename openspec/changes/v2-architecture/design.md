# Design — v2-architecture

## Pipeline framing

The pipeline descends a chain of nested containers. Each layer produces the **interior** of
the container created by the layer above it. A layer never decides another layer's children.

    container above     layer              produces (the interior)
    ─────────────────────────────────────────────────────────────────────────────
    (root)              universe           the universal constants
    the universe        superclusters      the superclusters          (collection)
    a supercluster      galaxies           the galaxies               (collection)
    a galaxy            solar_systems      the solar systems          (collection)
    a solar system      solar_system       stars, orbits, moons, belts (structure)
    a planetary body    planetary_body     atmosphere, hydrosphere, geosphere (structure)

Two kinds of interior: a **collection** of child entities (plural layer names) or a
**structure** of fields of one entity (singular layer names).

Canonical roster (16 layers, generation order):

    universe superclusters galaxies solar_systems solar_system planetary_body
    tectonics geology elevation climate hydrology biomes
    resources hazards settlements kingdoms
    [+ culture, religion — reserved]

`planetary_body` is the generic non-stellar body: terrestrial planet, gas giant, ice giant,
dwarf planet, moon, asteroid, comet. Its interior (atmo/hydro/geo) is kind-dependent. The
surface pipeline may target any body kind, not only a rocky world.

## World container

A world is a directory; names are unique within a parent, not globally, so the folder tree is
the address:

    world/
      world.toml
      cosmology/
        universe.toml
        superclusters/<c>/
          cluster.json
          galaxies/<g>/
            galaxy.json
            solar_systems/<s>/
              system.json
              bodies/<b>/
                body.json
                layers/<field>.png + <field>.json
                entities/settlements.json
                entities/kingdoms.json

Rasters are PNG + JSON sidecar (explicit linear colour space); structs are JSON/TOML. Every
file carries `schema_version` and `kind`, validated on load (hard error on mismatch).

## Types

All data structures and common interfaces live in one crate, `geoforge-types` (the
dependency sink); logic crates depend on it and never on each other. The law is dependency
direction, not centralization: shared-contract types are centralized, module-private types
stay with their module.

## Injection

Any layer may be replaced from a hand-authored file: a raster PNG + sidecar, or a struct
JSON/TOML document. Injection replaces the whole layer; downstream layers cannot tell the
difference. Injected files become part of the reproducibility key.

## Determinism

Hard contract over `(seed, injected layers, code version)`. The seed is the sole entropy
source; wall-clock time, unseeded RNG, and hash-map iteration order are forbidden; parallelism
is limited to order-independent work or fixed-order reductions.

## Alternatives considered

- **Parent decides the child roster ("Design B").** Rejected: it splits each object's
  definition across two layers (slot by the parent, content by the child). The chosen framing
  gives each layer one owner for its interior.
- **HEALPix / equal-area mesh grid.** Rejected for v2: invalidates v1's spherical unwrapping
  and complicates hand-editing; recorded as a possible v3 direction.
- **Raw/NPY sidecar instead of PNG.** Rejected: hand-editability is a requirement.
- **Merge-injection instead of replace.** Rejected: merge semantics differ per layer type;
  replacement is clean and partial injection can be layered on later.
- **One crate per stage.** Rejected as ceremony for a POC; stage logic groups into a few
  crates with a module per layer.
- **Nesting whole child entities inside parent files.** Rejected: file sizes explode; children
  live in their own files, referenced by name + location.

## Open questions

- Cross-platform bit-exact determinism vs same-toolchain only.
- `planetary_body` interior fields per body kind (gas giant vs rocky vs asteroid).
- Belt membership: sampled members + statistical population vs full materialization.
- Roster/index storage for collection layers (parent document vs per-level index).
