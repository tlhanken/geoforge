# Tasks — v2-architecture

## 1. Infrastructure specs
- [x] 1.1 `cli` — command surface, seed rule, overwrite protection, ranges, exit codes
- [x] 1.2 `world-container` — directory layout, manifest, file formats, validation
- [x] 1.3 `determinism` — reproducibility contract, entropy rules, parallel rules
- [x] 1.4 `layer-applicability` — declared preconditions, skip semantics, skip markers

## 2. Cosmology layers
- [x] 2.1 `layer-universe` — universal constants
- [x] 2.2 `layer-superclusters` — the superclusters in the universe
- [x] 2.3 `layer-galaxies` — the galaxies in a supercluster
- [x] 2.4 `layer-solar-systems` — the solar systems in a galaxy
- [x] 2.5 `layer-solar-system` — stars, orbits, moons, belts (the system interior)
- [x] 2.6 `layer-planetary-body` — atmosphere, hydrosphere, geosphere of a body

## 3. Surface layers
- [x] 3.1 `layer-tectonics` — plate id, crust type/age, boundaries, per-plate velocity
- [ ] 3.2 `layer-geology` — geologic provinces
- [ ] 3.3 `layer-elevation` — elevation, ocean flag, slope
- [ ] 3.4 `layer-climate` — temperature, precipitation, wind
- [ ] 3.5 `layer-hydrology` — river flow, lakes, watersheds
- [ ] 3.6 `layer-biomes` — biome classification
- [ ] 3.7 `layer-resources` — resource kind and abundance
- [ ] 3.8 `layer-hazards` — hazard kind and intensity
- [ ] 3.9 `layer-settlements` — settlement entities
- [ ] 3.10 `layer-kingdoms` — kingdom regions and entities

## 4. Reserved
- [ ] 4.1 `layer-culture` — placeholder requirement only
- [ ] 4.2 `layer-religion` — placeholder requirement only

## 5. Reconciliation
- [x] 5.1 Fold `docs/proposals/0001` and `docs/specs/cli.md` into this change and remove them;
      `docs/proposals/0002` (bloomery) retained on disk, deferred to a later change
- [x] 5.2 Validate: `openspec validate v2-architecture --type change --no-interactive`
