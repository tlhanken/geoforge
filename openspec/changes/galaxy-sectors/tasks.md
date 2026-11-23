# Tasks — galaxy-sectors

## 1. Galaxy components
- [ ] 1.1 Component shapes and densities: exponential disk with age-dependent height, Sérsic spheroid, bar, flattened power-law halo, logarithmic arms, nucleus
- [ ] 1.2 The Milky Way as components, its numbers checked against sources (local densities, thick disk and halo shares, scale lengths and heights, bar, arms)
- [ ] 1.3 Morphologies (Sa–Sd barred or not, S0, E, dE/dSph, dIrr, starburst) replacing the presets; galaxy documents load
- [ ] 1.4 Luminosity tables per component and brightness bin

## 2. Sectors
- [ ] 2.1 Cube levels per bin from its reach; cube keys
- [ ] 2.2 Stars per cube: Poisson counts on the integrated density, positions by the density, kinds from the tables
- [ ] 2.3 Star identities and their seeds
- [ ] 2.4 Queries: within a radius to an absolute magnitude; brighter than an apparent magnitude from a point

## 3. Young groups and dust
- [ ] 3.1 Birth events: rate from star formation, masses, ages, bound share, spreading, members per bin on demand; field tables over 300 Myr
- [ ] 3.2 Dust layer and clouds; extinction along a ray

## 4. The night sky
- [ ] 4.1 The sky from sectors, the system's own stars removed, in the system's frame
- [ ] 4.2 The glow from all components through the dust, cached; a sum of sources
- [ ] 4.3 Checks at the Sun's place: naked-eye count, first magnitude, giants, the Gould Belt, the 10 pc census, the band and bulge

## 5. A system's place
- [ ] 5.1 `GalacticSite` with position and star identity
- [ ] 5.2 Standalone systems draw a position; the habitable recipe keeps a Sun-like primary
- [ ] 5.3 A system made at a star from its phenotype and seed
- [ ] 5.4 Existing worlds placed at their stored radius

## 6. Edits
- [ ] 6.1 Edit types and `galaxy/stars.json` in the store
- [ ] 6.2 Edits applied in every query
- [ ] 6.3 `geoforge star add|move|change|remove`; `geoforge stars export|import`
- [ ] 6.4 Constellations: stored per observer, drawn on the chart, `geoforge constellation add|remove`
