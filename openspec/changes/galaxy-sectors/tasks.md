# Tasks — galaxy-sectors

## 1. Galaxy components
- [x] 1.1 Component shapes and densities: exponential disk with age-dependent height, exponential ellipsoid (bulge, bar), Hernquist spheroid, flattened broken power-law halo, Plummer nucleus, logarithmic arms (`galaxy/{shape, arms}.rs`)
- [x] 1.2 The Milky Way as components, its numbers from Bland-Hawthorn & Gerhard (2016) and Vallée (2015); the thin disk's mass set by the 10 pc census (`galaxy/presets.rs`)
- [x] 1.3 Morphologies (Sa–Sd barred or not, S0, E, dSph, dIrr, starburst) replacing the presets; galaxy documents load (`new --galaxy-file`, `galaxy/galaxy.json`)
- [x] 1.4 Luminosity tables per population and brightness bin (`sectors/tables.rs`); the galaxy-integrated IMF (α = 2.6 above 1 M☉)

## 2. Sectors
- [x] 2.1 Cube levels per bin from its reach; cube keys (`sectors/bins.rs`)
- [x] 2.2 Stars per cube: Poisson counts on the integrated density, positions by the density, kinds from the tables (`sectors/field.rs`); a disk's plane split once for all its age slices
- [x] 2.3 Star identities and their seeds (`geoforge-core::galactic::StarId`)
- [x] 2.4 Queries: within a radius to an absolute magnitude (`near`); brighter than an apparent magnitude from a point (`visible`); a star by identity (`star`); cubes in parallel

## 3. Young groups and dust
- [x] 3.1 Birth events of 100 M☉ and more (smaller groups as young field stars): rate from star formation, masses, ages, bound share, spreading, members per bin on demand (`sectors/groups.rs`)
- [x] 3.2 Dust layer and clouds; extinction along a ray; an angular index of clouds per observer (`sectors/dust.rs`)

## 4. The night sky
- [x] 4.1 The sky from sectors, the system's own stars removed, in the system's frame (`sky.rs`); crowded places listed to a brighter limit
- [x] 4.2 The glow from all populations' unresolved light through the dust
- [x] 4.3 Checks at the Sun's place: about 9,000 naked-eye stars, 20–25 of the first magnitude, giants the commonest, B stars in groups, 0.08 systems per pc³, the band bright toward Sagittarius and Cygnus on the real Earth

## 5. A system's place
- [x] 5.1 `GalacticSite` with position, galaxy seed, star identity and orientation
- [x] 5.2 Standalone systems draw a position from every component by mass; the habitable recipe a Sun-like place; the real Earth the Sun's, its ecliptic turned as it is
- [x] 5.3 A system made at a star from its phenotype and seed (`galactic::system_at`, `geoforge visit`)
- [x] 5.4 Existing worlds placed at their stored radius

## 6. Edits
- [x] 6.1 Edit types (`StarEdits`, `StarRecord`, `ClusterEdit`, `Constellation`) and `galaxy/stars.json`, `constellations.json` in the store
- [x] 6.2 Edits applied in every query (`starfield.rs`)
- [x] 6.3 `geoforge star add|move|change|remove|cluster`; `geoforge stars near|export|import`
- [x] 6.4 Constellations: stored per observer, drawn on the chart, `geoforge constellation add|remove|list`

## 7. Next
- [ ] 7.1 Brown dwarfs: the IMF below 0.08 M☉ and their cooling (item 14)
- [ ] 7.2 Other galaxies in the night sky, from the universe layers (TODO 12b)
- [ ] 7.3 A Gould-Belt-like local structure: the Sun's sky sits in a local excess of young stars
