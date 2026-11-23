## ADDED Requirements

### Requirement: Stars live in nested cubes by brightness
Stars MUST be generated in cubes of the galaxy's frame, each brightness bin in cubes sized
to how far its stars can be seen, and each cube's stars MUST be a pure function of the
galaxy's seed and the cube's key.

#### Scenario: the same star from anywhere
- **WHEN** two observers a few parsecs apart query the stars around them
- **THEN** every star both see has the same identity, position and kind for both

#### Scenario: zooming in
- **WHEN** a region is queried for fainter stars than before
- **THEN** the brighter stars are unchanged and only the cubes of the fainter bins near it
  are generated

### Requirement: Stars have stable identities
Every generated star MUST have an identity (its bin, cube and index) that is the seed of the
system made at it.

#### Scenario: visiting a star in the sky
- **WHEN** a system is made at a star seen in a sky
- **THEN** its primary has that star's mass, age and metals, and making it again gives the
  same system

### Requirement: Counts follow the galaxy
The number of stars in a cube MUST be a Poisson draw on the components' density integrated
over it, and their positions MUST follow the density within the cube.

#### Scenario: counts over many cubes
- **WHEN** many cubes of one region are generated
- **THEN** their mean count matches the integrated density, and stars sit nearer the plane
  in proportion to the layer's scale height

### Requirement: Young stars are born in groups
Stars younger than about 300 Myr MUST belong to birth events (clusters and associations)
placed where stars form, a small share staying bound as open clusters and the rest
spreading with age; older stars MUST be smooth field stars.

#### Scenario: young stars clump
- **WHEN** the B stars within a kiloparsec are generated
- **THEN** they cluster in groups rather than spreading evenly

### Requirement: Dust dims along the line of sight
Extinction to a star MUST be integrated along the line of sight through a smooth dust layer
and discrete clouds.

#### Scenario: a star behind a cloud
- **WHEN** a star lies behind a dust cloud
- **THEN** it is dimmer than a star at the same distance beside the cloud
