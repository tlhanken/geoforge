## ADDED Requirements

### Requirement: The sky is read from the galaxy's stars
A system's night sky MUST show the galaxy's generated (and edited) stars around the system's
place, brighter than a limiting magnitude, through the dust, without the system's own stars.

#### Scenario: a sky like the Sun's
- **WHEN** the sky is drawn at the Sun's place in the Milky Way model
- **THEN** it has about 9,000 stars to magnitude 6.5 and about 20 brighter than 1.5, giants
  the commonest of them

#### Scenario: neighbours see each other
- **WHEN** two systems are made at stars a few parsecs apart
- **THEN** each appears in the other's sky at the distance between them

### Requirement: The glow is the galaxy's unresolved light
The glow MUST be the components' unresolved starlight along each line of sight through the
dust, and MUST be a sum of sources so other galaxies can be added later.

#### Scenario: the band and the bulge
- **WHEN** the glow is read from the Sun's place
- **THEN** it is brightest along the plane, brightest of all toward the bulge, and cut by
  dust lanes

### Requirement: A system has a place in its galaxy
Every system MUST have a position in its galaxy: a standalone system's drawn from its
world's seed, a system made at a star the star's.

#### Scenario: a standalone system
- **WHEN** a system is made on its own
- **THEN** it has a position in its world's galaxy, and its sky is drawn from there
