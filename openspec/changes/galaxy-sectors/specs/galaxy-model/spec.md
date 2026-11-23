## ADDED Requirements

### Requirement: A galaxy is made of components
A galaxy MUST be described as components (nucleus, bulge and bar, thin disk, spiral arms,
thick disk, halo, dust), each giving its stellar density at a point, the age and metals of
the stars born in it, and its share of stars by brightness.

#### Scenario: the Sun's neighbourhood
- **WHEN** the Milky Way model is read at the Sun's place
- **THEN** it holds about 0.09 stars per cubic parsec, the thick disk about 2% of the thin
  disk's density and the halo about 1% of the disk's

#### Scenario: old stars sit higher
- **WHEN** the thin disk's stars are sorted by age
- **THEN** their scale height grows with age, from tens of parsecs for the youngest to about
  300 pc for the oldest

### Requirement: Morphologies
The model MUST offer spiral (Sa to Sd, barred or not), lenticular, elliptical, dwarf
spheroidal, dwarf irregular and starburst galaxies as sets of components, and a galaxy MUST
be loadable from an editable document.

#### Scenario: an elliptical has no young stars
- **WHEN** an elliptical galaxy is drawn from
- **THEN** none of its stars is younger than a few gigayears and it has no spiral arms

#### Scenario: a world's own galaxy
- **WHEN** a world holds a galaxy document
- **THEN** its stars and skies are drawn from that galaxy
