## ADDED Requirements

### Requirement: Addresses
Every galaxy, star's system and body of a universe MUST have an address that reads back to
itself and from which the item is generated, with no other input but the edits kept over
it.

#### Scenario: read back
- **WHEN** an address is written and parsed again
- **THEN** it is the same address, and nonsense (a Milky Way in a seeded universe, a level
  below a body) is refused

#### Scenario: the same item from the same address
- **WHEN** a body is resolved twice from its address
- **THEN** its galaxy, system and body are the same both times

### Requirement: One universe per seed
All the bodies made from one seed MUST be in one universe: a star's system there is drawn
from the same galaxy's stars, under the same star edits.

#### Scenario: visiting a star
- **WHEN** a world goes to a star of its sky
- **THEN** the new body is in the same universe directory, its galaxy the same, its system
  the star's

### Requirement: Default body
A universe's default body MUST be a wet world of a Sun-like star in its galaxy most like
the Milky Way, found by a fixed rule.

#### Scenario: universe 7
- **WHEN** universe 7's default body is made
- **THEN** its star is a single main-sequence star of 5,200–6,200 K, the body has surface
  water, and its galaxy is a barred Sbc like the Milky Way's model

### Requirement: Galaxy models from galaxies
A generated galaxy's model MUST be its type's preset scaled to its stellar mass and size.

#### Scenario: a Milky Way-like galaxy
- **WHEN** a barred Sbc's model is built
- **THEN** its stellar mass is the galaxy's and its lengths follow its half-light radius
