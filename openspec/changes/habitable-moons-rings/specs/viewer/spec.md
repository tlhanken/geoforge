## ADDED Requirements

### Requirement: Habitable systems flagged
The viewer MUST mark, in the stars around a system, those whose most Earth-like world is
temperate or habitable, and count them.

#### Scenario: the default star
- **WHEN** the stars around universe 7's default body are graded
- **THEN** its own star is flagged temperate, and asking again gives the same flags from cache

### Requirement: Rings drawn
The viewer MUST draw a planet's rings in the system and body views, tilted with its equator.

#### Scenario: a ringed giant
- **WHEN** a giant with bright rings is opened
- **THEN** its scene carries the rings' extent in its radii, their optical depth and tilt
