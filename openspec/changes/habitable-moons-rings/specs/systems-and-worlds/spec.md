## ADDED Requirements

### Requirement: Habitability before generating
The system MUST estimate a solid world's climate from the latitude-band energy balance before
any surface layer is generated, and grade it temperate, habitable or not habitable by its
surface water, mean temperature and thawed share.

#### Scenario: known worlds
- **WHEN** Earth, Mars, Venus, the Moon and Europa are graded
- **THEN** Earth is temperate with a mean of 8–20 °C and the others are not habitable

#### Scenario: a cold Earth
- **WHEN** an Earth-like world gets 0.55 of Earth's sunlight, its carbon cycle settled
- **THEN** it is habitable but not temperate

### Requirement: A temperate default world
The default body of a seeded universe MUST be a temperate world of a single Sun-like star
when one lies within 240 pc of the default point; otherwise the earlier rule applies.

#### Scenario: seeds
- **WHEN** the default body of seeds 1, 2, 3, 42, 257, 525, 820 and 12345 is resolved
- **THEN** each is temperate

### Requirement: Moon systems
Giants MUST get regular moon systems holding the share of their mass their disks set, small
inner moons and irregular moons; ice giants sometimes a captured backward moon instead;
rocky planets beside an asteroid belt sometimes captured specks.

#### Scenario: ice giants
- **WHEN** many ice giants are drawn
- **THEN** about a third hold a backward captured moon and one small moon, the rest 3–5
  regular moons and 1–2 small ones

### Requirement: Rings
Every giant MUST have rings: faint dust, or on about a quarter of those under 150 K bright
ice inside the Roche limit with an optical depth of 0.2–1.6.

#### Scenario: a warm giant
- **WHEN** a giant's equilibrium temperature is 300 K
- **THEN** its rings are dust with an optical depth under 10⁻⁴

### Requirement: Zodiacal light
A world's sky MUST carry the zodiacal light along its star's path, scaled by its system's
belts against the Solar System's.

#### Scenario: Earth
- **WHEN** Earth's sky is drawn at the March equinox
- **THEN** the light is brightest toward the Sun along its path, fades off it, is gone within
  10° of the Sun and rises again opposite it
