## ADDED Requirements

### Requirement: Continents part into Earth-sized basins
Stable continental interiors MUST carry basins and swells, so that an Earth-like world's largest
drainage basin holds no more than twice the Amazon's share of the land.

#### Scenario: the Earth archetype
- **WHEN** the Earth archetype (seed 42) is generated at 1024×512
- **THEN** its largest basin drains under 9.4% of its land and 15–40% of the land stands above
  1000 m

### Requirement: Peaks held by rock strength
Final heights MUST stay below a limit of 9 km at Earth's gravity, scaled as 1/g, bending toward
it from three-quarters of it.

#### Scenario: a rising range
- **WHEN** a range would stand 9.5 km high at Earth's gravity
- **THEN** it stands between 8.5 and 9 km, and lower peaks keep their order

### Requirement: Deserts drain to no sea
A dry basin MUST drain to no sea when it is at least the smallest basin's area and its land
yields under the dead-river runoff on average.

#### Scenario: the Earth archetype's closed land
- **WHEN** the Earth archetype (seed 42) is generated at 1024×512
- **THEN** 8–30% of its land drains to no sea, and its highest point stands 7–9 km high
