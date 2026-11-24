## ADDED Requirements

### Requirement: Halos by mass
A universe MUST hold dark-matter halos as many as the halo mass function gives, clustered
more strongly the heavier they are.

#### Scenario: abundance
- **WHEN** the halos above 10¹² M☉ in a large volume are counted
- **THEN** there are as many as the mass function gives to within 15%

#### Scenario: clustering
- **WHEN** Milky Way halos are counted in 10 Mpc cells
- **THEN** their counts scatter more than a Poisson process's

### Requirement: Galaxies in halos
Each halo MUST hold a central galaxy whose stellar mass follows its halo's mass, and
satellites in subhalos that heavy ones sink from.

#### Scenario: a Milky Way halo
- **WHEN** halos of 1.3×10¹² M☉ are filled
- **THEN** the central holds 3–6×10¹⁰ M☉ of stars, about a dozen satellites hold over
  10⁵ M☉, and roughly one halo in eight has an LMC

### Requirement: Galaxies generated where looked at
Every galaxy MUST have a stable identity that regenerates it, and a query MUST generate only
the cubes it reaches.

#### Scenario: an identity read back
- **WHEN** a galaxy's identity is written out and looked up again
- **THEN** the same galaxy comes back

### Requirement: A world's galaxy in its universe
A world's galaxy MUST be placed in a universe as a generated galaxy like it; the Solar
System's universe MUST be the real one.

#### Scenario: a Milky Way model
- **WHEN** a world in the Milky Way model is placed
- **THEN** its home galaxy is a disk within 0.25 dex of the model's stellar mass, with
  satellites of its own
