## ADDED Requirements

### Requirement: Solar-systems layer produces a galaxy's interior
The `solar_systems` layer MUST produce the solar systems that form the interior of a galaxy:
how many, where in the galaxy, and their identifying properties.

#### Scenario: solar systems produced
- **WHEN** the `solar_systems` layer runs on a galaxy
- **THEN** it emits that galaxy's solar systems, each with a name, a location in the galaxy,
  and a phenotype (multiplicity, primary spectral class)

### Requirement: Solar-systems layer inputs
The `solar_systems` layer MUST depend only on the `universe`, `superclusters`, and `galaxies`
layers and the seed.

#### Scenario: declared inputs
- **WHEN** the layer runs
- **THEN** it reads only `universe`, `superclusters`, and `galaxies`

### Requirement: Solar-systems layer storage
The layer MUST write a roster of a galaxy's solar systems (name, location, phenotype) at the
galaxy node. It MUST NOT write a system's interior.

#### Scenario: roster written
- **WHEN** the layer completes for a galaxy
- **THEN** `<g>/solar_systems/index.toml` lists every solar system by name, location, and
  phenotype, and no system interior is produced

### Requirement: Solar-systems layer injection
The `solar_systems` layer MUST be replaceable from editable documents.

#### Scenario: injected systems
- **WHEN** a hand-authored solar-system roster is inserted
- **THEN** downstream layers read it in place of the generated roster

### Requirement: Solar-systems layer determinism
The `solar_systems` layer MUST be deterministic from its declared inputs.

#### Scenario: reproducible
- **WHEN** the layer runs twice with identical inputs
- **THEN** it produces identical rosters

### Requirement: Solar-systems layer performance
The `solar_systems` layer MUST scale linearly with the number of solar systems.

#### Scenario: linear cost
- **WHEN** the number of solar systems doubles
- **THEN** the layer's work roughly doubles
