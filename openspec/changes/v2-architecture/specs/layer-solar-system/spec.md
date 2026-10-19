## ADDED Requirements

### Requirement: Solar-system layer produces a system's interior
The `solar_system` layer MUST produce the interior of a single solar system: its stars and its
orbital hierarchy of bodies, moons, belts, rings, asteroids, and comets.

#### Scenario: interior produced
- **WHEN** the `solar_system` layer runs on a solar system
- **THEN** it emits that system's stars and a set of orbital nodes, each with a kind, a parent,
  and orbital elements

### Requirement: Solar-system layer emits body parameters
The `solar_system` layer MUST emit, for each body, the basic parameters body generation
consumes: kind, mass, radius, rotation, and orbital elements (semi-major axis, eccentricity,
inclination, longitude of ascending node, argument of periapsis, mean anomaly).

#### Scenario: body parameters emitted
- **WHEN** the `solar_system` layer runs
- **THEN** each body's orbital node carries the basic parameters needed to generate its interior

### Requirement: Orbital hierarchy
The system interior MUST be expressible as a tree: bodies orbit stars, moons and rings orbit
bodies, belts orbit a star or a body, and asteroids and comets are members of belts.

#### Scenario: moon around a planet
- **WHEN** a moon is generated
- **THEN** its orbital node names a body as its parent
#### Scenario: belt around a star or a body
- **WHEN** a belt is generated
- **THEN** its parent is either a star or a body
#### Scenario: entities within a belt
- **WHEN** a belt is generated
- **THEN** it carries a population estimate and a bounded set of sampled members

### Requirement: Solar-system layer inputs
The `solar_system` layer MUST depend only on the `solar_systems` layer (the system's roster
entry) and the seed.

#### Scenario: declared inputs
- **WHEN** the layer runs
- **THEN** it reads only `solar_systems` and the seed

### Requirement: Solar-system layer storage
The layer MUST write a system interior document at `<s>/system.json`.

#### Scenario: interior document
- **WHEN** the layer completes for a system
- **THEN** `<s>/system.json` holds the stars and orbital nodes

### Requirement: Solar-system layer injection
The `solar_system` layer MUST be replaceable from an editable document.

#### Scenario: injected interior
- **WHEN** a hand-authored system document is inserted
- **THEN** downstream layers read it in place of the generated interior

### Requirement: Solar-system layer determinism
The `solar_system` layer MUST be deterministic from its declared inputs.

#### Scenario: reproducible
- **WHEN** the layer runs twice with identical inputs
- **THEN** it produces identical interiors

### Requirement: Solar-system layer performance
The `solar_system` layer MUST scale linearly with the number of orbital nodes.

#### Scenario: linear cost
- **WHEN** the number of nodes doubles
- **THEN** the layer's work roughly doubles
