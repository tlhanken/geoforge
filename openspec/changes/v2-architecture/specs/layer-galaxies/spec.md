## ADDED Requirements

### Requirement: Galaxies layer produces a supercluster's interior
The `galaxies` layer MUST produce the galaxies that form the interior of a supercluster: how
many, where in the supercluster, and their properties.

#### Scenario: galaxies produced
- **WHEN** the `galaxies` layer runs on a supercluster
- **THEN** it emits that supercluster's galaxies, each with a name, a location in the
  supercluster, and its properties (morphology, mass, size)

### Requirement: Galaxies layer inputs
The `galaxies` layer MUST depend only on the `universe` and `superclusters` layers and the seed.

#### Scenario: declared inputs
- **WHEN** the layer runs
- **THEN** it reads only `universe` and `superclusters`

### Requirement: Galaxies layer storage
The layer MUST write a roster of a supercluster's galaxies (name + location) at the
supercluster node, and a document per galaxy.

#### Scenario: roster and documents
- **WHEN** the layer completes for a supercluster
- **THEN** `<c>/galaxies/index.toml` lists every galaxy by name and location, and
  `<c>/galaxies/<g>/galaxy.json` holds each galaxy's properties

### Requirement: Galaxies layer injection
The `galaxies` layer MUST be replaceable from editable documents.

#### Scenario: injected galaxies
- **WHEN** a hand-authored galaxy set is inserted
- **THEN** downstream layers read it in place of the generated set

### Requirement: Galaxies layer determinism
The `galaxies` layer MUST be deterministic from its declared inputs.

#### Scenario: reproducible
- **WHEN** the layer runs twice with identical inputs
- **THEN** it produces identical galaxies

### Requirement: Galaxies layer performance
The `galaxies` layer MUST scale linearly with the number of galaxies.

#### Scenario: linear cost
- **WHEN** the number of galaxies doubles
- **THEN** the layer's work roughly doubles
