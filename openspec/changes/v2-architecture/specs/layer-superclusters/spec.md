## ADDED Requirements

### Requirement: Superclusters layer produces the universe's interior
The `superclusters` layer MUST produce the superclusters that form the interior of the
universe: how many, where in the universe, and their properties.

#### Scenario: superclusters produced
- **WHEN** the `superclusters` layer runs
- **THEN** it emits a set of superclusters, each with a name, a location in the universe, and
  its properties (mass, radius, type)

### Requirement: Superclusters layer inputs
The `superclusters` layer MUST depend only on the `universe` layer and the seed.

#### Scenario: declared inputs
- **WHEN** the layer runs
- **THEN** it reads only `universe`

### Requirement: Superclusters layer storage
The layer MUST write a roster of its superclusters (name + location) at the universe node,
and a document per supercluster.

#### Scenario: roster and documents
- **WHEN** the layer completes
- **THEN** `cosmology/superclusters/index.toml` lists every supercluster by name and location,
  and `cosmology/superclusters/<c>/cluster.json` holds each cluster's properties

### Requirement: Superclusters layer injection
The `superclusters` layer MUST be replaceable from editable documents.

#### Scenario: injected superclusters
- **WHEN** a hand-authored supercluster set is inserted
- **THEN** downstream layers read it in place of the generated set

### Requirement: Superclusters layer determinism
The `superclusters` layer MUST be deterministic from its declared inputs.

#### Scenario: reproducible
- **WHEN** the layer runs twice with identical inputs
- **THEN** it produces identical superclusters

### Requirement: Superclusters layer performance
The `superclusters` layer MUST scale linearly with the number of superclusters.

#### Scenario: linear cost
- **WHEN** the number of superclusters doubles
- **THEN** the layer's work roughly doubles
