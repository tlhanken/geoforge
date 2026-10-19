## ADDED Requirements

### Requirement: Universe layer produces the universal constants
The `universe` layer MUST produce the universe's global constants — the environment every
later layer reads — and nothing else.

#### Scenario: constants are produced
- **WHEN** the `universe` layer runs
- **THEN** it emits a `UniverseConstants` value holding the global cosmological parameters

### Requirement: Universe layer inputs
The `universe` layer MUST depend on the root seed only.

#### Scenario: seed-only input
- **WHEN** the `universe` layer runs on a fresh world
- **THEN** it requires no upstream layer

### Requirement: Universe layer storage
The universe constants MUST be stored at `cosmology/universe.toml`.

#### Scenario: file round-trips
- **WHEN** the layer completes
- **THEN** `cosmology/universe.toml` exists and reloads to the same values

### Requirement: Universe layer injection
The `universe` layer MUST be replaceable from an editable JSON or TOML document.

#### Scenario: injected constants
- **WHEN** a hand-authored universe document is inserted
- **THEN** downstream layers read it in place of the generated constants

### Requirement: Universe layer determinism
The `universe` layer MUST be deterministic from the seed.

#### Scenario: reproducible
- **WHEN** the layer runs twice with the same seed
- **THEN** it produces identical constants

### Requirement: Universe layer performance
The `universe` layer MUST run in constant time and space.

#### Scenario: constant cost
- **WHEN** the layer runs
- **THEN** its cost does not scale with world size
