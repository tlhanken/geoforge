## ADDED Requirements

### Requirement: Layers declare applicability
Every body-construction layer MUST declare the preconditions under which it applies — typically
the presence flags of the body it operates on.

#### Scenario: declared preconditions
- **WHEN** a layer is defined
- **THEN** it declares the body properties (for example `has_lithosphere`, `has_hydrosphere`)
  that must hold for it to apply

### Requirement: Non-applicable layers skip
A layer whose preconditions are not met MUST skip generation. It MUST NOT run, and it MUST NOT
error.

#### Scenario: tectonics on a gas giant
- **WHEN** a layer requiring `has_lithosphere` runs on a body without one
- **THEN** it skips and produces no layer data

### Requirement: Skips are recorded
A skipped layer MUST record a not-applicable marker, distinguishable from a missing input.

#### Scenario: absence is intentional
- **WHEN** a downstream layer reads an upstream layer that was skipped
- **THEN** it reads a not-applicable marker, not a missing-input error

### Requirement: Skip determinism
Skipping MUST be a deterministic function of the declared preconditions and the inputs.

#### Scenario: reproducible skips
- **WHEN** the pipeline runs twice with identical inputs
- **THEN** exactly the same layers skip
