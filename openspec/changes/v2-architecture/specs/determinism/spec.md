## ADDED Requirements

### Requirement: Reproducibility
The same seed, the same injected layers, and the same code version MUST produce identical
outputs.

#### Scenario: repeat run is identical
- **WHEN** the pipeline is run twice with identical inputs
- **THEN** the outputs are identical

### Requirement: Sole entropy source
The seed MUST be the sole source of entropy; wall-clock time, unseeded randomness, and
hash-map iteration order MUST NOT affect output.

#### Scenario: no ambient entropy
- **WHEN** generation runs under differing system time or environment variables
- **THEN** the outputs are unchanged

### Requirement: Frozen sub-seed derivation
Sub-seeds MUST be derived from the root seed through a frozen, documented function.

#### Scenario: stable derivation
- **WHEN** a layer derives a sub-seed from the root seed
- **THEN** it uses the frozen derivation, and changing that derivation is a breaking change

### Requirement: Ordered parallelism
Parallelism MUST be limited to order-independent work, or to reductions whose combination
order is fixed and specified.

#### Scenario: deterministic parallel reduction
- **WHEN** a layer performs parallel iteration
- **THEN** the result is independent of thread scheduling
