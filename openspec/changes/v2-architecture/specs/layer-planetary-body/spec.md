## ADDED Requirements

### Requirement: Planetary-body layer determines presence and interior
The `planetary_body` layer MUST, from the body's basic parameters and the seed, determine the
body's presence flags — whether it has an atmosphere, a hydrosphere, and a lithosphere — and
produce the interior for whichever are present.

#### Scenario: presence determined
- **WHEN** the `planetary_body` layer runs on a body
- **THEN** it emits `has_atmosphere`, `has_hydrosphere`, and `has_lithosphere` flags, the
  material composition of each present sphere, and the interior for whichever are present

### Requirement: Sphere composition
The layer MUST record each present sphere's material composition — the lithosphere's shell
material (rock, ice, or mixed) and the hydrosphere's liquid (water, methane, and so on) —
rather than a separate cryosphere flag.

#### Scenario: icy shell
- **WHEN** the body is an ice world whose shell is ice
- **THEN** `has_lithosphere` is true with an ice composition, and no separate cryosphere flag is
  emitted
#### Scenario: surface ice on a rocky body
- **WHEN** a rocky body carries surface ice
- **THEN** the ice is recorded as a frozen fraction of its hydrosphere

### Requirement: Sphere chemistry
The layer MUST express each present sphere's makeup as a set of chemical species drawn from a
shared chemical enum, each with a fractional abundance, and the fractions MUST sum to 1 for the
atmosphere, the hydrosphere, and the lithosphere alike.

#### Scenario: atmosphere makeup
- **WHEN** a body has an atmosphere
- **THEN** its makeup is a set of (chemical species, fraction) pairs summing to 1
#### Scenario: hydrosphere makeup
- **WHEN** a body has a hydrosphere
- **THEN** its makeup is a set of (chemical species, fraction) pairs summing to 1
#### Scenario: lithosphere makeup
- **WHEN** a body has a lithosphere
- **THEN** its makeup is a set of (chemical species, fraction) pairs summing to 1

### Requirement: Presence flags gate downstream layers
The presence flags MUST be the preconditions the body-construction layers use to decide whether
they apply.

#### Scenario: gas giant
- **WHEN** the body is a gas giant
- **THEN** `has_lithosphere` is false, so the solid-surface layers skip and only atmospheric
  construction applies

### Requirement: Body kinds
The `planetary_body` layer MUST handle every non-stellar body kind — terrestrial planet, gas
giant, ice giant, dwarf planet, moon, asteroid, comet — with interior fields appropriate to
the kind.

#### Scenario: gas giant interior
- **WHEN** the body is a gas giant
- **THEN** the interior is atmosphere-only (composition, pressure–temperature profile, cloud
  decks, zonal winds), with no lithosphere
#### Scenario: airless body
- **WHEN** the body is an airless asteroid
- **THEN** atmosphere and hydrosphere are absent, and the lithosphere carries the interior

### Requirement: Planetary-body layer inputs
The `planetary_body` layer MUST depend only on the `solar_system` layer (the body's orbital
node) and the seed.

#### Scenario: declared inputs
- **WHEN** the layer runs
- **THEN** it reads only `solar_system` and the seed

### Requirement: Planetary-body layer storage
The layer MUST write a body interior document at `<b>/body.json`.

#### Scenario: interior document
- **WHEN** the layer completes for a body
- **THEN** `<b>/body.json` holds the body's atmosphere, hydrosphere, and geosphere parameters

### Requirement: Planetary-body layer injection
The `planetary_body` layer MUST be replaceable from an editable document.

#### Scenario: injected interior
- **WHEN** a hand-authored body document is inserted
- **THEN** downstream layers read it in place of the generated interior

### Requirement: Planetary-body layer determinism
The `planetary_body` layer MUST be deterministic from its declared inputs.

#### Scenario: reproducible
- **WHEN** the layer runs twice with identical inputs
- **THEN** it produces identical interiors

### Requirement: Planetary-body layer performance
The `planetary_body` layer MUST scale linearly with the number of bodies.

#### Scenario: linear cost
- **WHEN** the number of bodies doubles
- **THEN** the layer's work roughly doubles
