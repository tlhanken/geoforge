## ADDED Requirements

### Requirement: Every product reaches the outputs
A new product MUST be declared in `ProductId` and `ProductData`, rendered with a legend
that states its units, listed in the report in pipeline order, and comparable by
`geoforge compare` when it is a quantity, a flow or a set of classes.

#### Scenario: a new product appears everywhere
- **WHEN** a generator adds a product
- **THEN** `geoforge render` draws it, `geoforge report` shows it in its layer's section,
  and `geoforge compare` sets it beside an observed map of the same product

### Requirement: Generators version their output
A generator MUST bump its `VERSION` whenever its output changes for the same inputs.

#### Scenario: stale worlds regenerate
- **WHEN** a generator's algorithm changes
- **THEN** its version rises, and worlds generated before regenerate that layer and those
  downstream

### Requirement: Kinds stay in their place
Stored, editable results MUST be products of a generator; results read from products on
demand MUST be derived views; results that depend on the date or a place over time MUST be
almanac features.

#### Scenario: a date-dependent result
- **WHEN** a feature changes with the date
- **THEN** it is an almanac feature, not a product

### Requirement: Checked against Earth
A layer that models something observable on Earth MUST have known cases with observed
values, and MUST be run before and after a change; a case it misses MUST be recorded as a gap
whose text states what was measured, never hidden by widening the expected range.

#### Scenario: a change loses a case
- **WHEN** a change makes a known case fail
- **THEN** the case becomes a gap stating the measured value and its cause, and the change
  reports what it lost as well as what it gained

### Requirement: Draws are seeded from the world
Every random draw MUST be seeded from the world's root seed and the index of the step,
month or cell it belongs to, never from the query or the wall clock.

#### Scenario: the same date asked two ways
- **WHEN** a day is asked for alone and as part of a year
- **THEN** both give the same weather
