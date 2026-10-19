## ADDED Requirements

### Requirement: Tectonics layer produces the plate model
The `tectonics` layer MUST produce, for the focus body, a per-pixel plate assignment and crust
rasters, a boundary-type raster, and a per-plate velocity vector.

#### Scenario: plate model produced
- **WHEN** the `tectonics` layer runs
- **THEN** it emits `plate_id`, `crust_type`, `crust_age_ma`, and `boundary_type` rasters, plus
  a `TectonicMeta` holding each plate's centre and velocity vector

### Requirement: Velocity is primary, convergence is derived
Plate motion MUST be stored as a per-plate velocity vector. Convergence rate MUST be derived
from relative plate velocity at a boundary and MUST NOT be stored as a primary layer.

#### Scenario: derived convergence
- **WHEN** a boundary's convergence is needed
- **THEN** it is computed from the two plates' velocity vectors and the boundary normal, not
  read from a stored field

### Requirement: Tectonics layer inputs
The `tectonics` layer MUST depend only on the `planetary_body` layer (the focus body's grid and
parameters) and the seed.

#### Scenario: declared inputs
- **WHEN** the layer runs
- **THEN** it reads only `planetary_body` and the seed

### Requirement: Tectonics visualisation exports
The layer MUST support these exports: a distinct-colour plate map, a velocity colour wheel
(hue = direction, saturation = speed), crust-type and crust-age maps, and a boundary-type map.

#### Scenario: plate map export
- **WHEN** the plate map is exported
- **THEN** every plate is rendered in a distinct colour
#### Scenario: velocity wheel export
- **WHEN** the velocity export is produced
- **THEN** hue encodes direction and saturation encodes speed

### Requirement: Tectonics layer injection
The `tectonics` layer MUST be replaceable from editable rasters and documents.

#### Scenario: injected plates
- **WHEN** a hand-authored plate model is inserted
- **THEN** downstream layers read it in place of the generated model

### Requirement: Tectonics layer determinism
The `tectonics` layer MUST be deterministic from its declared inputs.

#### Scenario: reproducible
- **WHEN** the layer runs twice with identical inputs
- **THEN** it produces identical plate models

### Requirement: Tectonics layer performance
The `tectonics` layer MUST scale linearly with the grid's pixel count.

#### Scenario: linear cost
- **WHEN** the grid's pixel count doubles
- **THEN** the layer's work roughly doubles
