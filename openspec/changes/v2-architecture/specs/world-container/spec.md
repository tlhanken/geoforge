## ADDED Requirements

### Requirement: World is a directory
A world MUST be a directory containing a manifest and the layer tree.

#### Scenario: world layout
- **WHEN** a world is created
- **THEN** it contains `world.toml` and a `cosmology/` tree with per-body layer directories

### Requirement: Path-namespaced cosmology
Cosmology entities MUST be stored in per-parent folders, and an entity name MUST be unique only
within its parent, not globally.

#### Scenario: duplicate names in different parents
- **WHEN** two solar systems share a name but live in different galaxies
- **THEN** they are stored at distinct paths and do not collide

### Requirement: Parent references children
Each parent entity MUST carry a roster of its children by name and location — in its own
document or an index at its node.

#### Scenario: galaxy references its systems
- **WHEN** a galaxy's children are written
- **THEN** a roster at the galaxy node lists its solar systems by name and location

### Requirement: Manifest contents
The manifest MUST record `schema_version`, seed, code version, grid, calibration, the layer
index, and provenance hashes.

#### Scenario: manifest round-trip
- **WHEN** a world is loaded
- **THEN** the manifest values are available to the pipeline

### Requirement: Layer file formats
Raster layers MUST be stored as PNG with a JSON sidecar; struct layers MUST be stored as JSON
or TOML.

#### Scenario: raster storage
- **WHEN** a raster layer is written
- **THEN** a PNG and its JSON sidecar are produced

### Requirement: Version and kind validation
Loading a layer MUST validate its `schema_version` and `kind` and MUST hard-error on mismatch.

#### Scenario: version mismatch
- **WHEN** a layer file with an incompatible `schema_version` is loaded
- **THEN** the load fails with a hard error
