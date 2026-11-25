## ADDED Requirements

### Requirement: Fictional materials listed per universe
A universe directory MUST read its fictional materials from `exotics.toml` (none if
absent), checking names, chances, kinds and colours, and MUST place each from its own
stream keyed by its name, so that adding, removing or editing one changes no other material's
places and nothing else that is made.

#### Scenario: a new universe
- **WHEN** a universe directory is created
- **THEN** it holds a commented `exotics.toml` that lists nothing, and uncommenting its
  example gives one material

#### Scenario: another material added
- **WHEN** a second material is added to a list
- **THEN** the first is found in the same places as before
