## ADDED Requirements

### Requirement: Nebulae and globular clusters in the sky
A world's sky MUST list the nebulae and globular clusters it sees through the dust, with
their brightness and size, and its chart MUST draw them toned by surface brightness.

#### Scenario: a sky chart
- **WHEN** `geoforge sky` draws a world's chart
- **THEN** nebulae and clusters brighter than magnitude 9 are listed and drawn
