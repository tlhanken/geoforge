## ADDED Requirements

### Requirement: Stars can be edited
A world MUST be able to remove, change and add stars and add clusters, the edits kept beside
the generated stars and applied over them in every query.

#### Scenario: a moved star
- **WHEN** a generated star is moved
- **THEN** every sky shows it at its new place and none at its old one

#### Scenario: an added cluster
- **WHEN** a cluster is added with a centre, radius, number of members and age
- **THEN** its members appear in skies that reach them, the same every time

#### Scenario: edits outlive the generator
- **WHEN** the generator changes after a star was changed
- **THEN** the changed star keeps its stored record, and an edit naming a star no longer
  generated is reported

### Requirement: Constellations belong to an observer
A world MUST be able to name constellations as lines between stars of its focus's sky, and
the sky chart MUST draw them.

#### Scenario: a custom constellation
- **WHEN** a constellation is added from stars of the sky
- **THEN** the chart draws its lines and names it

### Requirement: Stars can be exported and imported
The stars of a sky MUST be exportable with their identities, positions, distances and
magnitudes, and an edited export MUST import as edits.

#### Scenario: a round trip through a spreadsheet
- **WHEN** a sky's stars are exported, one is moved in the file and it is imported
- **THEN** the world holds one change edit for that star
