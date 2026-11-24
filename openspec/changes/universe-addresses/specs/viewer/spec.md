## ADDED Requirements

### Requirement: Zoom through a universe
The viewer MUST show a universe directory at four levels — the universe's galaxies, a
galaxy's stars, a system, a body — generating each from its address when it is opened.

#### Scenario: from the universe to a body
- **WHEN** a galaxy, then a star, then a body is picked
- **THEN** each level is shown, and the body's layers can be drawn on a sphere or a map

### Requirement: Search
The viewer MUST find items by name among those made or generated near what is shown, and
go to any address typed.

#### Scenario: an address
- **WHEN** a valid address is typed
- **THEN** the viewer shows that item, generating it if needed

### Requirement: True colour
Every body MUST have a true-colour appearance for the system view.

#### Scenario: the system
- **WHEN** a system is shown
- **THEN** its star, planets and moons are drawn in colours from their surfaces, clouds,
  bands or temperature
