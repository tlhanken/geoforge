## ADDED Requirements

### Requirement: Other galaxies in the sky
The night sky MUST list the other galaxies brighter than magnitude 8 and draw them into its
glow by their surface brightness, dimmed by the home galaxy's dust.

#### Scenario: Earth
- **WHEN** Earth's sky is drawn
- **THEN** the Magellanic Clouds and Andromeda are at their real positions, the Large Cloud
  more than twice as bright as the sky 40° from it, Andromeda brighter than the sky beside
  it

#### Scenario: a generated Milky Way
- **WHEN** skies of worlds in the Milky Way model are drawn
- **THEN** most hold a naked-eye galaxy or a few, its satellites or a neighbour

### Requirement: Automatic constellations
A world without constellations of its own MUST be shown figures drawn from its brightest
stars, which it can save as its own.

#### Scenario: a Sun-like sky
- **WHEN** the automatic figures of a Sun-like sky are drawn
- **THEN** there are 40–120 of 2–10 stars each, no line is longer than 15° and every star
  brighter than magnitude 1.5 is in one
