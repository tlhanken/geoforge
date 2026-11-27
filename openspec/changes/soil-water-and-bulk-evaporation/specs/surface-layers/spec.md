## ADDED Requirements

### Requirement: Land evaporates its soil water
Land MUST evaporate into the air the Budyko share of its year's rain and evaporative demand,
spread over the seasons by their demand, so that a continental interior keeps evaporating
through a summer its own rain could not feed.

#### Scenario: a desert and a wet forest
- **WHEN** land gets 50 mm a year against a demand of 2000 mm, or 3000 mm against 600
- **THEN** it evaporates 48–50 mm, or 570–600 mm, and never more than its rain or its demand

#### Scenario: a continental summer
- **WHEN** a place's evaporative demand falls all in its summer
- **THEN** all of its year's evaporation is spread over its summer season and none over its
  winter, and the two seasons together give the year's

### Requirement: The sea evaporates by the bulk formula
Sea evaporation MUST follow E = C_E |U| (q_sat − q) on the air column, with the wind no
calmer than gusts allow, so that air crossing a narrow sea takes up only what the crossing
allows.

#### Scenario: Earth's own map
- **WHEN** the climate runs on Earth's topography with observed winds and sea temperatures
- **THEN** its rain's log pattern correlation with GPCP is at least 0.68, and Riyadh gets
  under 600 mm a year (GPCP 124)

### Requirement: Sinking air over land exports its vapour
Under large-scale sinking over land the air column MUST lose vapour to the divergent
outflow, at about 0.6 × 10⁻⁶ s⁻¹ under a desert's sinking.

#### Scenario: the Sahara on the model's winds
- **WHEN** the Earth climate suite runs at one degree on the model's own winds
- **THEN** Cairo and In Salah get under 105 and 87 mm a year, and the Sahara's biome is hot
  desert
