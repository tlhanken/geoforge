# Open items, in order

The working queue, highest payoff first. Detail and history live in `V2_ROADMAP.md` (section
named in brackets) and the scorecards in `V2_REALISM.md`. Accuracy comes first for now:
slower is fine. Speed is item 15.

## 1. Earth A/B in Rust, every layer against a real baseline [8f]
- `geoforge compare REAL MODEL`: paired maps with one key, metrics and the page, written by
  Rust like `geoforge report`. Python stays only for the one-time step that turns observed
  data (ERA, GPCP, ERSST, ETOPO) into world layers.
- Every section shows the real result beside the model, with all other layers set to real
  Earth. Missing baselines today:
  - biomes (a WWF ecoregions map);
  - ocean currents (OSCAR or drifter climatology);
  - sea ice (NSIDC);
  - hazards (IBTrACS cyclones, USGS quakes, Smithsonian volcanoes);
  - the sky (a real star catalogue and the JPL ephemeris);
  - rivers drawn on the map against HydroSHEDS, not only discharge numbers.

## 2. Winds that feel the terrain [8e]
Only the monsoon part of the wind is blocked by high ground. The cells, westerlies and
subtropical highs ignore it: on the Earth run, North America, Europe, Russia and Tibet show
no deflection.

Needed:
- surface friction over land, which slows the wind and turns it toward low pressure;
- flow steered around ranges and blocked by them (the Rockies, the Alps, the Tibetan
  Plateau splitting the westerlies);
- lee troughs and cold-air damming;
- a check of wind direction by region against ERA-Interim.

## 3. Ocean currents that conserve mass [7]
Currents are diagnostic sums: a Sverdrup streamfunction (divergence-free) plus wind drift,
equatorial currents, coastal jets, overturning and the circumpolar flow (not
divergence-free).

On the Earth run, 13% of coastal sea cells push more than 0.1 m/s into land (water piling
up).

Needed:
- project the total onto a divergence-free field with no flow through coasts (solve for a
  streamfunction on the sea);
- a test that the net flux through every closed coastline and strait is about zero, and
  that the flow through the Drake Passage and Indonesia is continuous.

## 4. A seasonal climate [8e]
One continentality serves both halves of the year. Needed: winter and summer computed apart,
with winter highs (Siberian, Canadian), surface inversions, and sea ice stored and read back.

It fixes:
- Yakutsk (−24 °C coldest month vs −38) and Siberia's dry winters, which make Dw where
  Df belongs;
- Moscow too cold, Utqiagvik too warm;
- no ice in the Okhotsk Sea or the Baltic;
- the frozen-sea effect on land, tried twice without it.

## 5. Monsoons and storm tracks [8e]
- The winter storm track sits too far equatorward: Phoenix, Tehran and the Great Basin are
  drenched.
- No Somali jet, and Tibet is not an elevated heat source, so Mumbai is dry, rain reaches
  Arabia and the West African monsoon goes too far north.

## 6. An almanac for a place, with droughts [8d]
- A year day by day at a place: weather, moons, sunrise and sunset, and the day's events.
- Slow rain anomalies (AR(1) plus an ENSO-like cycle), SPI drought index and soil moisture;
  drought frequency in `place`.

## 7. Rain and temperature misses [8e, V2_REALISM]
- London: no rain shadow.
- Land rain is short overall (Buenos Aires and Alice Springs get half their rain).
- Lima: no stratus deck, 7 °C warm.
- The East African highlands and the Kazakh steppe are too wet.
- Antarctic plateau ~17 °C warm, with double the seasonal range.

## 8. Ocean leftovers [8e]
- No Indonesian Throughflow or Leeuwin Current: the seas off western Australia and the
  Timor Sea upwell.
- The Brazil Current runs ~3 °C warm; the Kuroshio is slightly warm.
- Upwelling that comes in one season (Somalia, Oman) cannot show in an annual anomaly.
- Cyclones ignore wind shear, so the South Atlantic gets them.

## 9. Bugs seen on Thumi
- A straight-edged wedge in the ocean-floor age map.
- Blocky cyclone patches on the hazard map.

## 10. Mountains [8e]
- Arc height and width from slab age and dip, not noise.
- A volcanic gap over a flat slab.
- Some sutures as sag basins rather than belts (the Amazon).
- Island arcs maturing to 2–2.5 km.
- Great escarpments (optional).

## 11. Rivers and lakes
- Breach the lowest pass instead of filling.
- The wetland rule keys on discharge, so every big-river cell is wetland.
- The Great Basin drains to the sea.

## 12. Sky and almanac refinements [6, 8b]
- A system's local environment: clusters, nebulae, supernova remnants.
- Constellations, tides by date, proper motions, the Moon's periodic terms.

## 13. Star systems [5]
- Period ratios nearer Kepler's.
- Hot Jupiters.
- Obliquity stability with and without a big moon.
- Trojans.

## 14. Other worlds [Still open]
- Ice-shell worlds (Europa).
- Locked worlds carrying heat to the night side.
- Asynchronous rotation.
- Brown dwarfs; planets of companion stars.
- Giant storms.
- Streaks on fast rotators.

## 15. Speed: a lite generator
Once the accurate model is settled, any of:
- lighter algorithms;
- smaller grids (generation time grows ~4× per doubling of width: 2.8 s at 256, 4.2 s at
  512, 11 s at 1024, 41 s at 2048 on Thumi);
- models trained on the accurate generator's output.

Also:
- cache the ground and views so commands start in under a second;
- caps for derived views.

Resolution consistency is tested (`a_seed_is_the_same_world_at_any_resolution`): on
Thumi, land and sea at 256, 512 and 1024 agree with 2048 on 95%, 97% and 98% of cells.

## 16. Names
Languages and cultures in `geoforge-names`.

## 17. A visualization tool
An interactive viewer for a world:
- pan and zoom the layers;
- toggle keys and overlays;
- inspect a cell's values through every layer;
- step through the almanac's days.
