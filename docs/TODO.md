# Open items, in order

The working queue, highest payoff first. Detail and history live in `V2_ROADMAP.md` (section
named in brackets) and the scorecards in `V2_REALISM.md`. Accuracy comes first for now:
slower is fine. Speed is item 15.

## 1. Earth A/B in Rust, every layer against a real baseline [8f] — command done
Done: `geoforge compare REAL MODEL… --probe … --source …` writes the page. Every product
the real world holds as observed (authored) data gets a section: the observed map, then each
model's, one renderer and key; bias, RMSE and pattern correlation for quantities, direction
agreement and speed ratio for winds and currents, class agreement and shares for biomes and
Köppen (read from the observed climate); values at named stations, ranges and current boxes;
the numbers embedded as JSON. A model's copy of an observed input is skipped. `ab.sh` calls
it; `ab_page.py` is gone. Python stays only to turn observed files into world layers.

Still needed, each an observed layer for the real world (sources in `REFERENCE_DATA.md`):
- biomes (RESOLVE Ecoregions 2017), which the command already compares once imported;
- ocean currents (OSCAR);
- sea ice, as a product or view the real world can hold (NSIDC);
- hazards (IBTrACS cyclones, USGS quakes, Smithsonian volcanoes);
- the sky (HYG and the JPL ephemeris);
- rivers drawn against Natural Earth or HydroRIVERS, and the discharges at the mouths, which
  `check.py` still prints;
- the observed Köppen map itself (Beck et al.), beside the one read from observed climate.

## 2. Winds that feel the terrain [8e] — first step done
Done: the surface winds are a layer of stable air up to 3.5 km above sea level, thinning
over higher ground. A Helmholtz projection (`atmosphere/flow.rs`, `Keep::FlatLayer`) sends
the air round ranges and through gaps, keeping the rising and sinking the cells and monsoons
set over flat ground. The highs' rising and sinking in the moisture model is read from the
flat-ground winds. On the real-Earth chain, January direction agreement with ERA-Interim
went from 78% to 80% of the ocean and from 72% to 75% of land below 1.5 km.

Cost on the climate suite (now gaps with causes): Mongolia gets no rain and runs 4 °C warm,
the Pampas dry out, the Great Plains are drenched, the central Sahara gets rain (not yet
traced).

Still needed:
- an East Asian summer monsoon in the winds: July winds over China are easterly, so the
  turned westerlies dry the interior;
- surface friction over land, which slows the wind and turns it toward low pressure;
- lee troughs and cold-air damming;
- wind direction checked by region in `geoforge compare` (item 1).

## 3. Ocean currents that conserve mass [7] — done
- The gyres are drawn from the Sverdrup streamfunction on the cells' corners, with ψ = 0 on
  every corner touching land (`flow::from_streamfunction`): exactly divergence-free, nothing
  through a coast, and their east–west legs back.
- The jets (coastal, separating, overturning) are projected so none runs into a coast; only
  the overturning's sinking stays.
- The wind drift is the Ekman layer and keeps its divergence: water it piles against a coast
  sinks there, water it draws off one wells up.

On the real-Earth chain the SST error fell from 2.55 to 2.47 °C, the Kuroshio runs north
again (0.26 m/s), and the Labrador Sea is cold (−3 °C). Costs, now gaps: the Canary Current
runs 6 °C cold, the Bering Sea stays open, the central Arctic opens in summer.

Still needed:
- a test that the flow through the Drake Passage and Indonesia is continuous (ψ is 0 on every
  coast, so no gyre transport passes between continents; the circumpolar current is drift);
- the flip that stops the western boundary current past its separation breaks the gyre's
  closure locally: draw it into ψ instead.

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
