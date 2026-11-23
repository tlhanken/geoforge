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

## 4. A seasonal climate [8e] — temperature done
Done: the warmest and coldest months are computed apart (`climate::summer_winter`). Each
half-year's continentality comes from that season's winds; sea air reaches three quarters as
far inland in summer as in winter; a sea frozen in a season is barely maritime in it (the
frozen-sea effect, for land; the sea's own air keeps the open-water reckoning the sea-ice
view is set on); and a high-latitude interior gets the winter high's cold on top.

On observed winds (the A/B page), against the earlier climate:
- land temperature error 3.84 → 3.67 °C, range correlation 0.78 → 0.81, Köppen 36.5% →
  38.4% of land in full class, 62.4% → 63.5% by group;
- the range north of 50°N is within 1 °C (it was 6 °C short at 60–80°N);
- Moscow 5.5 °C (observed 5.4, was 2.7), Yakutsk −7.9 °C with a 50 °C range (observed −8.7
  and 57; was −5.3 and 46), Utqiagvik −9.7 °C (observed −11).

On the model's own winds the climate suite went 190 → 185: its January winds bring little
Atlantic air to Moscow and its July winds bring sea air into the Sahara.

Still needed:
- sea ice stored as a product and read back by the climate, not a view set on its numbers;
- surface inversions in basins (Yakutsk, Verkhoyansk, Fairbanks);
- Siberia's dry winters (Dw where Df belongs);
- the Southern Hemisphere's continental ranges, 5–11 °C too wide at 20–50°S (Australia,
  southern Africa, Argentina; before this change too);
- Greenland's top, whose range is 51 °C on observed winds (observed 27).

## 5. Monsoons and storm tracks [8e]
- Done: the monsoon's pull follows a smooth heat low (three box blurs, close to a Gaussian,
  and a pull that eases off instead of a hard cap), so its winds no longer come in rectangles
  with straight edges; the summer inflow and winter outflow blend across the equator instead
  of meeting at a seam; an equinox or a locked world has no monsoon.
- The winter storm track: on observed winds the climate gives Phoenix 233 mm, the Great Basin
  436 and Tehran 339 (GPCP 344, 314, 446); on the model's winds 468, 626 and 548. The
  model's January winds over the dry subtropical land are 5 m/s westerlies where ERA has
  1 m/s or calm, carrying Pacific and Mediterranean air inland.
- Slowing the wind over land by 40% for friction cut rain over every interior (the Congo
  turned savanna, Alice Springs got 55 mm) and turning it 15° toward low pressure moved the
  directions off ERA's: the moisture transport is set on frictionless winds, so friction
  needs the moisture model retuned with it.
- Raising the monsoon pull's ceiling from 1 to 2.5 strengthened the July winds over the
  Arabian Sea (+5.8 m/s at 15°N 65°E) but not into westerlies at 8°N 55°E, overdid the Sahel
  (+10 m/s) and lowered the rain correlation (0.59 → 0.57): the jet needs the
  cross-equatorial flow turned east and held against the East African highlands, not a
  stronger pull.
- No Somali jet: in July ERA has 16–18 m/s westerlies over the Arabian Sea (8°N 55°E,
  15°N 65°E), the model 0.5–3; its July winds over southern China are easterly (−7.5 m/s,
  ERA +2 from the south). Arabia gets 732 mm at Riyadh on the model's winds (513 on observed
  winds, GPCP 124), and the Sahel 960 mm at 12–16°N (550 on observed winds, GPCP 550).

## 6. An almanac for a place, with droughts [8d] — done
- `geoforge almanac`: a year day by day at a place (sunrise and sunset, moons, weather,
  events), its months' rain, SPI-3 and soil water; `place` reports drought frequency.
- Slow rain anomalies (regional AR(1) plus an ENSO-like cycle) in all daily weather.
- Still needed: SPI-6/12, drought onset and end, eclipses and meteor showers in the
  almanac, a printable month page, tying the cycle to an ocean basin, and a faster drought
  record (20 years of daily weather take ~8 s on a 1024-wide Earth).

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

## 9. Bugs seen on Thumi — done
- The straight-edged wedge in the ocean-floor age map was a plate with no ridge, aged from
  its own edges and jumping against its neighbours; its floor now carries on from the aged
  floor beside it (20 km/Myr half rate inward), and ages past 180 Ma are capped rather than
  restarted. The edges left are plate boundaries, where real ages jump too.
- The cyclone map's blocks were the 2° storm grid read cell by cell; it is now read between
  cell centres.

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
