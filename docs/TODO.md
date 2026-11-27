# Open items, in order

The working queue, highest payoff first. Detail and history live in `V2_ROADMAP.md` (section
named in brackets) and the scorecards in `V2_REALISM.md`. Accuracy comes first for now:
slower is fine. Speed is item 15.

## Where to pick up (27 November 2026)
The last Earth A/B run (`scripts/real-earth/ab.sh`, page in `V2_REALISM` → *Rain*) puts the
largest remaining gaps in the model's own circulation. In order:
1. **The circulation's divergence** (item 5): the model's winds converge 2.5–3× as hard as
   ERA's at 850 hPa, and into the desert heat lows. Fixing it unlocks moisture flux
   convergence (item 7), the biggest rain gain measured (land rain correlation 0.65 → 0.74
   on observed winds), and with it the land/sea rain split and too many big rivers (item 11).
2. **A Walker circulation** (item 7): East Africa and the Arabian Sea rise where Earth's air
   sinks.
3. **Temperatures that follow** (item 4): dry subtropical land 5–9 °C cold, tropical
   highlands' seasons too wide, Lima warm, Reykjavik cold on the model's winds.
4. **Inland mountains** (item 10) and **river routing** (item 11).
5. **Bloomery all the way** (item 18), after `nix flake lock` and a first `bloomery check`.
Items further down that say "done" are kept as history.

## 18. Bloomery all the way [proposal 0002]
The flake is a Bloomery workspace (a33f6a1), but with its opinions turned down. Take them all:
- **Specifications**: move the OpenSpec changes (15 of them, ~150 requirements, never merged
  into a living spec) into `.bloomery/specs/`, as areas → features → requirement TOML in
  EARS form with design notes; tie the ~300 tests to requirement IDs with
  `#[cfg_attr(any(), bloomery("ID"))]`; map every crate to the feature it owns; turn on
  `checks.specifications` with `requireCrateOwnership` and untied tests denied. Retire
  OpenSpec and rewrite the workflow in `CLAUDE.md` and `AGENTS.md` around `bloomery spec`
  and `bloomery review`.
- **Layout**: `crates/{bins,libs}` → `packages/rust/{bins,libs}`, as Bloomery scaffolds it.
- **Code health**: no crate over 5000 production lines (`geoforge-cosmos` and
  `geoforge-surface` are ~8000 each; split them by layer), then work the soft limits
  (1500-line crates, cyclomatic 10, 250-token functions, duplication) down and commit the
  baseline.
- **Clippy**: no suppressions. Rework the functions behind `too_many_arguments` and
  `too_many_lines` (parameter structs, smaller steps) and empty the allow list.
- **Features**: `default-features = false` on every dependency, features named where used;
  turn `noDefaultFeatures` back on. Dependencies declared once in `[workspace.dependencies]`.
- **Toolchain and formatting**: nixpkgs' `rustc` only (drop `rust-toolchain.toml`); `nix fmt`
  (Bloomery's treefmt set) instead of `cargo fmt` and `rustfmt.toml`'s long lines.
- **Commands**: `bloomery check` is the gate (vulnerabilities, licences, bans, unused
  dependencies, lock freshness, per-crate tests); `just test`/`just check` call it.
- First: lock the flake (`nix flake lock`) and run `bloomery check` once as it stands, to see
  what each of these costs.

## 17. A visualization tool
An interactive viewer for a world:
- pan and zoom the layers;
- toggle keys and overlays;
- inspect a cell's values through every layer;
- step through the almanac's days.

## 12e. Compact objects, asteroids, fictional materials, nebulae, galaxy centres — done
(OpenSpec `compact-objects`, `asteroids-and-fictional-materials`, `nebulae-and-galaxy-centres`)
- Done: brown dwarfs (a quarter as many as stars), neutron stars and black holes in the
  starfield, their spins, horizons and feeding told, and their systems open in the viewer;
  belt members with spectral types, meteorite elements, porosity and families; a universe's
  fictional materials (`exotics.toml`) as deposits and markers only; nebulae of five kinds
  where their makers are (the just dead as stars of their own bins); galaxy centres (black
  holes, active nuclei, jets) and globular clusters; nebulae and clusters on the star chart.
- Also done since: belt sizes steepen past a knee (Sol's Kuiper belt has ~10 bodies over
  1000 km, Trojans ~20 over 100 km; `PowerLaw::with_knee`); M dwarfs' bolometric correction
  from nearby stars (they had come out 1–4 magnitudes too faint below ~3200 K).
- Next, in the viewer: a "next habitable system" jump; a body's kind told as icy, ocean or
  rocky (an icy moon now reads "rocky"); nebulae drawn on the star view's backdrop, not only
  listed; the galaxy centre's black hole drawn so it shows against the bulge.
- Next, in the cosmos: globular clusters' and the galaxy centre's stars generated one by
  one; an active nucleus's light in other galaxies' brightness and a quasar's in the sky;
  field core-collapse supernovae outside birth events counted once; fictional materials in
  a system's star or rings; nebulae given an address if they ever gain stored edits (see the
  change's design); the main belt's own knee near 100 km (its largest few are now slightly
  too many).
- Next: the climate items — winds, seasons, monsoons and the rain and temperature misses
  (2, 4, 5, 7): land runs wet (rivers carry 1.6× Earth's water to the sea) and deserts small —
  then the rest of mountains and rivers (10, 11), then speed (15).

## 12d. Temperate defaults, habitable systems, moons and rings — done (OpenSpec `habitable-moons-rings`)
- Done: a climate estimate before generating and habitability grades; default worlds are
  temperate; the viewer circles habitable systems and marks bodies; moons (regular, inner,
  captured, irregular, specks); rings on every giant; zodiacal light on the star chart.
- Next: moons of dwarf planets and belt members; resonant moon chains; narrow ringlets; a
  rocky planet's passing ring; the planet's own rings as an arc in its sky and their shadow;
  a "next habitable system" jump in the viewer.

## 12c. One universe per seed, and the viewer — done (OpenSpec `universe-addresses`)
- Done: addresses for every galaxy, star's system and body, generated from the seed and
  address alone; galaxy models from generated galaxies; default body by rule; universe
  directories (`new --seed/--at`, `focus`, `visit` in the universe); true colour for every
  body (`appearance` view); the viewer (`geoforge-viewer`), zoom universe → galaxy → stars → system →
  body, search, layers on a globe or map, make and generate from the page.
- Next: clouds over the true colour; stars' and galaxies' names searched farther out (an
  index of what has been looked at); the real Earth's true colour from its data; the
  viewer's galaxy view at the stars' own scale (sectors drawn as they come into view).

## 16. Names
Languages and cultures in `geoforge-names`.

## 12. Galaxy sectors and the night sky — done (OpenSpec `galaxy-sectors`)
- A galaxy of components (nucleus, bulge and bar, thin disk by age, arms, thick disk, halo,
  dust) and morphologies (Sa–Sd barred or not, S0, E, dSph, dIrr, starburst); the Milky Way
  from Bland-Hawthorn & Gerhard (2016) and Vallée (2015).
- Stars that live in the galaxy: nested cubes sized by how far their stars are seen, each
  cube generated from its seed, each star with an identity that seeds its system; young
  stars in clusters and associations of 100 M☉ and more; dust clouds where stars form.
- The sky read from them: neighbouring systems see each other; the real Earth's band runs
  through Sagittarius and Cygnus; the Sun's place holds 8,500–9,200 naked-eye stars (9,110
  observed), 18–26 brighter than 1.5 (22) and 320–360 systems within 10 pc (339), with the
  galaxy-integrated IMF (α = 2.6 above 1 M☉). Systems are drawn from every component by mass
  (a quarter of the Milky Way's in its bulge), with a place, a frame and a galaxy seed.
- Edits: `geoforge star add|move|change|remove|cluster`, `stars near|export|import`,
  `constellation add|remove|list` (drawn on the chart), `visit` (a new world at a star).
- Still open: a Gould-Belt-like local excess of young stars around the Sun (the real sky's
  young B stars); proper motions; a faster sky in a crowded core (cut to a brighter limit
  today); generated galaxy-recipe seeds now draw different systems than before (the
  realism tables' s-worlds predate this).

## 12b. Above one galaxy — first part done (OpenSpec `universe-and-figures`)
- Done: a universe of halos and galaxies in nested cubes, clustered, each galaxy with an
  identity; a world's galaxy placed in it (the Sun's is the real one); other galaxies in
  the night sky and its glow, rated by surface brightness; automatic constellations;
  `geoforge universe`, `geoforge constellation auto`.
- Open: visiting another galaxy's stars (a galaxy model scaled to a generated galaxy and
  its sectors, so the zoom runs universe → galaxy → star → planet); superclusters'
  filaments and voids beyond the lognormal field; a world holding more than one system.

## 12c. Sky and almanac refinements [6, 8b]
- A system's local environment: nebulae and supernova remnants (clusters and associations
  come with 12).
- Tides by date, proper motions, the Moon's periodic terms.
- A tidal range map for worlds with moons: the forcing of each moon and the star (mass and
  distance), the semidiurnal and diurnal parts, and the basins' response (resonance, wide
  shelves amplifying it, as the Bay of Fundy and the Severn), giving mean and spring tidal
  range along every coast; then the water level by date in the almanac.

## 10. Mountains [8e] — partly done (OpenSpec `basins-swells-and-dying-rivers`)
- Done: peaks held below a strength limit (9 km × 1/g); basins and swells on old stable
  crust; dynamic lift 1200 m. High ground averages 22% of land (Earth 27% with its ice) but
  varies 8–42% between worlds.
- Arc height and width from slab age and dip, not noise.
- A volcanic gap over a flat slab.
- Some sutures as sag basins rather than belts (the Amazon).
- Island arcs maturing to 2–2.5 km.
- Great escarpments (optional).
- Mountains away from plate edges: on PB2002 plates the Rockies average 190 m (observed
  1720), the Colorado Plateau 255 (2010), the Brazilian highlands 120 (770) and the East
  African plateau 600 (1260); land above 1 km is 3% of the surface (observed 7%). Uplift
  comes only at plate boundaries (flat-slab and intraplate uplift, plateaus).

## 11. Rivers and lakes — partly done (OpenSpec `basins-swells-and-dying-rivers`)
- Done: basins near Earth's sizes (the five largest drain 19% of land, Earth 14%; was 28%);
  dry basins close where their rivers die out (16.5% of land endorheic, Earth 18%; was 7%);
  Earth's reference for large rivers corrected to 16.
- Too many large rivers remain (21–39 over 10 000 m³/s on five worlds at 512): the climate
  gives land too much of the rain (item 7).
- Breach the lowest pass instead of filling.
- The wetland rule keys on discharge, so every big-river cell is wetland.
- The Great Basin drains to the sea.
- On Earth's own map (20 km ETOPO1) gorges block big rivers: the Congo carries 3 000–3 700
  m³/s (observed 41 000) and the Yangtze 2 500–3 000 (30 000) even on observed rain; the
  Nile 6 500–7 700 (2 800) and the Ob far too much.

## 14. Other worlds — done but for giants' storms (OpenSpec `systems-and-worlds`, `other-worlds-2`)
Done: brown dwarfs (a fifth of systems, cooling, M/L/T/Y); planets of wide companions;
thermal-tide rotation of locked worlds, thick air included (Venus's 117-day day, backward
half the time); binary evolution of close pairs (stable transfer, common envelope, mergers);
a companion's light on the other star's planets; ice shells (Europa-like and Ganymede-like
generators); locked oceans that superrotate and carry heat to the night side; fast spinners
without an ice cliff or rain lines at their fronts, and a stack of gyres per basin, one pair
per jet; giants' polar cyclones and hexagons; Titan's methane cycle (rivers, ~1% lakes near
the poles) and a dunes view for any world with rivers. Open:
- Storms that drift by date (a dated view of the bands).
- Titan's dune belt reaches ±50°, not ±30°: the climate gives no methane rain, so the dunes
  view sees only its lakes and rivers; Titan's pole-to-equator contrast is 12 K, not 3.
- Titan's seas lie in both hemispheres equally (the orbit's asymmetry is not modelled), and
  its rain does not follow the seasons.
- Binary evolution: a second transfer (cataclysmic variables, type Ia supernovae), donors of
  8 M☉ or more, and contact binaries are not modelled; the sky's stars still live as if
  alone.
- Mars's and Venus's dunes (sand from older water and wind abrasion).

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
- Dry subtropical land is cold: Riyadh 19.9 °C on the model's winds, 22.3 on observed winds
  (observed 26.5), Tehran 9.0 / 11.6 (18.5); on observed temperature and rain the
  derived layers get Riyadh right, so it is the climate's own temperature.
- Tropical highlands' seasons are too wide: the Brazilian highlands swing 14 °C (observed
  4), the Altiplano 19 (6). The Cerrado leaves the tropical classes for it (a known gap in
  `known_earth_climate`).
- Lima 28 °C on the model's winds (observed 19): no stratus deck, a weak Humboldt current.
- Reykjavik −4.6 °C on the model's winds (1.8 on observed winds, observed 4.7): too little
  Atlantic air reaches it.

## 5. Monsoons and storm tracks [8e]
- Done: the monsoon's pull follows a smooth heat low (three box blurs, close to a Gaussian,
  and a pull that eases off instead of a hard cap), so its winds no longer come in rectangles
  with straight edges; the summer inflow and winter outflow blend across the equator instead
  of meeting at a seam; an equinox or a locked world has no monsoon.
- The divergence the winds carry is 2.5–3× ERA's at 850 hPa and correlates only 0.3–0.5
  with it; the monsoon inflow converges into the desert heat lows. Fix this before the
  climate takes moisture convergence from the winds (item 7).
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
  almanac, a printable month page, and a faster drought record (20 years of daily weather
  take ~8 s on a 1024-wide Earth).
- The ENSO-like cycle comes from the widest equatorial sea, measured along the equator
  (real Earth: the Pacific, 18,000 km from 118°E to 80°W): its two periods are Earth's
  (2.5–5 and 4–7 years, drawn) scaled by that width, the delayed oscillator's basin-crossing
  time; it weakens below 12,000 km and is gone under 4,000 km or with no shore on the
  equator. Its warm phase wets the basin's east end and coast and dries its west end and
  the lands beyond either end (on Earth: Guayaquil wetter; Jakarta, Darwin and Manaus
  drier); elsewhere the sign stays a seeded pattern. Still to do: the far teleconnections
  (southern Africa, India's monsoon, the Horn of Africa) from the circulation rather than
  seeded, and the cycle's amplitude from the basin's thermocline tilt.
- Soil water by day — done: a 150 mm bucket with a snowpack (3 mm/°C/day melt) runs in the
  daily weather, started from two normal years; each day, the almanac, `weather`, `place`
  (as the window starts) and the report's weather tables give it, and the drought months
  read it. Finding: Moscow lies under 300 mm of snow water at the March equinox (580 in a
  snowy year), where 35 cm of late-winter snow (Climate of Moscow) is about 100 mm; the
  climate puts too much of Moscow's rain in winter (Jan–Mar 235 mm in that year, observed
  136) and too little in summer.

## 6b. Time types — done
- `geoforge-core::time`: `Instant` and `Span` in SI seconds; a world's `Clock` (its mean
  solar day and its year, constants); a `Calendar` trait with `Twelfths`, `Lunar` (a moon's
  lunations from the ephemeris's new moons) and `Custom` (a world's `calendar.json`: named
  months of whole days, a week, an era). Dates are whole days from the epoch, each in the
  year and month it begins in, so calendar years differ by a day (leap years) while the
  clock's year stays one orbit. OpenSpec change `time-and-calendars`.
- `geoforge almanac --calendar world|twelfths|lunar|lunar:<moon>`; the report dates
  fictional worlds' weather in their own calendar. Drought statistics stay on twelfths.
- Still to do (tasks 3.x of the change): weather, storms, sky and the ephemeris take
  `Instant` instead of Earth days, and commands take dates in a calendar.

## 7. Rain and temperature misses [8e, 8l, V2_REALISM]
- Done (OpenSpec `soil-water-and-bulk-evaporation`): land evaporates its soil water
  (Budyko), the sea by the bulk formula; storm-track eddies at 8 × 10⁵ m²/s; the column's
  capacity on 60% of the ground's swing; sinking air over land exports its vapour. Rain
  pattern correlation 0.59 → 0.68 on the model's winds, Köppen by group 61% → 68%; Siberia,
  Moscow and the Urals get their rain.
- Next, the largest gain in sight: **moisture flux convergence**. On observed winds it lifts
  the land's rain correlation 0.65 → 0.74 and the sea's 0.63 → 0.75; on the model's winds it
  drowns Arabia, because the circulation's divergence is 2.5–3× ERA's and converges into the
  desert heat lows (item 5). Fix the circulation's divergence first.
- The land gets too much of the rain (Earth A/B land 1070 mm, observed 860; generated
  worlds 1090 vs 780): the subtropical sea's column is 20–30% too moist, so the sea
  evaporates too little (the missing convergence again).
- No Walker circulation: the East African highlands (Nairobi 5200 mm, observed 840) and the
  Arabian Sea rise where Earth's air sinks.
- The western interior and Great Plains on the model's winds (the Rockies 1170 mm, observed
  410; Dodge City 1690 at one degree, observed 560).
- London: no rain shadow. Lima: no stratus deck, 7 °C warm.
- The Ob drains far too much (109 000 m³/s on the model's climate, observed 12 500): West
  Siberia is now wet, and its basin takes in too much.
- Antarctic plateau ~17 °C warm, with double the seasonal range.

## 8. Ocean leftovers [8e]
- No Indonesian Throughflow or Leeuwin Current: the seas off western Australia and the
  Timor Sea upwell.
- The Brazil Current runs ~3 °C warm; the Kuroshio is slightly warm.
- Upwelling that comes in one season (Somalia, Oman) cannot show in an annual anomaly.
- Cyclones ignore wind shear, so the South Atlantic gets them.
- On the model's winds the SST pattern correlates 0.49 with observed (0.97 on observed
  winds): the Canary current runs −5.6 K (observed −3.0).

## 9. Bugs seen on Thumi — done
- The straight-edged wedge in the ocean-floor age map was a plate with no ridge, aged from
  its own edges and jumping against its neighbours; its floor now carries on from the aged
  floor beside it (20 km/Myr half rate inward), and ages past 180 Ma are capped rather than
  restarted. The edges left are plate boundaries, where real ages jump too.
- The cyclone map's blocks were the 2° storm grid read cell by cell; it is now read between
  cell centres.

## 13. Star systems [5] — done (OpenSpec `systems-and-worlds`)
- Kepler-like compact chains (period ratios log-normal about 1.8, peas in a pod); packed
  giants that scatter by Petrovich's criterion, merging close in and ejecting far out
  (Safronov number); hot Jupiters by high-eccentricity migration, lonely.
- A rocky world's tilt wanders unless a big moon steadies it.
- Trojans where no heavier giant clears them.

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
