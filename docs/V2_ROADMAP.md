# Geoforge roadmap

What comes next, in order, and why. Sources: our own gaps (`V2_REALISM.md`), and a
comparison with Artifexian's *Worldbuilder's Log* series (Logs 0–48), from two NotebookLM
summaries of it (fact-checked below). Each item says what it adds, what it reads, what it
writes, and which Logs it follows.

Status: **done**, **next** (agreed, not started), **later** (agreed, not scheduled), **open**
(a gap in what exists).

## Principles

- **Bottom-up causation**, as now: star → system → body → tectonics → geology → elevation
  → circulation → climate → erosion → hydrology → biomes → (new) resources, hazards.
- **Stored products are things you can edit.** A classification or summary that is a pure
  function of stored products (Köppen climate, sea-level pressure) is a **derived view**:
  computed on export or render, never stored, never stale, never edited.
- **Checked against Earth** where Earth is the answer, as the known-case suites do now.

## 1. Climate at a glance: Köppen–Geiger view — done

The `koppen` derived view: 30 classes by Peel et al. (2007)'s rules, worked out from
`temperature`, `temperature_range`, `precipitation` and `precipitation_warm_share` (see
`V2_PIPELINE.md`). The C/D boundary is a parameter, 0 °C by default as in Artifexian's
Logs 36–40. Months are rebuilt as sinusoids, so two rainy seasons and short tropical dry
seasons are smoothed out. Checked against Earth in `known_earth_climate`: Köppen groups, full
classes, and the classifier fed station normals (`V2_REALISM.md`, Köppen view).

Artifexian (Log 40) decides B against D from a lookup table, because his rain is
qualitative. We have rain in mm, so Köppen's own threshold applies directly.

## 2. Derived views — done

Done: the mechanism and the January and July sea-level pressure views (`pressure_jan`,
`pressure_jul`), recovered from the winds plus land domes; 10 of 12 ERA5 centres within
8 hPa (`V2_REALISM.md`, Pressure). The plan follows.

`views::View` names the products a view reads; `render` and `export` take view names
beside products (`render <dir> koppen`). Categorical views export an exact-colour PNG, a
legend with land shares, and world files.

**Sea-level pressure (Jan, Jul), hPa**, from the circulation pattern and land heating:
the equatorial trough (~1008–1010), subtropical highs (~1020–1025, strongest over cool
ocean in summer), subpolar lows (~990–1000), polar highs; winter highs over cold
continents (Siberia ~1035) and summer heat lows (~1000–1005). Rendered with isobars and
H/L marks. Check: surface winds should cross isobars toward low pressure at 20–40° (as
they do in the stored winds), and against ERA5's seasonal pressure maps. Artifexian (Logs
29–32) draws pressure first and the winds from it: summer subtropical highs at ~35°,
winter ones at ~25°, over the oceans toward the cold eastern currents. Continental winter
highs and summer heat lows sit in the interiors, subpolar lows at ~60° in winter, and the
ITCZ hugs the equator at sea but swings poleward over hot land. Our winds come from the
cells, so the pressure view is drawn from them, and should land in the same places.

## 3. Resources and soils — done; refinements open

Done as derived views (`geoforge_analysis::resources`; see `V2_PIPELINE.md`): prospectivity
for iron, copper, tin, lead–zinc–silver, gold, coal, salt, peat and gemstones, soils with
fertility, and timber; `place` reports what is within reach of a point and draws how many
deposits a local map holds. Next:

- Done since: the **giant-deposit catalog** (the resources layer, `deposits`), the one
  stored and editable part.
- An Earth check against USGS deposit data on the ETOPO-forced run (porphyry belts on the
  Andes, banded iron on the Pilbara, salt in the Permian and Zechstein basins).
- Loess, mercury, sulfur, obsidian, flint, clay, building stone, aquifers; organic
  resources (fisheries on upwelling coasts, furs, fibre plants, guano).
- Peat on ground glaciated in the last ice age (needs an ice-age extent).

The original plan follows.

Global maps cannot place individual mines, so the layer stores
**prospectivity**. For each deposit kind it keeps an expected count per 10,000 km², λ,
from the geology, climate and history already stored. Exported, that becomes a grayscale
map of P = 1 − e^(−λ), the chance a region holds at least one deposit (0–255). When a
region is zoomed into for local generation, deposits are drawn there as a Poisson process,
seeded by place, so every zoom agrees. A few giant deposits per kind (the Pilbaras,
Potosís, Zechsteins) are drawn globally and stored as a catalog: kind, place, grade,
size, depth, age. Fields: soil class and fertility, timber, fishery productivity.

Metals:

| Resource | Where (deposit types) | Reads |
|---|---|---|
| Iron | banded iron formations in Archean–Paleoproterozoic cratons (2.5–1.8 Ga); bog iron in cold wet lowlands; laterite crusts over mafic rock in tropical wet–dry climates; oolitic ironstones on old shelves; magnetite skarns | lithology, rock age, provinces, climate, wetlands |
| Copper | porphyry copper in arcs (most of the world's copper); skarns where arc intrusions meet limestone; volcanogenic massive sulfides in ophiolites and arcs; native copper in flood basalts (Keweenaw); sediment-hosted copper in rift basins | provinces, lithology |
| Tin | greisens and pegmatites in collisional (S-type) granites; placers downstream | provinces, lithology, rivers |
| Lead, zinc, silver | Mississippi Valley-type in carbonate platforms; sedimentary-exhalative in rift basins; epithermal veins in arcs (Potosí silver) | provinces, lithology |
| Gold | orogenic gold in greenstone belts and old orogens; epithermal and porphyry gold in arcs; placers in rivers below them | provinces, rock age, rivers |
| Also | mercury (cinnabar, arcs and hot springs: gold amalgamation); arsenic (arsenical bronze); nickel, chromium, platinum (ophiolites, layered intrusions) | |

Fuel: timber (by biome, below the tree line: hardwood/softwood), charcoal from it; peat
(cold, wet, flat, poorly drained, often post-glacial); coal (swamp basins of the last
~360 Myr, since land plants, kept in foreland and intracratonic basins); bitumen and oil
seeps in sedimentary basins.

Salt and chemicals: halite (marine evaporites in old restricted basins and rifts, salt
domes, playas in closed basins, sea-salt pans on arid coasts); gypsum; natron and soda
(alkaline lakes: glass); saltpetre (deserts, caves: gunpowder); sulfur (volcanoes);
alum (dye mordant).

Stone and clay: building stone (limestone, granite, marble, sandstone); flint and chert
(chalk, limestone: toolstone); obsidian (young rhyolitic volcanics); clay (floodplains:
pottery, brick); glass sand; gemstones (pegmatites; diamonds in kimberlites in cratons).

Water: aquifers (sandstone, karst), springs, oases.

**Soils** (fertility for plants and farming), from parent rock, climate, relief and age:
young volcanic soils (very fertile), loess (wind-blown silt downwind of deserts and ice
margins), floodplain alluvium (renewed by floods), black grassland soils (chernozems:
breadbaskets), leached tropical laterites (poor), boreal podzols (acidic, poor), arid
soils (saline), peat soils. Soils should feed back into the biome layer's vegetation.

**Organic resources**: fisheries (upwelling coasts, shelves, cold currents, salmon rivers),
whales and seals, grazing and game (grasslands), furs (boreal), timber species, fibre
plants (flax, cotton, hemp), wild grains (the domestication argument: Mediterranean,
winter-rain climates), spices and incense (tropical and arid), dyes (murex, indigo, woad),
guano (arid islands under upwelling: fertiliser), amber (old conifer coasts), pearls
(warm shallow seas), reeds and papyrus (wetlands).

Needs a little history the pipeline does not have yet (when land plants arose, which
basins were swamps or restricted seas); until item 9, ages come from province and rock
ages.

**Artifexian's placement rules (Logs 45–48)**, as NotebookLM gives them. They key on what
we store, and on deep-time windows that must come from this world's own history, not
Earth's dates:

| Resource | Rule | Our inputs |
|---|---|---|
| Peat | subarctic, and non-dry-summer continental and temperate climates, on ground glaciated at the last glacial maximum | Köppen view; an ice-age extent (not stored: from today's ice and a cooling) |
| Timber | below a tree line of 4000 m to 30° latitude, falling 130 m per degree to 50°, then 75 m per degree to sea level (~69°); none in B, ET or EF | elevation, Köppen view |
| Coal | wet lowland basins at 0–45° latitude in a "Carboniferous" window (340–280 Ma on Earth), near eroding mountains. High grade where later folded into mountains, low grade in flat sandstone and carbonate. Destroyed under later flood basalts | rock age, lithology, provinces, LIPs |
| Copper | basements of young and old orogens. Native copper in carbonates around flood basalts and hotspots (mostly poor). Enriched VMS in collisional belts, as sediment beside andesite and rhyolite, only where winters stay above 6 °C (snow stops the enrichment). Arsenical copper where copper belts cross coal basins and carbonates. Sedex in old back-arc basins later collided. Skarns where intrusions meet limestone. IOCG in rhyolitic cratonic shields | provinces, lithology, rock age, climate |
| Tin | rhyolites of collisional orogens; placers downstream; natural tin-bronze where tin belts cross copper zones | provinces, lithology, rivers |
| Iron | bog iron in wetlands downhill of metal belts (bogs in subarctic climates; riverbanks and marshes in wet tropical, subtropical and oceanic ones); skarns; laterites on basalt plateaus in tropical wet–dry climates; oolitic ironstone in warm shallow seas at 0–45° (290–120 Ma on Earth) fed by mountains; meteoric iron as rare single finds in deserts and tundra | wetlands, provinces, lithology, Köppen view |
| Salt | sea salt on every shore (sun-dried where arid, boiled elsewhere); salt flats in arid basins with clastic rock (up to 10,000–43,000 km²); rock salt from evaporites laid at 20–33° in warm seas (290–120 Ma on Earth), later uplifted into domes | coasts, Köppen view, closed basins, lithology |

NotebookLM mixes up fahlore with native copper here too; where its rules disagree with the
deposit geology above (porphyry copper, the largest source, is missing from them), the
geology wins.

Check: Earth's own deposits by type and setting (USGS mineral resource data) against the
ETOPO-forced Earth run: do porphyry belts sit on the Andes and Rockies, tin on the
Southeast Asian and Cornish granites, banded iron on the Pilbara and Minas Gerais cratons,
salt in the Permian and Zechstein basins?

## 4. Hazards — done as views; refinements open

Done as derived views (`geoforge_surface::hazards`): earthquakes (with magnitudes),
tsunamis, volcanoes, cyclones, tornadoes, floods, wildfire, landslides, blizzards, dust
storms, fog, auroras and lightning tiers; `place` gives each one's chance over a window,
its season counted, and the chance of magnitude 7 and 8 quakes. Checked on Earth's plates
and climate (`V2_REALISM.md`). Open: wind shear for cyclones (the South Atlantic), storm
intensity, droughts, volcanic explosivity, weather on a given day (item 6).

The original plan follows.

A new layer after hydrology. Rates and seasons per cell (or region), plus a **date and
region calculator**: the chance of each hazard in a window, P = 1 − e^(−λT).

- **Earthquakes**: Gutenberg–Richter, log₁₀ N(≥M) = a − bM with b ≈ 1, the rate set by each
  boundary's slip rate × length × seismogenic depth (seismic moment, Kostrov 1974), and a
  corner magnitude by boundary type (Bird & Kagan 2004: subduction zones M9+,
  collisions ~M8+, transforms ~M8, ridges rarely above M7), with intraplate quakes rare.
  Reads plates, boundaries, provinces. Tsunamis from subduction coasts.
- **Volcanoes**: arcs, hotspots, rifts; eruption rates and styles (explosive arcs,
  effusive hotspots).
- **Tropical cyclones**: basins where sea surface exceeds 26.5 °C, 5–30° from the equator,
  in the warm season; tracks west with the trades, curving poleward. Reads SST (climate),
  winds.
- **Severe storms and tornado alleys**: flat mid-latitude plains where warm moist air
  from a tropical sea meets cold polar air under dry air off a plateau.
- **Others**: river floods (monsoon, snowmelt), droughts (rain variability), wildfire
  (Mediterranean and savanna), landslides (steep, wet, seismic), blizzards, dust storms,
  fog coasts (cold upwelling), lightning (convective rain).
- **Auroras**: from the magnetic field model (done: `geoforge_core::magnetism`): ovals
  where field lines reach the magnetopause (Earth 71°), around a pole tilted from the spin
  axis; none where the dynamo has died (Mars, Venus).

Artifexian's numbers (Logs 42–43), as NotebookLM gives them:

- **Lightning**: four zones of influence: warm-month ≥ 18 °C (not deserts or polar), big
  landmasses, low pressure, windward slopes. Under 1 flash/km²/yr outside them, 1–5 in
  any, 5–15 where warm moist air blows onshore into one, over 15 where three or four
  meet; none poleward of the 18 °C isotherm; one tier less in rain shadows.
- **Tornado alleys**: summer only, 30–50° latitude, flat plains where dry air off a
  plateau meets warm moist sea air under low pressure.
- **Auroras**: the magnetic pole up to ~10° from the spin pole; the main oval 10–20° from
  it, a sub-auroral band 20–30°, storm-only displays 30–40°.
- **Fog**: upwelling and cold-water coasts, as winter, summer or year-round fog.
- **Tropical cyclones**: summer, over sea above 27 °C, not within 5° of the equator; form
  in the east of warm oceans, run west with the trades, recurve poleward, and die past the
  18 °C sea isotherm; none in basins fed by polar currents (the South Atlantic).
- **Winds**: gale corridors where open-ocean fetch hits coasts beyond 30°; katabatic
  winds off ice sheets and high plateaus; föhn on the lee of ranges over 2000 m.
- **Blizzards** in D and E climates in winter with strong wind; **sandstorms** on flat B
  plains in summer; **dust storms** on dry flat plains outside B and E in summer.

Check: Earth's earthquake catalog (USGS) by region and magnitude, IBTrACS cyclone tracks,
the Smithsonian volcano list, against the ETOPO-forced Earth run.

## 5. Star system fixes — mostly done

Done (`geoforge_cosmos::dynamics`):

- **Orbit spacing by Hill radii.** Neighbours are kept at least 8 mutual Hill radii apart,
  ((m₁+m₂)/3M)^(1/3)·(a₁+a₂)/2 (Chambers et al. 1996; Jupiter–Saturn sit at ~8; two
  planets alone are safe beyond 3.5, Gladman 1993). Too close, the outer planet moves out;
  in the habitable recipe the focus stays put and its neighbours give way.
- **Unstable systems** are not simulated; what is left is. Pairs closer than 5 scatter:
  the lighter planet is thrown out and the survivor turns eccentric (e 0.1–0.5). Orbits
  widen as a star loses mass but Hill spheres widen more, so white dwarf systems resettle
  (some scatter, as rock falling onto white dwarfs suggests).
- **Belts set by resonances**: the asteroid belt between the first giant's 4:1 and 2:1
  (Jupiter's 2.06–3.28 AU), the Kuiper belt between the outermost giant's 3:2 and 2:1
  (Neptune's 39.4–47.8 AU); a = a_giant·(P/P_giant)^(2/3).
- **Moons** sit between the Roche limit and 0.4 of the planet's Hill radius (prograde
  orbits are stable to ~0.5, Domingos et al. 2006); a big moon drawn beyond is brought in.
- **`show`** prints breathability (O₂ 0.16–0.5 bar, CO₂ under 0.005 bar, N₂ under 3 bar,
  gravity) and tides from every moon and the star, in Earth's lunar tides, with spring
  and neap ranges.

Also done (OpenSpec `systems-and-worlds`): small planets in Kepler-like chains (period
ratios log-normal about 1.8; masses alike within a system, Weiss et al. 2018); giants born
packed between the 3:2 and 2:1 resonances, scattering when Petrovich's (2015) criterion
fails, merging when the Safronov number is under 1 and ejecting otherwise (survivors'
eccentricities Rayleigh 0.3, Jurić & Tremaine 2008); hot Jupiters from 7% of giant
scatterings by tidal circularisation at twice a periastron of 1–3 Roche radii, clearing the
inner system; tilts that wander (10–45° from birth) unless a moon's pull outweighs the
star's; Trojans unless a heavier giant lies within 15 mutual Hill radii. Checked: giants'
median e 0.17 (Kipping 2013), hot Jupiters 0.8% of Sun-like stars, 1.1 small planets within
200 days per Sun-like star and 1.8 per red dwarf.

Also done: the `positions` command and export (every star, planet and moon's place at a
date, by Kepler orbits), and close binaries now share one orbit about their barycentre
(each star had drawn its own tilt and phase).

## 6. Sky, calendar and almanac — done; refinements open

Done: the focus's solar day, year and synodic months (`show`), its local date and moons'
phases (`positions`), and for any date (`sky`, `weather`):

- its equatorial frame (the present is its northern spring equinox), so the star's
  declination through the year drives the seasons, eccentric orbits included;
- **eclipses** from the ephemeris: solar (total, annular, partial) and lunar, checked on
  Sun–Earth–Moon (2.3 solar and 1.4 lunar a year against NASA's 2.38 and 1.45);
- a **star map**: the stellar neighbourhood drawn from the same galaxy population at the
  system's radius, apparent magnitudes and colours, the galaxy's band with a dust lane; a
  chart and a catalogue; the system's stars, planets and moons on the date;
- **weather on that day** at a place: seeded, adding back up to the climate's means,
  rain and seasons.

Also done: a global storm field (lows with fronts, tropical cyclones) that the daily
weather follows, and a weather map for any date (`weather-map`).

Also done: **the real Solar System** for archetype worlds (`geoforge_cosmos::sol`): the
Sun, eight planets and Pluto from JPL's approximate-position elements, the Moon from
Meeus's mean elements with its node and perigee turning, 20 major moons from JPL's mean
elements, the main and Kuiper belts and Saturn's rings, set at the body's northern spring
equinox (Earth: 2026-03-20 14:46 UTC) with real dates in `show` and `sky`. Checked against
the 2026–27 almanac (`known_solar_system`): the five eclipses from February 2026 to
February 2027 with the right kinds, within 4 hours (the Moon's periodic terms are left
out); six conjunctions and oppositions of 2026 within 2 hours; Earth's equinox from the
elements within 13 minutes and Mars Year 38's start within half a day. The neighbourhood
stars are still drawn from the galaxy model, not a catalogue.

Open: constellations (named asterisms from the brightest stars); sunrise, sunset and day
length by place; tides by date; the Moon's periodic terms (evection, variation) for
eclipses to the minute; drawn moons' precession; seasonal phenomena on a date (cyclone season, monsoon onset, meteor showers); a check
that no supernova progenitor lies within ~10 pc; the Gould Belt (the Sun's sky has 9,110
naked-eye stars, a smooth disk gives ~6,500: its young B and A stars are a local
excess); proper motions.

Artifexian (Log 6): tides scale as the moon's mass over its distance cubed, plus the
star's (spring and neap tides). The calendar comes from the year (Kepler), the day, and the
synodic month, 1/P_syn = 1/P_orbit − 1/P_year. His neighbourhood (Log 3) is a 12-light-year
sphere at 0.004 stars per cubic light year (34 stars, 20 systems), with Raghavan et al.
(2010)'s multiplicity for Sun-like stars (56% single, 33% binary, 8% triple, 3% more). We
draw 0.09 stars per cubic parsec (0.0026 per cubic light year), the RECONS 10-pc census of
stars and white dwarfs.

## 7. Oceans — mostly done

Done as derived views (`geoforge_surface::oceans`): sea ice, coral reefs and atolls,
productivity, checked on the Earth run (`V2_REALISM.md`, Oceans). Open: deep-water
formation; sea ice feeding back on albedo and coastal climate (it would need to be stored
and read by the climate layer); reefs on finer grids than the shelf shows. The plan
follows.

- **Sea ice**: winter and summer extent from SST and air temperature; feeds albedo and
  coastal climate.
- **Deep-water formation** (thermohaline): where cold salty water sinks (subpolar seas
  by sea ice), and its effect on the warm currents.
- **Coral reefs and atolls**: water 18–29 °C, under ~50 m, clear (away from big river
  mouths); atolls on subsiding volcanic islands of hotspot tracks.
- **Upwelling and productivity** as a map (fisheries, guano, fog deserts).
- Check: NSIDC sea-ice extent, Allen Coral Atlas, against the ETOPO-forced Earth run.

## 8. Editing near the poles — done

Edit views are equirectangular PNGs with world files, which QGIS reads and Blender maps
onto a UV sphere as is (Artifexian's workflow: paint on the sphere, save the
equirectangular image, import it back). Done: rasters export at N times the grid's size and
import at any 2:1 size; optional polar cap views (azimuthal equidistant, poleward of 60°,
georeferenced) merge back where painted; a Blender round-trip guide (`V2_PIPELINE.md`,
Editing a layer). Open: caps for derived views; a cap edge other than 60°.

## 8a. Weather and hazards as one — done

`place` counts tropical cyclones from the storm field's tracks (and the cyclone map is their
track density); the daily weather draws thunderstorms, tornadoes, floods, landslides,
blizzards, wildfires, dust, fog, quakes, eruptions, tsunamis and auroras from the same
yearly rates, on the days that favour them. Open: lows merging or repelling rather than
piling up; hail; storm surge with cyclones; magnitudes for drawn quakes.

## Refactoring — done

- `geoforge-almanac` holds everything in time: ephemeris and calendar (from core), eclipses,
  the sky, the storm field, daily weather and hazards near a place (from the pipeline);
  `geoforge-pipeline` is orchestration, views and the place report.
- `geoforge-analysis` holds the derived readings (Köppen, ground, hazards, oceans,
  pressure, resources); `geoforge-surface` the generators of stored layers.
- Store: `render/{mod, maps, sky, weather}`, `edit/{mod, import, patch}`, `cache`.
- CLI: a module per command family (`world`, `edit`, `system`, `almanac`), full help.
- v1 (`crates/bins/geoforge`, examples) removed; the binary is `geoforge`.
- `rustfmt.toml` matches the codebase's long-line style; a `justfile`; `GETTING_STARTED.md`.

Open: cache the ground and views so commands start in under a second.

## 8c. Names — crate done; languages later

Naming lives in `geoforge-names`, a crate with no other dependency: `name(seed)` gives the
same names as before (one sound set: 16 onsets, 8 nuclei, 2–3 syllables). Later:
phonologies — consonant and vowel inventories, syllable shapes, stress, spelling rules, as
in Artifexian's conlang videos — seeded per world, culture or region, so stars, planets,
storms and places each get a consistent language.

## 8d. An almanac for a place — first version done

A calendar for a whole year at a place, worked out day by day on its own calendar: the
date (local day, month by its moons), sunrise and sunset and day length, the moons'
phases, the day's weather (high, low, rain or snow, sun or cloud, wind, the storm
passing), and the natural phenomena of that day: eclipses, auroras, thunderstorms and
tornadoes, cyclones, earthquakes and eruptions drawn from their rates, meteor showers.
Output as a table (JSON/CSV) and a printable calendar page per month.

**Droughts** (and wet spells) belong here too, for harvests and societies later: the
weather today has no year-to-year swings beyond the storms, so a dry year never comes.
Add slow rain anomalies — a seeded monthly-to-yearly process (AR(1) over months, plus an
ENSO-like oscillation of a few years where the ocean allows) that scales the storms' rain
and the day's odds — then track a standardised precipitation index (SPI over 3, 6 and 12
months) and soil moisture (rain against evaporation from temperature), report drought
onset, severity and length in the almanac, and summarise drought frequency in `place`
for crop risk.

Done (November 2026): `geoforge almanac DIR --lat --lon [--year N]`
(`almanac/src/year.rs`, `drought.rs`): each local day's sunrise and sunset (the star's
declination from the season, 0.83° for the sun's radius and refraction; polar night and
midnight sun), each moon's lit share, the weather and the day's events; then the year's
months — rain summed from those days, SPI-3 read against the 20 years before, soil water
— and the place's droughts. Rain now has slow wet and dry years: a factor whose logarithm
is a regional AR(1) over months (0.6 a month, on a 15° lattice so neighbours share it)
plus an ENSO-like cycle of two seeded periods of 2.5–7 years, felt in the tropics with a
regional sign; spread 0.2 (wet) to 0.5 (desert). SPI-3 is a log-normal fit per calendar
month over the place's own daily weather of the 20 years before, so it is set on the same
storms as the days it reads (a monthly stand-in without the storms read real years at a
mean SPI of −1.0 with a spread of 2.1; on the daily record six Thumi years read 0.2 and
1.0); about 16% of months read moderate and 7–8% severe, as a standardised index should.
`place` reports drought frequency and the index as its window starts (the record costs
about 8 s on a 1024-wide Earth).

Then (November 2026): the cycle comes from the world's widest equatorial sea (periods
scaled by its width over the Pacific's 18,000 km, the delayed oscillator's crossing time;
weak under 12,000 km, gone under 4,000 or with no equatorial shore), its warm phase wetting
the basin's east and drying its west and the lands beyond (Guayaquil wetter; Jakarta,
Darwin, Manaus drier on real Earth). The soil bucket runs by day with a snowpack, so every
day, `place` and the almanac give the soil's water and the snow lying. Time has types
(`geoforge-core::time`: `Instant`, `Span`, a world's `Clock`, and calendars — twelfths,
lunar, the world's own `calendar.json`; OpenSpec change `time-and-calendars`), and the
almanac dates its days in any of them (`--calendar`).

Open: SPI-6 and -12, onset and end of each drought, eclipses and meteor showers in the
almanac, a printable month page, the far teleconnections from the circulation, and the
internals (weather, storms, sky) on `Instant`.

## 8a. Galaxy sectors — done (November 2026)

OpenSpec change `galaxy-sectors`. A galaxy is a list of components
(`geoforge_cosmos::galaxy`): shapes (exponential disks whose stars thicken with age,
exponential ellipsoids for bulges and bars, Hernquist spheroids, Plummer nuclei, flattened
broken power-law halos), histories and metallicity relations, logarithmic arms that gather
young stars, and dust. Stars are generated where they are looked at
(`geoforge_cosmos::sectors`): each brightness bin in cubes half as wide as it is seen from,
each cube a pure function of the galaxy's seed and its key, each star an identity that seeds
its system; young stars come from birth events (clusters and associations of 100 M☉ and more,
a twentieth bound, the rest spreading at 1.5 pc/Myr), dust from a smooth layer and molecular
clouds placed with them. The sky reads them through the dust (`sky`, `starfield` with the
world's edits). Calibration: the galaxy-integrated IMF (α = 2.6 above 1 M☉; with the canonical
2.3 the Sun's sky held 11,400 stars and twice the third- and fourth-magnitude ones) and the
thin disk's mass set by the 10 pc census. A sky takes under a second on two cores.

## 8h. One universe per seed, and the viewer — done (November 2026)

OpenSpec change `universe-addresses`. Every item of a universe has an address
(`geoforge_core::address`) and is generated from the universe's seed and that address alone
(`geoforge_cosmos::resolve`): a galaxy (generated, with a model of its own — its type's preset
scaled to its mass and size), the system born with one of its stars, a body of that system.
A universe's default body is an address chosen by a fixed rule (its galaxy most like the Milky
Way nearest its origin; the nearest Sun-like star with a wet world near that galaxy's
reference circle; its focus). A universe directory keeps what is made or edited
(`universe.json`, shared `galaxies/<galaxy>/` documents, a world per body at
`bodies/<galaxy>/<star>/<body>/`). Every body has a true-colour face (`appearance`: a generated
surface's seas, ice, land cover and snow; a body without one drawn from what it is; stars by
temperature). The viewer, `geoforge-viewer DIR` (its own program and crate: scenes, search and maps in
Rust behind a local server; the page draws them with three.js).

## 8g. The universe and the sky's figures — done (November 2026)

OpenSpec change `universe-and-figures`. A universe from a seed
(`geoforge_cosmos::universe`): dark-matter halos by the z = 0 mass function in thirteen
half-decade bins, each bin in cubes of its own size, clustered by a lognormal field with
mass-dependent bias; in each halo a central (Moster 2013 stellar masses, types by mass and
environment, Shen 2003 sizes) and satellites in subhalos, heavy ones sunk by friction and
inner ones torn apart by tides. Each galaxy has an identity (`h{bin}:{i},{j},{k}:{n}.{m}`)
that regenerates it. A world's galaxy is the generated galaxy like it nearest the origin;
the Sun's universe is the real one (the Local Group and the bright galaxies to magnitude
7.5). The sky lists galaxies to magnitude 8 with how they are seen (by magnitude and
central surface brightness) and adds them to the glow through the home galaxy's dust.
Worlds without constellations of their own get figures drawn from the brightest stars
(`geoforge_cosmos::constellations`), which `constellation auto` saves to edit.
Open: visiting another galaxy's stars, superclusters' filaments, more than one system in a
world.

## 8i. Temperate defaults, habitable systems, moons and rings — done (November 2026)
OpenSpec `habitable-moons-rings`. A climate estimate before generating
(`worlds::climate_estimate`: the band energy balance the climate layer is built on, 2.5 K
colder for the land; within ~2 K of the generated means) grades every solid world temperate,
habitable or not (`worlds::habitability`). The default star rule asks first for a Sun-like
star whose best world is temperate (all eight test seeds find one within 2 s); the Earth-like
score weighs habitability. The viewer circles systems with open water in the star view
(green temperate, blue other; `/api/habitable`, ~1 s for 800 systems, cached), marks bodies in
the system view and gives each body's estimate; `geoforge show` prints it. Moons come from one
generator (`moons`): 3–6 regular moons for gas giants and 3–5 for ice giants, shared over 2.5
decades, small inner moons, Triton-like captures for a third of ice giants, irregular moons as
a population, captured specks beside asteroid belts. Every giant has rings, a quarter of cold
ones bright ice (`Population::optical_depth`: Saturn's 0.4), drawn banded in the viewer. The
star chart shows the zodiacal light and gegenschein, scaled by the system's belts.

## 8b. A system's local environment — next

A system's site is only its distance from the galactic centre. Young systems should often
sit in or near where they were born: an open cluster, an OB association, the nebula that
made them, or a superbubble like the Sun's Local Bubble and Gould Belt. Draw a local
environment from the system's age and the galaxy's star formation (cluster membership
falling with age; nearby OB associations by the local young-star density); put its
clusters, nebulae (emission, reflection, dark clouds) and supernova remnants into the sky
chart and catalogue; let it raise the sky's young-star count (the Gould Belt gap) and the
supernova risk.

## 8f. Report — done

`geoforge report` writes one page (`index.html` and `img/`) to take a world's pulse:
maps grouped as rock, air and sea, water and life, resources, hazards, each with its key
(`render::legend` for products, the view's classes or ramp for derived views); numbers
against Earth's; three places picked by Köppen class (or named) with hazards and weather;
the weather map; the sky chart with a key to magnitudes, spectral colours and planets.
Open: a verdict section written by hand stays outside it.

`geoforge compare REAL MODEL…` (November 2026) writes the Earth A/B page the same way: for
each product the real world holds as authored observations, the observed map and each
model's with one key, the numbers between them and the values at named probes
(`pipeline/src/compare.rs`). It replaced `ab_page.py`.

## 8e. Mountains, gyres and the weak climates — todo

From the Thumi pulse (November 2026). The Earth checks in brackets are what each item should
reproduce on the ETOPO-forced run or a generated Earth-like world.

**Mountains** (none needs a full history; the mock-history case is the Appalachians):

- **Paleo-orogens around cratons** — done. Belts run through the accreted ground between
  the craton and the coast, parallel to both, at a share of the way out, on the side with
  the most room toward a passive coast (an active margin's ground is its foreland); a
  large continent's second on another side, older and closer in. They wander, swell and
  pinch, vary in height along their length, thin out at their ends, and never cross a
  craton or each other (the X across Australia is gone).
- **Sutures between cratons** — done: a continent over 20 million km² is welded from two
  cratons along its long axis, over 40 million from three, and an old belt 150–300 km wide
  runs across it midway between each neighbouring pair, wandering ±250 km, 300–1,900 Ma
  old (the Urals, the Trans-Hudson). Open: a suture is always a belt; some are sag basins
  instead (the Amazon between the Guiana and Brazilian shields).
- **Arcs vary along strike** — partly done: cordilleras ±20% and island arcs more in
  ~900 km stretches, with volcano-free gaps and larger islands on the strongest stretches
  (peaks to 2.4 km). Open: tie the variation to slab age and dip rather than noise.
- **Dynamic topography** — done: continents far behind any active margin's arc or orogen
  and on slow plates stand high in patchy swells (Africa, Brazil); those just behind an arc
  (the Amazon) or on fast plates (Australia, India) stay low. On Earth's real plates land
  height now correlates 0.57 with ETOPO1 (0.27 before) and 39% of land is above 500 m
  (ETOPO1 39%, before 13%).
- **Flat slabs lift land far inland** — done. Where an oceanic hotspot track or plateau
  reaches an ocean–continent trench, block ranges of old basement rise 450–1,000 km
  inland along ~700 km of the margin (the Sierras Pampeanas over the Juan Fernández
  Ridge). A margin over 2,500 km long with no ridge arriving has an even chance of a past
  flat slab's ranges 700–1,500 km inland, 50–80 Ma old: the mock history for the Rockies.
  The ranges are ~180 km blocks up to ~3.4 km on a 1.2 km rise; their rocks are basement
  cores with sedimentary basins between, and they count as orogens for ore. Open: no
  volcanic gap over the flat stretch of the arc.
- **Plume domes** — done for generated worlds: every continental hotspot and flood basalt
  under 45 Ma sits on a 1.5 km dome ~650 km across; continental tracks run to 800 km.
  Earth's real plumes (Afar, Kenya) are not inputs to the real-plates run.
- **Island arcs.** Peaks are capped at 1.2 km; mature arcs reach 2–2.5 km and build bigger
  islands as their crust thickens with age. Back-arc rifting can split a continental sliver
  off an Andean margin. [Aleutians vs Marianas; Japan as a rifted sliver.]
- **Great escarpments** on rifted margins (optional). [Drakensberg, Western Ghats, Serra
  do Mar.]

**Oceans** — mostly done. The interior flow is right to be mostly zonal; the gyres now
close:
- an equatorward jet runs along shores where the winds upwell (the eastern legs);
- upwelling counts every season with equatorward wind (a summer of it chills California
  all year) and needs an ocean-wide sea;
- the trades upwell a cold tongue along the equator in a wide basin's east, which leaves
  the warm pool warm;
- western boundary currents leave the coast midway between where the westerlies begin
  and the gyre front, then climb across the basin, and poleward of that a cold current
  comes down the coast to meet them (the Labrador Current, Oyashio, Malvinas);
- a current carries heat across the gradient only ~1,500 km before it mixes;
- the current map's colours stretch the slow end, so the eastern legs show.

On the ETOPO-forced Earth, the SST departure boxes (observed in brackets) are:

| Box | Before | Now | Observed |
|---|---|---|---|
| California | −1.9 | −3.8 | −4.2 |
| Benguela | −2.5 | −3.8 | −4.4 |
| Humboldt | −3.9 | −4.4 | −4.9 |
| Canary | −3.0 | −3.9 | −3.0 |

On the A/B chain, the western Pacific warm pool goes from +0.1 to +1.1 (+1.3) and the
pattern correlation from 0.33 to 0.46.

Mass conservation (November 2026). The gyres are drawn from ψ on the cells' corners with
ψ = 0 on every corner touching land, so they are exactly divergence-free and cross no coast
(`atmosphere/flow.rs`). The jets are projected off the coasts, keeping only the overturning's
sinking; the wind drift keeps its Ekman divergence. Projecting the whole sum instead turned
every current that met a shore into a radial outflow from it (the Kuroshio stopped:
−0.02 m/s) and cost 8 climate cases. The Kuroshio now runs north at 0.26 m/s and the
Labrador Sea is cold (−3 °C); on the A/B chain the SST error is 2.47 °C (2.55 before) and
the anomaly correlation 0.49 (0.46).

**Winds round ranges** (November 2026): the surface winds are a layer of stable air up to
3.5 km above sea level, thinner over higher ground, projected (`Keep::FlatLayer`) so it goes
round the Rockies, the Alps and Tibet and keeps the flat-ground pattern's rising and
sinking, scaled by the layer that is left (unscaled, a plateau's thin layer had to draw in a
whole column's air and the air off it reached 25 m/s). Direction agreement with ERA-Interim
in January: 80% of the ocean (78%) and 75% of land below 1.5 km (72%). The climate suite
lost Mongolia's rain, the Pampas and the Great Plains (gaps with causes); the base winds
lack the East Asian summer monsoon, so the turned westerlies dry China's interior.

**Seasons apart** (November 2026): the warmest and coldest months are computed apart, each
with its own season's continentality, a shorter summer reach for sea air, frozen seas barely
maritime, and a winter high over high-latitude interiors. Judged on observed winds with
`geoforge compare`, the land temperature error fell from 3.84 to 3.67 °C and the range
north of 50°N is within 1 °C; on the model's own winds the suite lost 5 cases to its
January winds over Europe and its July winds over the Sahara.

Open:
- The NE Atlantic — done: a basin that opens onto a polar sea through straits wider than
  ~100 km (not the Bering Strait), without crossing a circumpolar band, draws an
  overturning branch poleward along its eastern side, which keeps giving up heat past
  where the zonal gradient flattens at freezing. The NE Atlantic box went from +2.4 to
  +4.0 °C (observed +4.4). The climate lets that open water stand above a zonal mean set
  by frozen land, so the Norwegian Sea stays open and the Bering Sea freezes.
- The Brazil Current runs warm (+4 to +5 °C off southern Brazil vs +1.4).
- Western Australia and the Timor Sea upwell where the Indonesian Throughflow keeps them
  warm (the Leeuwin Current).
- Somalia and Oman upwell only in the summer monsoon, which an annual anomaly cannot
  show.

**Weak climates**: the layer-by-layer Earth A/B (`scripts/real-earth/ab.sh`, November 2026)
fed the climate real winds and sea temperatures, which splits the misses:

- *Circulation's*: Phoenix 806 → 237 mm (real 200) and Tehran 406 → 224 mm (230) on real
  winds, so the dry Southwest US and Iran fail upstream (winter storms reach too far
  south). Mumbai 1,003 → 1,707 mm (2,200): the weak monsoon is mostly the winds too. The
  cold eastern currents ran at half strength (California −1.8 K vs −4.2, Benguela −2.4 vs
  −4.4); see Oceans above.
- *Climate's own*: Lima stays wet on real winds and sea (1,025 → 766 mm, real 15); Yakutsk's
  winter stays mild (range 45 vs 58 K) and its rain near 0 (240); Antarctic seasonal ranges
  are about double; typical land rain is half of GPCP's while some coasts (London 1,949 mm
  vs 620) are far too wet.
- *Routing*: on GPCP rain 11 of 18 rivers land within a factor of two (7 on model rain);
  the Congo stays blocked at its rapids (1,025 m³/s vs 41,000).
- *Mountains from real plates*: land above 1 km 7% vs 25%; the Andes and Tibet come out,
  the Alps do not (52 m mean vs 972), Iran too high (2,229 vs 1,351); intraplate ranges
  absent as expected.
- The Köppen classifier on the real climate agrees with the Köppen–Geiger map on 87% of
  land (the model climate 58–60%), so it is not the classifier.

Findings since (November 2026):
- *Thumi's weak seasons are physics, not a bug*: its 168-day year and heavy air column
  (1.49 bar at 0.8 g, 1.8× Earth's air mass per m²) give the column 4× Earth's thermal
  inertia against the seasons; with an 18.9° obliquity its swing is ~0.27× Earth's. Its
  summers stay too cool for boreal forest, and its ice matches its land poleward of 60°
  (21% of land).
- *The Congo is a resolution artifact*: on 10′ ETOPO at a third of a degree, the gorge
  through the Crystal Mountains (a few km wide) averages to a 430–470 m barrier, so the
  basin spills over the lower divide into Gabon. Breaching sinks by least cost might help,
  but a gorge this narrow is below the data's resolution.
- *Frozen seas* (tried, not kept): letting a sea frozen in winter send continental air
  over land fixed Utqiagvik (mean −10.0 °C, range 26; observed −11.0, 31). But it chilled
  Iceland further through the model's then too-cold East Greenland water, halving
  Reykjavik's rain, and the colder Arctic coast rains more (336 mm vs 140). Worth a retry
  now that the NE Atlantic is warm.
- *Siberian winters*: the missing winter cold (Yakutsk −24 °C in the coldest month vs −38,
  Astana −7 vs −15, Winnipeg −7 vs −16) does not follow continentality: Moscow is more
  continental than Astana in the model and already too cold. Nor does it follow the
  season's sinking air: the model's winter mid-latitudes converge everywhere (1.5–2.2 at
  every station), so it cannot tell a winter high from a storm track. It wants a winter
  pressure field (the Siberian and Canadian highs) to key inversions on, with Europe's
  maritime reach fixed alongside.
- *Frozen seas, retried after the Atlantic fix*: counted for the whole year they still cost
  Reykjavik's rain; counted for winter alone they chill every Arctic coast 5–9 °C too far
  (Utqiagvik −16 °C, Summit −37). The climate's symmetric seasons (one continentality for
  both halves of the year) are the limit: a seasonal energy balance with sea ice is the
  real fix.

The first guesses, before the A/B:

| Miss | Likely cause | Fix to try |
|---|---|---|
| West-coast deserts wet (Lima 1,021 mm vs 15, Phoenix 1,011 vs 200) | weak cold current and upwelling; 1° grid smears mountain rain onto the coast | stronger eastern legs; rain on the windward cell only; test with real SST |
| Siberian winters mild | sea air carried too far inland in winter; no snow-cover cooling or inversion | longer winter decay of maritime reach; snow cooling (Yakutsk −39 °C in January) |
| Indian monsoon July jet a third of real | no cross-equatorial jet along the East African highlands; Tibet not an elevated heat source | elevated heating; let the cross-equatorial flow form a jet against the highlands |
| West Pacific too cool | no warm pile-up from the trades | equatorial SST tilt tied to trade strength |
| Rivers blocked by gorges | gorges below grid scale; basins filled, not breached | breach the lowest pass through a barrier |

**Sky**: neighbourhood stars are drawn from the galaxy's statistics, not its systems; tie
the brightest ones to systems that can be generated and visited.

## 9. Deep-time tectonics — later

Artifexian's procedure (Logs 10–21), as NotebookLM gives it: a supercontinent of ~25% of
the surface, round old cratons. Rifts start at triple junctions between cratons, never
through them: two arms spread and the third fails (an aulacogen). Child plates split off
by Euler poles; ridges make ocean floor in 10-Myr steps; subduction builds arcs and
back-arc basins. Arcs, microcontinents and slivers ride until they hit a margin, where
they are welded on as terranes. Failed rifts become lowlands that steer big rivers inland.
Accreted arcs and sutures host the VMS, Sedex, tin and gold. Old cratonic sutures hold the
greenstone belts whose tropical weathering makes gold laterites.


Artifexian's core method: run the plates forward from a supercontinent (GPlates, 500–850
Myr), so sutures, old orogens, failed rifts, accreted terranes and ocean ages come from a
history. We infer the history from a present-day plate snapshot. A forward kinematic
history would place failed rifts that steer big rivers, terranes and their ophiolites,
and give resources their ages (coal basins, banded iron). Large; worth doing.

## Still open (from earlier work)

Done (OpenSpec `systems-and-worlds`): ice shells (`geoforge_surface::shell`: ridged plains,
chaos, bands and ridges on a thin shell, dark and grooved terrain on a thick one); locked
oceans carry the day side's heat (an eastward oceanic superrotation, Zeng & Yang 2021, and
the sea-surface anomaly on the substellar baseline); thermal-tide rotation of locked worlds
with 0.5–5 bar (Leconte et al. 2015; Valente et al. 2024); brown dwarfs (Burrows & Liebert
1993; a fifth of systems, Kirkpatrick et al. 2024) and planets of wide companions; giants'
polar cyclone rings and polygonal jets; fast spinners' ice cliff (the balance lets frost in
no faster than radiation and transport answer) and rain lines at their fronts.

Done (OpenSpec `other-worlds-2`): thermal tides of thick air (the day lengthens as
(p/5 bar)² past 5 bar, which gives Venus its 117-day day; such worlds turn backward half the
time, Correia & Laskar 2001); binary evolution of a system's inner pair
(`geoforge_cosmos::binary`: Roche-lobe overflow by Eggleton, Hurley et al.'s critical mass
ratios, stable transfer with a rejuvenated companion, common envelopes with αλ = 0.5,
mergers); the other stars' light on each star's planets; a gyre front under each band of
westerlies on fast spinners; Titan's methane cycle (`MethaneCycle`: the cold trap's rain and
evaporation, no sea, a 0.85 m inventory filled coldest first: ~1% lakes poleward of 55°)
over a drainage core shared with water (`hydrology::drain`); a `dunes` view for any world
with rivers (Fryberger's RDP/DP for the form).

The substellar point of a locked world is 0°E by the IAU's convention for synchronous bodies
(the mean sub-parent point); its continents fall anywhere relative to it, as the tectonics
seed places them.

Open:
- Land is assumed even for a locked world's weathering (a continent under the star would
  hold less CO₂): the air is decided before the surface exists.
- Giants: storms do not drift by date; colour by class only.
- Titan: the climate gives no methane rain, so the dunes view bounds its dune belt by lakes
  and rivers alone (±50° against ±30°), and its pole-to-equator contrast is 12 K, not 3; its
  seas lie equally north and south, and its rain ignores the seasons.
- Binary evolution: second transfers, donors of 8 M☉ or more, contact binaries; the sky's
  stars live as if alone (a pair's members are worked out only where its system is built).

## Notes on the Artifexian summaries

Most of the summaries' numbers hold up; these do not:

- Frost line 4.85·√L AU puts the Sun's near 5 AU; the usual value is ~2.7 (Hayashi 1981),
  which we use.
- "Cowler's gravity relief law" (h ∝ 1/g^1.6–2): no source found; height scales about as
  1/g. Summary 1 contradicts itself (its formula gives 7.3 km at 1.2 g; its table says ~5).
- The circumpolar current flows eastward, not westward.
- 3 light years is ~190,000 AU, not 10,000.
- "Faolore / Mississippi Valley-type" native copper is a mix-up: fahlore is a copper
  sulfosalt ore, Mississippi Valley-type deposits are lead–zinc. Neither summary mentions
  porphyry copper.
- A giant planet as an impact shield is overstated: Horner & Jones (2008–2012) found
  Jupiter lowers long-period comet impacts but can raise short-period comet and asteroid
  impacts.
- Ocean crust "collapsing after 300–400 Myr": Earth's oldest ocean floor is ~180–200 Myr
  (a Mediterranean relic aside).
- Andean ranges "80–200 km wide": the central Andes with the Altiplano are ~700 km.
- Debated rather than settled: super-Earths tending to stagnant lids; moonless worlds'
  chaotic tilt; plant and sky colours ("Panoptes-5" could not be traced).
- Köppen's C/D boundary: the series uses 0 °C (Logs 36–40), as we do by default.

NotebookLM's answers to our follow-up questions are folded into items 2, 3, 4, 6 and 9
above. Their own checks:

- The 4.85·√L frost line is from *GURPS Space*, a tabletop game's rule.
- "Cowler" is Kawaler: S. D. Kawaler (1988)'s magnetic-braking law for how stars spin
  down. NotebookLM's "Benton R. Cowler" does not exist. So "Cowler's gravity relief law" in
  summary 1 is probably garbled too.
- Panoptes-5 is named as a web tool for sky and plant colours by spectral class and
  pressure; not traced.
- The 56/33/8/3% multiplicity split is for Sun-like stars only.
- Artifexian's tree line, auroral rings and cyclone cutoffs match measured values well.
