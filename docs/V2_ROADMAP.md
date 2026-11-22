# Geoforge v2 roadmap

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

Done as derived views (`geoforge_surface::resources`; see `V2_PIPELINE.md`): prospectivity
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

Open:

- The ladder steps distance by 1.45–2.1 (period ratios 1.75–3.0); Kepler's multi-planet
  systems peak at period ratios ~1.5–2. Tighten small planets toward that.
- Scattering could also make hot or eccentric Jupiters by high-eccentricity migration.
- Obliquity stability: a big moon steadies the tilt; without one, a wider tilt draw
  (debated: Laskar 1993 vs Lissauer et al. 2012).
- Trojans at the giants' L4/L5 points.

Also done: the `positions` command and export (every star, planet and moon's place at a
date, by Kepler orbits), and close binaries now share one orbit about their barycentre
(each star had drawn its own tilt and phase).

## 6. Sky, calendar and almanac — calendar done; the rest next

Done: the focus's solar day, year and synodic months (`show`), its local date for any day
and its moons' phases (`positions`). Next, for any date:

- the star and moons in the sky: angular sizes, phases, eclipse type (total, annular,
  none), **eclipses** (new or full moon near a node, with the nodes precessing);
- sunrise, sunset and day length by place; seasons; tides;
- a **star map**: the stellar neighbourhood drawn from the same galaxy population
  (~0.1 stars per cubic parsec near the Sun), apparent magnitudes and colours,
  constellations, the galaxy's band; a check that no supernova progenitor lies within
  ~10 pc;
- **weather on that day** at a place: a seeded weather generator (wet/dry spells,
  temperature anomalies, storms) consistent with the climate's means and seasons;
- natural phenomena on that date: cyclone season, monsoon onset, auroras, meteor showers.

Artifexian (Log 6): tides scale as the moon's mass over its distance cubed, plus the
star's (spring and neap tides). The calendar comes from the year (Kepler), the day, and the
synodic month, 1/P_syn = 1/P_orbit − 1/P_year. His neighbourhood (Log 3) is a 12-light-year
sphere at 0.004 stars per cubic light year (34 stars, 20 systems), with Raghavan et al.
(2010)'s multiplicity for Sun-like stars (56% single, 33% binary, 8% triple, 3% more). We
already draw multiplicity by star mass, which matters because red dwarfs are mostly single.

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

## 8. Editing near the poles — next

Edit views are equirectangular PNGs with world files, which QGIS reads and Blender maps
onto a UV sphere as is (Artifexian's workflow: paint on the sphere, save the
equirectangular image, import it back). To finish that:

- accept an edited image at any 2:1 size (resample on import);
- optional polar cap views (azimuthal, poleward of ~60°) for painting in flat apps,
  merged back on import where they were changed;
- a short guide to the Blender round trip.

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

- Ice-shell worlds (Europa) have no generator.
- Tidally locked worlds: ocean currents do not carry heat to the night side; the
  substellar point is always 0°E; land is assumed even for weathering (a continent under
  the star would hold less CO₂).
- Asynchronous rotation: planets with air around 0.5–0.7 M☉ stars may escape locking
  (Leconte et al. 2015), though perhaps not around smaller stars (Salazar & Wordsworth
  2024).
- No brown dwarfs; companion stars get no planets; no binary evolution.
- Fast rotators (5–7 cells) show thin row-by-row streaks in ocean rain; ocean gyres follow
  one jet.
- Giants: storms do not drift; no polar hexagon; colour by class only.

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
