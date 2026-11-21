# Geoforge v2 roadmap

What comes next, in order, and why. Sources: our own gaps (`V2_REALISM.md`), and a
comparison with Artifexian's *Worldbuilder's Log* series (Logs 0–48), from two NotebookLM
summaries of it (fact-checked below). Each item says what it adds, what it reads, what it
writes, and which Logs it follows.

Status: **next** (agreed, not started), **later** (agreed, not scheduled), **open** (a gap in
what exists).

## Principles

- **Bottom-up causation**, as now: star → system → body → tectonics → geology → elevation
  → circulation → climate → erosion → hydrology → biomes → (new) resources, hazards.
- **Stored products are things you can edit.** A classification or summary that is a pure
  function of stored products (Köppen climate, sea-level pressure) is a **derived view**:
  computed on export or render, never stored, never stale, never edited.
- **Checked against Earth** where Earth is the answer, as the known-case suites do now.

## 1. Climate at a glance: Köppen–Geiger view — next

Today the biome layer writes 16 classes (ocean, lake, ice sheet, tundra, boreal forest,
temperate rainforest, temperate forest, temperate grassland, Mediterranean, cold and hot
desert, savanna, tropical seasonal forest, tropical rainforest, alpine, wetland): a
Whittaker-style vegetation map that borrows Köppen's thresholds where they are physical
(the 10 °C warm-month tree line, the aridity threshold, dry summers).

Add the Köppen–Geiger climate classes (~30: Af/Am/Aw, BWh/BWk/BSh/BSk, Cs/Cw/Cf with
a/b/c, Ds/Dw/Df with a/b/c/d, ET/EF) as a derived view from `temperature`,
`temperature_range` (warmest and coldest month as mean ± range/2), `precipitation` and
`precipitation_warm_share`.

- Choose the C/D boundary on purpose: Köppen's −3 °C or the common modern 0 °C (Peel et
  al. 2007). Make it a parameter; default 0 °C.
- Use Köppen's own aridity threshold, which already scales with temperature (no separate
  "cold places are not deserts" correction).
- Gap to note: there are no monthly series, only half-year rain shares, so "driest month"
  rules (Af vs Am, Cs vs Cf) are approximated from the warm/cold split.
- Check against Earth: compare the ETOPO-forced Earth run with the Beck et al. (2018)
  Köppen map, class by class (area shares and per-cell agreement).

## 2. Derived views — next

The mechanism for item 1 and the pressure map: a view names the products it reads and
renders (and exports) from them; `export` and `render` offer views alongside products.

**Sea-level pressure (Jan, Jul), hPa**, from the circulation pattern and land heating:
the equatorial trough (~1008–1010), subtropical highs (~1020–1025, strongest over cool
ocean in summer), subpolar lows (~990–1000), polar highs; winter highs over cold
continents (Siberia ~1035) and summer heat lows (~1000–1005). Rendered with isobars and
H/L marks. Check: surface winds should cross isobars toward low pressure at 20–40° (as
they do in the stored winds), and against ERA5's seasonal pressure maps.

## 3. Resources and soils — next

A new layer after biomes. Writes a **deposit catalog** (each deposit: kind, area, grade,
size, depth/accessibility, age) and fields (soil class and fertility, timber, fishery
productivity). Placement follows the geology, climate and history already stored, with
probabilities rather than certainties.

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

Check: Earth's own deposits by type and setting (USGS mineral resource data) against the
ETOPO-forced Earth run: do porphyry belts sit on the Andes and Rockies, tin on the
Southeast Asian and Cornish granites, banded iron on the Pilbara and Minas Gerais cratons,
salt in the Permian and Zechstein basins?

## 4. Hazards — next

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
- **Auroras**: ovals at 65–75° magnetic latitude; needs a magnetic field model (dynamo
  from core size, rotation and age) — which would also matter for air loss.

Check: Earth's earthquake catalog (USGS) by region and magnitude, IBTrACS cyclone tracks,
the Smithsonian volcano list, against the ETOPO-forced Earth run.

## 5. Star system fixes — next

- **Orbit spacing by Hill radii.** A planet's Hill radius r_H = a·(m/3M)^(1/3) is where its
  gravity beats the star's tide (Earth 0.01 AU, Jupiter 0.36). Neighbours are measured in
  mutual Hill radii, ((m₁+m₂)/3M)^(1/3)·(a₁+a₂)/2. Two planets are safe beyond 2√3 ≈ 3.5
  (Gladman 1993); systems of several need ~8–10 to last billions of years (Chambers et
  al. 1996); Jupiter–Saturn sit at ~8. Small planets can pack tightly (two Earths need
  only a 1.13 distance ratio for 10), giants cannot (two Jupiters need ~2.05, a period
  ratio ~2.9, for 8). Our ladder steps distance by 1.45–2.1 whatever the masses: two
  Jupiters at 1.45 are only 4.3 mutual Hill radii apart. Enforce ≥ 8 for every
  neighbouring pair (push outward, or drop).
- **The ladder's spacing.** The comments call 1.45–2.1 "period ratios" but it is the
  distance ratio (period ratio 1.75–3.0). Kepler's multi-planet systems peak at period
  ratios ~1.5–2; tighten small planets toward that.
- **Belts set by resonances.** Where an asteroid's orbital period is a simple fraction of
  a giant's (half: 2:1), the giant tugs it at the same point every few orbits, the kicks
  add up, and it is thrown out: Jupiter's 4:1, 3:1 and 2:1 resonances (2.06, 2.50 and
  3.28 AU) bound and gap our asteroid belt (the Kirkwood gaps). Outside Neptune the 3:2
  and 2:1 (39.4 and 47.7 AU) hold the Kuiper belt's plutinos and outer edge. Distance from
  a period ratio: a = a_giant·(P/P_giant)^(2/3). Use these instead of today's fixed
  fractions (inner belt to 0.72 of the first giant; outer belt 1.3–1.9 of the last planet).
  Trojans at the giants' L4/L5 points would follow naturally.
- **Breathability** in `show`: O₂ partial pressure 0.16–0.5 atm, CO₂ under 5000 ppm,
  N₂ under ~3 atm (narcosis), gravity. Today a 0.7-bar Earth-like world gives 0.15 atm
  of O₂, and a locked Earth 6200 ppm of CO₂.
- **Moons**: Roche-limit check; tidal range from moon and star; obliquity stability (a
  big moon steadies the tilt; without one, a wider tilt draw — debated: Laskar 1993 vs
  Lissauer et al. 2012).

## 6. Sky, calendar and almanac — next

A world clock (epoch, year, day, months from the moons, a calendar) and, for any date:

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

## 7. Oceans — next

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
- Köppen's C/D boundary: summary 1 uses −3 °C, summary 2 uses 0 °C; both conventions exist.

## Questions out to NotebookLM

1. Logs 45–48: what inputs and numbers place each ore, coal, salt and peat.
2. Logs 42–43: placement rules, with numbers, for cyclones, tornado alleys, fog, dust and
   auroras.
3. Logs 10–21: the supercontinent breakup steps, and how sutures, failed rifts and
   terranes feed rivers and ores later.
4. Logs 29–32: how the seasonal pressure maps are drawn before the winds.
5. Logs 36–40: which Köppen convention (−3 or 0 °C), which aridity formula, and the
   a/b/c/d subtype rules.
6. Logs 6 and 3: the tide and calendar method; how the stellar neighbourhood is made.
7. Where the 4.85·√L frost line, "Cowler" and "Panoptes-5" come from.
