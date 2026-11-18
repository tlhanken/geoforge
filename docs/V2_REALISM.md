# Geoforge v2 realism check

Six worlds at 2048×1024 (about 18 km cells at the equator), generated with
`geoforge-v2 new … && geoforge-v2 generate … && geoforge-v2 stats …`:

| World | Body |
|---|---|
| `earth42` | `--earth-like` reference body, seed 42 |
| `s2837` | Tevo IV: 5799 km, 0.89× Earth's starlight, 0.72 bar, day 30.5 h, obliquity 20.6° |
| `s42` | Kibedo IV: 5592 km, 0.88×, 0.84 bar, day 28.7 h, obliquity 18.4° |
| `s7` | Lisano IV: 6308 km, 0.93×, 0.74 bar, day 22.8 h, obliquity 23.3° |
| `s1234` | Vivi IV: 6278 km, 1.02×, 1.19 bar, day 36.8 h, obliquity 28.9° |
| `s99` | Shulor III: 6344 km, 0.91×, 1.46 bar, day 20.5 h, obliquity 29.1° |

Generation takes 20–30 s per world on one machine (tectonics 5 s, geology 4 s, elevation 8 s,
climate 3 s, hydrology 1 s, biomes <0.1 s). Editing elevation and regenerating takes ~12 s.

## Numbers

| Quantity | earth42 | s2837 | s42 | s7 | s1234 | s99 | Earth |
|---|---|---|---|---|---|---|---|
| Land fraction | 0.35 | 0.32 | 0.26 | 0.37 | 0.34 | 0.30 | 0.29 |
| Mean land elevation, m | 924 | 708 | 657 | 1292 | 870 | 954 | 840 |
| Land above 1000 m | 0.27 | 0.19 | 0.12 | 0.62 | 0.22 | 0.21 | 0.27 |
| Mean ocean depth, m | 4004 | 3918 | 3713 | 3746 | 3874 | 3945 | 3680 |
| Global mean temperature, °C | 14.0 | 10.6 | 10.7 | 10.9 | 16.1 | 13.1 | 14.5 |
| Equator (5°S–5°N), °C | 25.6 | 21.6 | 21.2 | 20.7 | 25.7 | 20.1 | 26.5 |
| 40–50°N, °C | 8.9 | 5.9 | 8.5 | 5.6 | 12.0 | 9.4 | 10 |
| 55–65°N, °C | −4.8 | −4.3 | −4.6 | −2.3 | 3.9 | 2.6 | 0.5 |
| Annual range, land 40–60°, °C | 32.7 | 24.6 | 23.0 | 30.6 | 42.9 | 27.2 | ~28 |
| Annual range, ocean 40–60°, °C | 5.9 | 4.9 | 3.9 | 5.4 | 7.1 | 6.5 | ~7 |
| Global precipitation, mm/yr | 997 | 896 | 943 | 899 | 1065 | 1004 | 1000 |
| Land precipitation, mm/yr | 807 | 596 | 1021 | 666 | 968 | 968 | 780 |
| Ocean precipitation, mm/yr | 1097 | 1039 | 915 | 1036 | 1113 | 1019 | 1100 |
| Land drier than 250 mm/yr | 0.25 | 0.44 | 0.36 | 0.46 | 0.31 | 0.22 | ~0.2 |
| Land wetter than 2000 mm/yr | 0.07 | 0.05 | 0.13 | 0.09 | 0.11 | 0.08 | ~0.11 |
| Precipitation 5°S–10°N, mm/yr | 1835 | 1739 | 2000 | 1377 | 1802 | 1469 | ~2000 |
| Precipitation 40–55°N, mm/yr | 729 | 531 | 716 | 823 | 542 | 1120 | ~1000 |
| Precipitation 70–90°N, mm/yr | 106 | 201 | 195 | 209 | 414 | 458 | ~250 |
| River discharge to the sea, km³/yr | 51 900 | 30 200 | 53 300 | 51 100 | 64 300 | 66 100 | 37 000 |
| Largest river, m³/s | 71 100 | 33 100 | 43 600 | 255 200 | 54 200 | 60 700 | 209 000 |
| Lakes over 10 000 km² | 26 | 3 | 13 | 23 | 17 | 35 | 16 |
| Lake area / land | 0.007 | 0.002 | 0.013 | 0.007 | 0.006 | 0.008 | 0.018 |
| Endorheic land | 0.03 | 0.14 | 0.14 | 0.22 | 0.07 | 0.03 | 0.18 |
| Largest plate, of surface | 0.12 | 0.12 | 0.09 | 0.15 | 0.12 | 0.13 | 0.20 |

Biomes, share of land, earth42 vs Earth: ice sheet 0.03 / 0.10 (no polar continent in this
world), tundra 0.06 / 0.06, boreal forest 0.15 / 0.11, temperate forest 0.16 / 0.08, temperate
grassland 0.22 / 0.08, Mediterranean 0.003 / 0.02, hot desert 0.04 / 0.14, cold desert
0.06 / 0.05, savanna 0.13 / 0.15, tropical rainforest 0.03 / 0.07, tropical seasonal forest
0.04 / 0.05.

The generated worlds are not Earth: most get 0.88–0.93× Earth's starlight and differ in
size, air pressure, day length and tilt, so cooler means and different zonal values are
expected. The weathering thermostat (default 0.75) keeps them habitable rather than snowballs.

## Edit check

`earth42` with a 2.6 km coastal range painted into the elevation view (2437 cells, imported
as a patch). Within ~12 s, climate, rivers and biomes follow:

| Band | Before, mm/yr | After, mm/yr |
|---|---|---|
| Windward coast | 2216 | 2738 |
| On the range | 1792 | 5389 |
| Leeward, first ~500 km | 1140 | 525 |

Boreal forest in the new rain shadow becomes grassland; the peaks become alpine.

## What matches

- Ocean depth from crust age, land fraction, hypsometry (for Earth-like water amounts).
- Zonal temperatures within ~2 °C of Earth from the equator to 50°; continental vs oceanic
  seasonal range; global precipitation, land/ocean split and the ITCZ.
- Rain shadows, wet windward coasts, dry continental interiors and subtropical deserts.
- Rivers, lakes and drainage basins of plausible size and count; endorheic basins appear in
  dry interiors; lakes persist only where rifting, depth or glaciation keep basins closed.

## Known gaps

1. **Plates are too uniform in size** (largest 9–15% of the surface vs Earth's 20%; seven
   largest 55–65% vs 84%). v1's electrostatic seeding is unweighted. Fix: size-weighted
   seeding or plate merging; it changes plates for every v1 seed, so it is a deliberate
   decision.
2. **High latitudes run cold** (55–65°N about 4 °C below Earth) and subtropics about 3 °C:
   there is no ocean heat transport by currents (no Gulf Stream), and hot deserts are
   under-represented (0.04 vs 0.14 of land) partly because of it.
3. **Mid-latitude continental interiors are dry** (40–55°N about 730 vs 1000 mm on the
   Earth-like world) and **polar precipitation is low** on cold worlds. Moisture transport is
   a two-season, one-layer model; no frontal cyclones or subtropical-high gyres.
4. **Too many large rivers** (rivers over 10 000 m³/s: 18–39 vs 8; runoff about 20% above
   Earth's per land area). Continents here are bigger and less broken up than Earth's.
5. **Closed basins depend on terrain**: endorheic land is 3–22% (Earth 18%) and lake area
   0.2–1.5% (Earth 1.8%, much of it lakes too small for 18 km cells).
6. **Province shapes**: wide collision orogens can have straight edges where their boundary
   runs straight, giving rectangular plateaus.
7. **Mediterranean biome is rare**: the climate products are annual, so summer drought is
   approximated from latitude, range and rainfall.
8. Elevation has no fluvial erosion (valleys are not carved); sediment fill stands in for it.
