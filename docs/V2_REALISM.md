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

A world takes about 40 s (tectonics 7 s, geology 7 s, elevation with erosion 17 s,
circulation 1 s, climate 5 s, hydrology 2 s, biomes <0.1 s).

## Numbers

| Quantity | earth42 | s2837 | s42 | s7 | s1234 | s99 | Earth |
|---|---|---|---|---|---|---|---|
| Largest plate, of surface | 0.25 | 0.22 | 0.18 | 0.24 | 0.25 | 0.24 | 0.20 |
| Seven largest plates | 0.81 | 0.81 | 0.82 | 0.85 | 0.83 | 0.81 | 0.81 |
| Land fraction | 0.31 | 0.23 | 0.28 | 0.30 | 0.30 | 0.31 | 0.29 |
| Mean land elevation, m | 804 | 1482 | 736 | 1137 | 773 | 940 | 840 |
| Land above 1000 m | 0.29 | 0.76 | 0.15 | 0.65 | 0.15 | 0.23 | 0.27 |
| Highest point, m | 8171 | 10682 | 9549 | 9629 | 10472 | 10501 | 8849 |
| Mean ocean depth, m | 3773 | 3465 | 3811 | 3357 | 3684 | 4010 | 3680 |
| Surface rock: Cenozoic / Precambrian | 0.24 / 0.38 | 0.28 / 0.34 | 0.23 / 0.38 | 0.25 / 0.36 | 0.27 / 0.39 | 0.29 / 0.38 | 0.28 / 0.32 |
| Global mean temperature, °C | 14.4 | 9.5 | 10.4 | 11.2 | 16.4 | 13.0 | 14.5 |
| Equator (5°S–5°N), °C | 26.1 | 21.3 | 22.1 | 21.2 | 25.3 | 22.1 | 26.5 |
| 40–50°N, °C | 8.7 | 7.5 | 5.3 | 7.2 | 12.1 | 9.6 | 10 |
| 55–65°N, °C | −2.2 | −3.9 | −8.1 | −5.9 | 4.9 | 2.8 | 0.5 |
| Annual range, land 40–60°, °C | 32.4 | 29.5 | 22.6 | 26.8 | 40.6 | 34.5 | ~28 |
| Annual range, ocean 40–60°, °C | 5.9 | 4.4 | 4.3 | 5.4 | 7.2 | 6.3 | ~7 |
| Global precipitation, mm/yr | 1033 | 965 | 928 | 961 | 1097 | 994 | 1000 |
| Land precipitation, mm/yr | 847 | 765 | 617 | 990 | 1088 | 801 | 780 |
| Land drier than 250 mm/yr | 0.27 | 0.46 | 0.40 | 0.25 | 0.13 | 0.22 | ~0.2 |
| Precipitation 5°S–10°N, mm/yr | 2099 | 1882 | 2035 | 1706 | 1851 | 1344 | ~2000 |
| Precipitation 40–55°N, mm/yr | 654 | 694 | 543 | 785 | 603 | 962 | ~1000 |
| River discharge to the sea, km³/yr | 47 400 | 41 900 | 25 000 | 68 500 | 56 800 | 43 000 | 37 000 |
| Rivers over 10 000 m³/s | 24 | 21 | 10 | 39 | 25 | 25 | 8 |
| Lakes over 10 000 km² | 7 | 5 | 7 | 3 | 8 | 15 | 16 |
| Endorheic land | 0.08 | 0.00 | 0.06 | 0.00 | 0.02 | 0.14 | 0.18 |
| Mediterranean biome, of land | 0.011 | 0.009 | 0.012 | 0.018 | 0.007 | 0.008 | 0.02 |
| Hot desert, of land | 0.06 | 0.01 | 0.05 | 0.04 | 0.07 | 0.06 | 0.14 |
| Temperate grassland, of land | 0.27 | 0.19 | 0.24 | 0.23 | 0.28 | 0.37 | 0.08 |

The generated worlds are not Earth: most get 0.88–0.93× Earth's starlight and differ in
size, air pressure, day length and tilt, so cooler means and different zonal values are
expected.

Drainage networks: on the Earth-like world at 1024×512, the exponent of the drainage-area
distribution (P(A ≥ a) ∝ a^−β) went from 0.66 without erosion to 0.54 with it; Earth's
networks give about 0.43–0.45.

## What matches

- Plate sizes (largest ~0.2–0.25, seven largest ~0.81–0.85); continents inside plates with
  passive margins.
- Ocean depth by crust age, land fraction and hypsometry for Earth-like water; surface rock
  ages close to Earth's mix.
- Rivers carve dendritic valleys; mountain belts stay high.
- Western boundary currents run warm and eastern margins cold with upwelling; zonal
  temperatures within ~2 °C of Earth from the equator to 50°; land/ocean seasonal ranges;
  global rainfall and the ITCZ; summer-dry west coasts in the subtropics (Mediterranean
  ~1–2% of land).

## Known gaps

1. **Mid-latitude interiors are dry and grassland is over-represented** (0.19–0.37 of land
   vs 0.08); hot deserts under-represented (0.01–0.07 vs 0.14) as subtropics run ~2 °C cool.
2. **Too many large rivers** (10–39 over 10 000 m³/s vs 8): continents here are bigger and
   less broken up than Earth's, and runoff is somewhat high.
3. **Drainage networks branch less than Earth's** (0.54 vs ~0.44); longer erosion flattens
   the land without improving it, so the remaining gap is likely resolution.
4. **Some worlds are very high** (s2837: 76% of land above 1000 m): a broad plateau term can
   lift a whole large continent. `plateau_m` controls it.
5. **Ocean currents are a diagnostic gyre model**, not a simulation: zonal flow comes from wind
   drift; currents do not exchange between basins except via the circumpolar band.
6. **Closed basins depend on terrain** (endorheic land 0–14% vs Earth's 18%).
7. Elevation now takes ~17 s at 2048×1024 because erosion reroutes water eight times.
