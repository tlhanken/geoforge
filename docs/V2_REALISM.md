# Geoforge v2 realism check

## Earth-like worlds

Six worlds at 2048×1024 (about 18 km cells at the equator), generated with
`geoforge-v2 new … && geoforge-v2 generate … && geoforge-v2 stats …`:

| World | Body |
|---|---|
| `earth42` | `--archetype earth`, seed 42 |
| `s2837` | Tevo IV: 5799 km, 0.89× Earth's starlight, 0.72 bar, day 30.5 h, obliquity 20.6° |
| `s42` | Kibedo IV: 5592 km, 0.88×, 0.84 bar, day 28.7 h, obliquity 18.4° |
| `s7` | Lisano IV: 6308 km, 0.93×, 0.74 bar, day 22.8 h, obliquity 23.3° |
| `s1234` | Vivi IV: 6278 km, 1.02×, 1.19 bar, day 36.8 h, obliquity 28.9° |
| `s99` | Shulor III: 6344 km, 0.91×, 1.46 bar, day 20.5 h, obliquity 29.1° |

A world takes about 37 s with the `iter` build (tectonics 7 s, geology 7 s, elevation with
erosion 17 s, circulation 1 s, climate 5 s, hydrology 2 s, biomes <0.1 s).

| Quantity | earth42 | s2837 | s42 | s7 | s1234 | s99 | Earth |
|---|---|---|---|---|---|---|---|
| Largest plate, of surface | 0.25 | 0.22 | 0.18 | 0.24 | 0.25 | 0.24 | 0.20 |
| Seven largest plates | 0.81 | 0.81 | 0.82 | 0.85 | 0.83 | 0.81 | 0.81 |
| Land fraction | 0.31 | 0.23 | 0.28 | 0.29 | 0.30 | 0.31 | 0.29 |
| Mean land elevation, m | 804 | 1482 | 736 | 1137 | 774 | 940 | 840 |
| Land above 1000 m | 0.29 | 0.76 | 0.15 | 0.65 | 0.15 | 0.23 | 0.27 |
| Highest point, m | 8171 | 10682 | 9549 | 9629 | 10472 | 10501 | 8849 |
| Mean ocean depth, m | 3773 | 3465 | 3811 | 3357 | 3684 | 4010 | 3680 |
| Surface rock: Cenozoic / Precambrian | 0.24 / 0.38 | 0.28 / 0.34 | 0.23 / 0.38 | 0.25 / 0.36 | 0.27 / 0.39 | 0.29 / 0.38 | 0.28 / 0.32 |
| Global mean temperature, °C | 14.1 | 8.4 | 9.5 | 10.1 | 16.2 | 12.6 | 14.5 |
| Equator (5°S–5°N), °C | 25.8 | 21.5 | 22.0 | 21.3 | 24.8 | 21.1 | 26.5 |
| 40–50°N, °C | 8.4 | 6.1 | 4.7 | 6.1 | 12.0 | 9.4 | 10 |
| 55–65°N, °C | −2.6 | −9.3 | −11.3 | −9.5 | 5.1 | 3.1 | 0.5 |
| Annual range, land 40–60°, °C | 32.5 | 37.9 | 26.8 | 35.5 | 38.0 | 20.1 | ~28 |
| Annual range, ocean 40–60°, °C | 5.9 | 5.7 | 5.1 | 7.1 | 6.7 | 3.6 | ~7 |
| Global precipitation, mm/yr | 1024 | 927 | 897 | 929 | 1089 | 980 | 1000 |
| Land precipitation, mm/yr | 837 | 740 | 592 | 974 | 1071 | 761 | 780 |
| Land drier than 250 mm/yr | 0.28 | 0.46 | 0.41 | 0.25 | 0.14 | 0.24 | ~0.2 |
| Precipitation 5°S–10°N, mm/yr | 2084 | 1865 | 2008 | 1690 | 1821 | 1321 | ~2000 |
| Precipitation 40–55°N, mm/yr | 646 | 620 | 505 | 740 | 611 | 971 | ~1000 |
| River discharge to the sea, km³/yr | 46 991 | 38 343 | 22 806 | 67 279 | 55 643 | 41 664 | 37 000 |
| Rivers over 10 000 m³/s | 24 | 17 | 8 | 38 | 24 | 23 | 8 |
| Lakes over 10 000 km² | 7 | 4 | 9 | 2 | 8 | 18 | 16 |
| Endorheic land | 0.08 | 0.00 | 0.06 | 0.00 | 0.01 | 0.15 | 0.18 |
| Mediterranean biome, of land | 0.011 | 0.007 | 0.010 | 0.015 | 0.008 | 0.014 | 0.02 |
| Hot desert, of land | 0.06 | 0.01 | 0.05 | 0.04 | 0.07 | 0.04 | 0.14 |
| Temperate grassland, of land | 0.27 | 0.22 | 0.24 | 0.24 | 0.28 | 0.36 | 0.08 |

The generated worlds are not Earth: most get 0.88–0.93× Earth's starlight and differ in
size, air pressure, day length and tilt. The carbonate–silicate thermostat (set in the solar
system step) cancels about 75% of the starlight difference, so they run 4–6 °C cooler, and
more at high latitudes where ice brightens the ground. Thick air damps the seasons (s99,
1.46 bar: 20 °C land range at 40–60°); thin air widens them.

Drainage networks: on the Earth-like world at 1024×512, the exponent of the drainage-area
distribution (P(A ≥ a) ∝ a^−β) went from 0.66 without erosion to 0.54 with it; Earth's
networks give about 0.43–0.45.

## Other worlds

Presets (`--archetype`, real sizes, orbits and air) at 2048×1024, seed 42. Elevation is
relative to the mean radius.

| Quantity | Mars | Venus | Moon | Mercury |
|---|---|---|---|---|
| Generators | lid ×3, circulation, climate | lid ×3, circulation, climate | lid ×3 | lid ×3 |
| Lowest / highest, km | −10.4 / +19.2 | −3.4 / +11.1 | −8.8 / +5.0 | −8.1 / +5.0 |
| Real, km | −8.2 / +21.2 | about −3 / +11 | 19.9 km range | −5.4 / +4.5 |
| Surface within 1 km of the mean | 0.30 | 0.82 (real: over 0.80) | 0.52 | 0.30 |
| Lowland (thin crust) share | 0.35 | 0.80 | 0.18 (maria ≈ 0.17) | 0.27 |
| Global mean temperature, °C | −66 (real ≈ −63) | 465 (real 464) | — | — |
| Equator / north pole, °C | −48 / −111 | 465 / 457 | — | — |
| Seasonal range, 40–60°, °C | 71 | 2 | — | — |

The climate model's greenhouse is fitted to Earth, Mars, Venus and Titan presets
(`cargo run -p geoforge-core --example radiation_check`, and
`-p geoforge-surface --example zonal_check` for zonal profiles).

Generated systems also give Venus-like inner worlds (runaway greenhouse, 35–360 bar CO₂),
cold super-Earths with thin CO₂, ice-shell ocean moons around giants, frozen cratered moons
and belt members; `geoforge-v2 show` lists each body with its regime, air and water.

## What matches

- Plate sizes (largest ~0.2–0.25, seven largest ~0.81–0.85); continents inside plates with
  passive margins.
- Ocean depth by crust age, land fraction and hypsometry for Earth-like water; surface rock
  ages close to Earth's mix.
- Rivers carve dendritic valleys; mountain belts stay high.
- Western boundary currents run warm and eastern margins cold with upwelling; zonal
  temperatures within ~2 °C of Earth from the equator to 50° on the Earth preset; land/ocean
  seasonal ranges; global rainfall and the ITCZ; summer-dry west coasts in the subtropics
  (Mediterranean ~1–2% of land).
- Mars, Venus, the Moon and Mercury: relief ranges, Venus's narrow hypsometry, Mars's
  dichotomy and volcanic rises, mean temperatures of Mars and Venus.

## Known gaps

1. **Mid-latitude interiors are dry and grassland is over-represented** (0.22–0.36 of land
   vs 0.08); hot deserts under-represented (0.01–0.07 vs 0.14) as subtropics run ~2 °C cool.
2. **Too many large rivers** (8–38 over 10 000 m³/s vs 8): continents here are bigger and
   less broken up than Earth's, and runoff is somewhat high.
3. **Drainage networks branch less than Earth's** (0.54 vs ~0.44); longer erosion flattens
   the land without improving it, so the remaining gap is likely resolution.
4. **Some worlds are very high** (s2837: 76% of land above 1000 m): a broad plateau term can
   lift a whole large continent. `plateau_m` controls it.
5. **Ocean currents are a diagnostic gyre model**, not a simulation.
6. **Closed basins depend on terrain** (endorheic land 0–15% vs Earth's 18%).
7. **The gray greenhouse warms the tropics more than the poles** for added CO₂ (real CO₂
   forcing is flatter), so dim, CO₂-rich worlds have cold high latitudes.
8. **Mercury's basins come out deep** (−8 km vs −5.4 km), and the Moon's highest highlands
   low (+5 vs +10.8 km): basin depth does not yet depend on the lithosphere's heat.
9. **No generator yet for ice-shell surfaces** (Europa, Titan) or giants; their layers report
   why they are skipped. Titan's methane cycle is out of scope.
