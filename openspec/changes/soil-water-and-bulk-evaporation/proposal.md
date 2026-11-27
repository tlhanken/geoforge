# soil-water-and-bulk-evaporation

## Why

On Earth's own map the rain was the weakest layer of the A/B page: a pattern correlation of
0.59 on the model's winds and 0.62 on observed ones, the same either way, so the fault was in
how rain is made. Checked against NCEP's column water as well as GPCP's rain, the moisture
model showed:
- land gave back only the rain that had just fallen on it, so a dry summer interior stayed
  dry (West Siberia's July air held 6 mm of water, observed 20; Yakutsk got 34 mm a year,
  observed 355);
- the sea relaxed passing air to 85% of saturation within hours, so narrow seas (the Red
  Sea, the Gulf) loaded every air mass that crossed them;
- the storm tracks mixed vapour at a tenth of the observed eddy diffusivity;
- the column's capacity followed the ground's full seasonal swing, so warm air over frozen
  winter land condensed as soon as it arrived;
- nothing carried vapour away from under sinking air over land, and with land returning all
  of a desert's rain the Sahara doubled its own.

## What Changes

- **Soil water** (`moisture::budyko_evaporation`, `moisture::soil_evaporation`): land
  evaporates the Budyko share (Fu, ω 2.6) of its year's rain and evaporative demand, spread
  over the seasons by their demand. The climate layer runs both seasons from a wet start,
  then twice more on the rain they made, each from its vapour.
- **Bulk sea evaporation**: E = C_E |U| (q_sat − q) on the column (C_E 1.3 × 10⁻³, vapour
  scale height 2.2 km, gusts 4 m/s).
- **Eddy diffusivity in m²/s**: 3 × 10⁴ + 8 × 10⁵ × storminess.
- **The column's capacity** reckoned on 60% of the ground's seasonal departure.
- **Orographic gain** at most twice the surroundings' rate (was four times).
- **Subsidence export over land**: 7 × 10⁻⁷ s⁻¹ per unit of ln(1/ascent multiplier).
- `transport::diffuse_rows` runs its rows in parallel, with sub-steps of up to 0.5.
- The climate generator's version goes 13 → 14.

## Impact

- `geoforge-surface`: `climate` (`moisture`, `mod`), `atmosphere::transport`. The climate
  layer takes about twice as long (6.9 s against 3.8 s for a 1024-wide Earth on two cores).
- Every wet world's rain changes, and with it its temperature (desert warming follows
  aridity), rivers, lakes, biomes and Köppen view.
- Earth A/B: rain correlation 0.59 → 0.68 (model winds), 0.62 → 0.71 (observed winds);
  Köppen by group 61% → 68%. The Earth climate suite: 186 → 190 of 253 cases pass; new
  known gaps for the Great Plains, land rain overall, Phoenix's warmth, the Cerrado and South
  Atlantic storms (`known_earth_climate`).
- Generated Earth-like worlds: land drier than 250 mm 0.34 → 0.24 of land (Earth 0.20);
  land rain 975 → 1090 mm (Earth 780); rivers over 10 000 m³/s 21 → 28 (16).
