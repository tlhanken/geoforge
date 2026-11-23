# systems-and-worlds

## Why

Systems were built from their star, but their planets did not yet look like the ones found:
small planets were spaced more widely than Kepler's compact systems, giants were never thrown
onto eccentric or hot orbits, every rocky world kept its birth tilt whether or not a moon
steadied it, giants had no Trojans, brown dwarfs and the planets of companion stars were
missing, every locked world was locked even with air enough for thermal tides, and some
kinds of world had no generator at all (ice shells like Europa's) or showed artefacts (a
locked world's ocean carried no heat; a fast spinner's ice edge jumped between two rows; a
giant's poles were plain).

## What Changes

- Inner small planets form as Kepler's compact chains (log-normal period ratios about 1.8,
  masses alike within a system); giants form packed and many scatter (Petrovich's criterion
  for eccentric pairs); a scattering ejects (giants far out) or merges (rocky planets close
  in), by the Safronov number; some giant scatterings end as hot Jupiters by tidal
  circularisation, which clear the inner system.
- A rocky world's tilt wanders unless a big moon steadies it; locked worlds with 0.5–5 bar
  of air are held off by their thermal tides beyond a_c ≈ (0.4 − 0.12 M)·(1 + 0.8e²) AU.
- Trojan swarms at giants with no heavier giant within 15 mutual Hill radii.
- Brown dwarfs (one system in five), cooling by Burrows & Liebert; wide companions get their
  own planets.
- Ice-shell generators (tectonics, geology, elevation) for thin Europa-like and thick
  Ganymede-like shells; erosion passes any solid surface through.
- A locked world's ocean superrotates and carries the day side's heat across the terminator;
  the sea-surface anomaly works on any baseline (latitude or substellar latitude).
- The energy balance lets frost in no faster than radiation and transport can answer, so a
  weak-transport world's ice edge does not jump between bands; the multi-jet bands spread
  their fronts over a storm track.
- Giants' poles: rings of cyclones or a polar vortex in a polygonal jet.

## Impact

- `geoforge-cosmos`: `planets`, `dynamics`, new `spin`, `galactic`, `system`, `population`,
  `stars`, `census`.
- `geoforge-core`: `StellarStage::BrownDwarf`; ice-shell province kinds; `ShellRenewal`;
  the energy balance's frost.
- `geoforge-surface`: new `shell` module; `circulation` (ocean jet, general baseline);
  `climate` (`Baseline`); `bands` (poles); `atmosphere::cells` (multi-jet fronts);
  `preconditions`; `erosion`.
- `geoforge-pipeline`: three ice-shell implementations in the registry.
- `geoforge-store`: the dry elevation tint fits the world's relief.
- Galaxy-drawn seeds draw different systems again (brown dwarfs, new planet draws); the
  habitable recipe keeps its planets, but rocky worlds without a big moon get a new tilt.
