# Design: systems-and-worlds

## Planets (`planets`, `dynamics`)
- Inner chain: 1 + Poisson(2) planets (2.6 for red dwarfs) on consecutive rungs from the
  innermost period (1.5–12 d), period ratios exp(ln 1.8 + 0.3·N(0,1)); the chain's own mass
  10^U(−1.2, 1.1)·√M M⊕ (and ∝ M for brown dwarfs), each planet ± 0.25 dex; heavier than
  3 M⊕ a sub-Neptune.
- Giants: 1 + Poisson(0.7), each after the first on the next rung at period ratio 1.5–2,
  born with e up to 0.2; 2% of giant systems migrate their first giant in the disk to the
  innermost rung.
- Settling: Hill spacing 8 (push apart), scatter under 5 or when Petrovich (2015) fails.
  Safronov number Θ = (m/M)(a/R) of the heavier planet: under 1 merge (e 0–0.1), else eject
  the lighter, survivor e = max(e, Rayleigh(0.3)). Giant–giant scatterings make a hot Jupiter
  7% of the time: periastron 1–3 Roche radii (2.16 R (M/m)^(1/3)), circularised at twice it;
  planets inside the giant's old orbit are lost.
- Checked: median adjacent period ratio 2.1 (Kepler peak 1.6, long tail); 1.08 small planets
  within 200 d per Sun-like star, 27% sub-Neptunes, neighbours' mass spread 0.36 dex against
  1.0 for strangers; giants' median e 0.17, 23% over 0.4 (Kipping 2013: 0.17, 18%); hot
  Jupiters 0.76% of Sun-like stars, 7% with an inner neighbour against 33% for warm ones;
  1.8 small planets per red dwarf within 200 d.

## Spin (`spin`)
- Big moon: G·m/a³ ≥ G·M/a_p³ (the Moon is 2.2 × the Sun on Earth). Without one, today's
  tilt is U(birth − R, birth + R) within 0–90°, R ~ U(10°, 45°).
- Thermal tides (0.5–5 bar): a_c as above; solar day 30 d·(a/a_c)^−4.7 (one day at 0.62 AU
  for 0.82 M☉), prograde. Thicker air (Venus) is left locked.
- Locked rotation: the year, or 2/3 of it when e > 0.1. The habitable recipe now locks its
  close siblings too; the focus's own draws are untouched (separate streams).

## Stars (`stars`, `population`)
- Brown dwarfs under 0.08 M☉: L = 4×10⁻⁵ (t/Gyr)^−1.3 (M/0.05)^2.64, T = 1550 K
  (t/Gyr)^−0.32 (M/0.05)^0.83 (from 10 Myr), types M/L/T/Y; disks lit by their glow at
  10 Myr. 20% of drawn systems (Kirkpatrick et al. 2024's 4:1); masses ∝ M^−0.6 below
  0.05, M^−0.25 above; companions of half to equal mass. Sector stars stay stars.
- Wide companions: their own disk inside the S-type limit with the primary as perturber,
  formed from a separate draw stream; planets named "B I", "B II". Light from the other star
  is ignored.

## Trojans
- Giants ≥ 10 M⊕ with e < 0.1 and no heavier giant within 15 mutual Hill radii (the Solar
  System: Jupiter and Neptune yes, Saturn and Uranus no); a ± a·√(8μ/3); 10^5–10^6 bodies
  (10^5.5–10^7 beyond 10 AU), largest 100–250 km.

## Ice shells (`geoforge-surface::shell`)
- Style by thickness: thin under 25 km, thick over 60 km, blended between. Thin: renewed 30%
  (chaos patches + bands), plains ~60 Ma, renewed ~25 Ma. Thick: renewed (grooved) 65%,
  dark terrain ~4 Ga, grooved ~2 Ga.
- Provinces: ridged plains / chaos or cratered ice / grooved terrain; basins; bands
  (5% of a thin shell, 300–1500 km long, 10–60 km wide); ridge lineae (~36 on a thin
  shell, 800–3000 km).
- Relief in the shell's own metres: regional ±400 m (×2 thick), chaos +100 m with ±250 m
  blocks, bands +100 m, ridges 200 m, grooves ±400 m; craters shallower beyond 8 km
  ((8/D)^0.6); palimpsests at 20% of a basin's depth.

## Locked oceans and fast spinners
- `Baseline`: per atmosphere cell, the sea-level temperature before currents and its band
  (row, or row of substellar latitude). The SST anomaly takes −(u − ū)·∂T/∂x −
  (v − v̄)·∂T/∂y with each reach capped as before and removes each band's mean: identical
  for rotating worlds.
- Ocean jet on locked worlds: 1 m/s eastward, e-folding 25°·√(P/30 d) (8–30°), at 20–100%
  on the night side by period (Zeng & Yang 2021).
- Frost width in the balance: max(10 K, 1.25·s·0.16/(B(frost) + C)); Earth unchanged
  (s ≈ 220 W/m² at its ice edge gives 8.8 K).
- Multi-jet winds averaged over ±10% of a band pair (weights 1-2-3-2-1).

## Giants' poles
- Gas giants: half the poles a ring of 5–8 cyclones ~4300 km across at 6–8° from the pole
  plus one at it; else a polar vortex in a jet at 10–14°, round (40%) or a 5–7-gon. Ice
  giants: a round polar vortex.

## Left open
- Storms that drift with their latitude's wind (the bands product is one moment).
- Ocean gyres of a fast spinner still follow one jet.
- Venus-like thick air's thermal tides; binary evolution; light from a companion star on
  the other star's planets; Titan's dunes and methane lakes.
