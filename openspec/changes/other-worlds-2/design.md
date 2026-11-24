# Design: other-worlds-2

## Thermal tides of thick air

Leconte et al. (2015) and Valente et al. (2024) worked out the thermal tide for ~1 bar. A
thicker atmosphere's radiative time grows with its mass, so the bulge lags less and the
torque balance falls at a longer day. One scaling past 5 bar, (p/5 bar)², is fixed by Venus:
92 bar at 0.72 AU gives a 117-day solar day from the 1-bar curve's ~8 h there. Which of the
two stable spins a thick-aired world falls into depends on things the model does not carry
(early spin, tilt, core–mantle friction: Correia & Laskar 2001, 2003), so it is a draw per
body from the system seed; backward spin is stored as the tilt past 90° (Venus 177°).

## Binary evolution

Scope: the inner pair of a system, donors under 8 M☉. The donor is the heavier star at birth.
The first time (sampled through the main sequence and finely through the giant branch) its
radius passes its Roche lobe at periastron:

| Donor | q = Md/Ma | Outcome |
|---|---|---|
| main sequence | > 3 | merger |
| main sequence | ≤ 3 | contact binary, not modelled: the stars live on apart |
| subgiant | ≤ 4 | stable transfer |
| giant | ≤ (1.67 − 0.3 + 2(Mc/M)⁵)/2.13 | stable transfer |
| otherwise | | common envelope, merger if the companion then overfills its own lobe |

Core mass from Refsdal & Weigert's core-mass–luminosity relation, capped at the IFMR's white
dwarf. Stable transfer: half the envelope accreted, the rest lost with the accretor's specific
orbital momentum, integrated in 200 steps; the accretor rejuvenated (its lived share of the
main sequence, scaled by mass, Hurley et al. 2002). Common envelope: αλ = 0.5 energy budget.
Orbits are circular afterwards, then widened by later winds. Planets keep forming around the
pair as it was born (`born_close`); the tree shows it as it is now (a post-common-envelope
pair is a barycentre).

## A companion's light

`host()` sums the light of the stars it is made of; the others shine at their separation
from it: the companion's own orbit when it orbits the host, the host's when the host orbits
it, otherwise the outer of the two orbits about their shared centre. Averaged over the orbit,
L/(4πa²√(1−e²)). The planets' own distances are ignored (stability keeps them inside a third
of the separation). The light goes into the year's mean flux, not into the star's colour or
the diurnal cycle.

## Gyres per jet

`currents` takes the jet count (`jets_for`). The zonal stress is smoothed over a share of a
jet's band (10°/jets); with more than one jet, every local maximum of westerly stress from 15°
to 85° is a front, its separation half-way from the westerlies' start; the jet, its rise and
the eastern-shore run scale with the band; the overturning stays with the most poleward
front. One jet keeps Earth's rule exactly.

## Methane cycle and dunes

`hydrology::drain(ctx, inputs, &Water, params)` is the old hydrology run, its inputs made
explicit: rain and evaporative demand (mm/yr of the liquid), an optional temperature for
frozen and glacial lakes, whether there is a sea, and an optional liquid inventory. Water's
generator passes Budyko's inputs as before (Earth unchanged). The methane cycle passes
R·C^(−z) and D·C^z by the ground's temperature against the world's mean (R 20, D 60 mm/yr,
C 3), no sea (the world drains to its lowest point), no freezing, and an inventory of 0.85 m
over the world, lakes filled where the demand is least. Calibrated on Titan: ~1% lakes, all
poleward of 55°.

The dunes view is generic (not keyed to the body): sand where rain is under 150 mm/yr, slope
under 1 in 50, and no lake or river within 150 km, on a world that has rivers at all; form by
RDP/DP of the January and July winds, each season's drift ∝ speed³. On Titan the climate
gives no methane rain, so its dune belt is bounded only by its liquids (to ~±50°, against the
observed ±30°): recorded as a gap.
