# other-worlds-2

## Why

Five gaps in worlds unlike Earth (TODO item 14) left them simpler than the physics already in
place allows: thick-aired worlds near their star stayed locked, close pairs of stars lived as
if alone, a companion star's light never reached the other star's planets, a fast spinner's
oceans had one pair of gyres whatever its jets, and Titan had no methane on its surface.

## What Changes

- **Thermal tides of thick air** (`geoforge-cosmos::spin`): air past 5 bar holds its heat
  longer, so the day the thermal tide settles to lengthens as (p/5 bar)², which gives Venus
  its 117-day solar day; such a world turns backward half the time (Correia & Laskar), its
  tilt then past 90°. `thermal_tide_solar_day_h`, `backward_rotation_h`.
- **Binary evolution** (new `geoforge-cosmos::binary`): a pair whose heavier star fills its
  Roche lobe (Eggleton) gives its envelope to the other (stable transfer: a white dwarf and a
  rejuvenated companion, the orbit widened), throws it off round the other (common envelope,
  αλ = 0.5: a white dwarf and a red dwarf days apart), or merges; Hurley et al.'s critical
  mass ratios decide. `galactic::build` builds the inner pair with it. Donors of 8 M☉ or
  more, second transfers and the sky's stars are left as before.
- **A companion's light** (`geoforge-cosmos::context`): a planet's flux adds the system's
  other stars' light at their separation from its host, averaged over that orbit.
- **Gyres on fast spinners** (`geoforge-surface::circulation`): one gyre front under each
  band of westerlies; the separated jets and subpolar coastal currents repeat per band, and
  the overturning stays beyond the most poleward. Earth (one jet) is unchanged.
- **Titan's methane** (new `MethaneCycle` generator for the hydrology layer, and a `dunes`
  view): methane's cold trap wets the coldest ground and dries the warmest; rivers, lakes
  and watersheds come from the same drainage as water's (refactored into a shared core),
  with no sea and a liquid inventory filled coldest first (Titan: ~1% lakes, poleward of
  55°). The `dunes` view maps sand seas on any world with rivers: dry, flat ground away from
  liquid, transverse, linear or star dunes by Fryberger's RDP/DP of the seasons' winds.
- Not done here: storms that drift by date on giant planets (left for later).

## Capabilities

### Modified Capabilities
- `systems-and-worlds`: close pairs evolve together; companions light the other star's
  planets; thick air frees locked worlds.
- `surface-layers`: a methane hydrology; stacked gyres; the dunes view.

## Impact

- `geoforge-cosmos`: `binary` (new), `spin`, `context`, `galactic`, `stars` (doc).
- `geoforge-surface`: `methane` (new), `hydrology` (`drain`, `Water`), `circulation`
  (`currents` takes the jet count; version 13), `preconditions::surface_methane`.
- `geoforge-analysis`: `dunes` (new). `geoforge-pipeline`: `Implementation::MethaneCycle`,
  `View::Dunes` in `climate_maps`, the report's hydrology group, stats count a dry world's
  river ends at lakes.
- Generated systems with close pairs, thick-aired locked worlds or companions change;
  stored worlds keep their system and regenerate only the circulation layer of fast
  spinners (version bump).
