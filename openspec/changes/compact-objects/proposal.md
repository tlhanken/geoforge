# compact-objects

## Why

The galaxies' stars had no brown dwarfs (the initial mass function started at 0.08 M☉) and
dropped their dark remnants (a neutron star or black hole had no brightness bin), so a
quarter of the objects near the Sun were missing and no pulsar or black hole could be found.
A system with no solid body could not be opened in the viewer.

## What Changes

- **A dark bin** after the 30 lit ones holds brown dwarfs and dark remnants. Brown dwarfs are
  a stratum of their own, one to every four stars born, by their own mass function (slopes
  0.6 and 0.25, breaks at 0.01, 0.05 and 0.08 M☉); their mass is left out of the
  populations' mass, so the lit bins' counts are unchanged. A brown dwarf's bolometric
  correction below 2000 K follows the L/T/Y sequence.
- **Remnants told** (`geoforge_cosmos::remnants`): a neutron star's spin from a birth
  period and field drawn per star (Faucher-Giguère & Kaspi 2006), dipole spin-down, the death
  line, magnetars while young and strong; millisecond pulsars recycled by a companion; a
  black hole's horizon; feeding by Roche-lobe overflow or a hot giant's wind.
- **Barren systems open**: `resolve::system_or_barren` gives a system with no solid body
  for viewing.
- **The viewer** colours stars by kind (L, T, Y dwarfs; white dwarfs; neutron stars; black
  holes ringed), counts them in the star view, draws a black hole dark with a halo and a
  pulsar's beams, and lists each star's notes; `geoforge show` prints them.

## Impact

- `geoforge-cosmos`: `sectors::bins` (`DARK`), `sectors::tables` (brown dwarf stratum),
  `population` (brown dwarf mass function), `sky` (bolometric correction), `resolve`,
  `remnants` (new).
- `geoforge-analysis`: `appearance::star_colour`.
- `geoforge-viewer`, `geoforge-cli`.
- Lit stars are drawn as before; the dark bin adds stars, so a neighbourhood's star list
  grows (universe 7's: 80,188 within 40 pc, 15,628 of them brown dwarfs).
