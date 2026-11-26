# size-knees-and-m-dwarfs

## Why

Two misses found while adding asteroids and nebulae. One power law per belt gave Sol's Kuiper
belt ~1,000 bodies over 1000 km (the real belt has ~10) and gas giants' Trojan swarms hundreds
over 100 km (Jupiter's: ~20). Torres's bolometric correction runs away below ~3200 K, so M
dwarfs — most stars — came out 1–4 magnitudes too faint, thinning every neighbourhood's
count by brightness and every sky's faint end.

## What Changes

- `PowerLaw` gains an optional knee (`Knee { km, exponent }`, `with_knee`) past which the
  differential slope steepens, with `count_above` and `diameter_at_rank`; belt members take
  their sizes from it. Sol's Kuiper belt and generated outer belts steepen to D⁻⁵ past 100 km,
  Trojans to D⁻⁵·⁵ past 50 km (Fraser et al. 2014).
- Below 3900 K the visual bolometric correction is interpolated from nearby M dwarfs'
  luminosities and magnitudes (Lacaille 8760 to TRAPPIST-1), joined to the fit above and the
  L dwarfs below.

## Impact

- `geoforge-core`: `PowerLaw` (a serde-defaulted field; stored systems read as before).
- `geoforge-cosmos`: `sol`, `galactic` (outer belts, Trojans), `context`, `sky`.
- M dwarfs move to brighter bins, so the stars drawn round any point change; outer belts'
  largest members shrink.
