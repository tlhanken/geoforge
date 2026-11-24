# Design: universe-and-figures

## Halos (`universe::halos`)
- Mass function, a fit to the Tinker (2008) abundance at z = 0 (h = 0.7):
  N(>M) = 8.9×10⁻⁴ (M/1.43×10¹²)^−0.9 exp(−(M/2.54×10¹⁴)^0.676) Mpc⁻³; checked at
  1.43×10¹², 10¹⁴ and 10¹⁵ M☉ to 5–15%.
- Thirteen half-decade bins from 10⁹ M☉. Each bin has its own cube size, large enough that
  a cube holds two or more halos on average; a cube's halos are a pure function of the
  universe's seed, the bin and the cube.
- Clustering: a lognormal field exp(σg − σ²/2), σ = 0.9, g from value-noise octaves of 32,
  16, 8 and 4 Mpc; a halo of mass M is kept with weight exp(b(M)σg − (b(M)σ)²/2) against
  the ceiling, b(M) = 0.7 + 0.35 (M/10¹²)^0.35 (heavier halos cluster more). Checked: the
  expected number within 15%, Milky Way halos' excess scatter in 10 Mpc cells 0.5–1.0.

## Galaxies in a halo (`universe::members`)
- Central: stellar mass from Moster et al. (2013) at z = 0 (M₁ = 10^11.59, N = 0.0351,
  β = 1.376, γ = 0.608), 0.15 dex scatter; a Milky Way halo of 1.3×10¹² M☉ holds
  3–6×10¹⁰ M☉ of stars.
- Satellites: subhalos above m, 0.05 (m/M)^−0.9 per host, each kept with chance
  exp(−(m/M)/0.15) (heavy ones sink by dynamical friction), placed in an NFW profile
  (c = 5) within R200; stellar mass from the subhalo's infall mass. Inside a tenth of R200
  the central's tides tear them apart: one survives at r with chance 1 − exp(−(r/0.1 R200)²)
  (Garrison-Kimmel et al. 2017). Checked: a Milky Way halo has 7–20 satellites above
  10⁵ M☉ (median distance 50–130 kpc, under 3% within 20 kpc) and an LMC in 8–25% of halos.
- Types: early-type share by stellar mass (half at 6×10¹⁰ M☉) raised in groups and for
  satellites; dwarfs below 10⁹ M☉ irregular or spheroidal (satellites under 10⁷ M☉ always
  spheroidal). Sizes from Shen et al. (2003); V-band mass-to-light by type.
- Identity `h{bin}:{i},{j},{k}:{n}.{member}` (member 0 the central), seeding everything
  about the galaxy.

## A world's galaxy (`universe`)
- `UniverseSite`: `Real` (the Solar System) or `Seeded`; absent, the universe is seeded
  from the galaxy's seed. The home galaxy is the central nearest the origin of the same
  class (disk, elliptical, dwarf) within 0.25 dex of the model's stellar mass, looked for
  in shells to 128 Mpc.
- The home galaxy's frame: its axis is the galactic north pole; the azimuth is drawn from
  its seed.

## Galaxies in the sky (`universe::sky`, `almanac::sky`)
- Each galaxy within 40 Mpc that could reach the magnitude limit is placed in the
  observer's frame: distance, magnitude, apparent axis ratio √(cos²i + t² sin²i), the
  major axis on the sky.
- Surface brightness: an exponential disk, Σ₀ = L/(2π h² q) L☉/pc² (independent of
  distance), dimmed by the home galaxy's dust out along that ray. The glow is in units of
  the galactic pole's brightness, so galaxies are added as Σ/Σ_pole.
- Listed size: 2 × 3.5 scale lengths (about the 25 mag/arcsec² isophote).
- How it is seen goes by its central surface brightness μ₀ = 26.4 − 2.5 log Σ₀ as well as
  its magnitude: the naked eye needs V ≤ 6.5 and μ₀ ≤ 22.5, binoculars μ₀ ≤ 23, a telescope
  μ₀ ≤ 24, else long exposures. A dwarf spheroidal of magnitude 4 at μ₀ 25 is invisible, as
  the Sagittarius dwarf is.
- The Sun: the LMC, SMC, Andromeda, Triangulum, Centaurus A, Bode's, Sculptor and the
  Southern Pinwheel at their catalogue positions, magnitudes and sizes. Checked: the LMC at
  l 280.5°, b −32.9°; Andromeda at l 121.2°, b −21.6°; on Earth's chart the LMC is more
  than twice as bright as the sky 40° away and Andromeda is a smudge; μ₀ 21.4 for the LMC
  (observed about 21), 20.2 for Andromeda's disk.

## Automatic constellations (`constellations`)
- Stars to magnitude 4; each joined to its 8 nearest; a minimum spanning forest with no
  line longer than 15°; trees over 10 stars split at their longest lines; a first-magnitude
  star left alone joined to its two nearest stars to magnitude 5. Each figure is named from
  the sky's seed.
- Checked: a Sun-like sky gets 40–120 figures of 2–10 stars, no line over 15°, every star
  brighter than 1.5 in one.
