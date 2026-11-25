# nebulae-and-galaxy-centres

## Why

Galaxies had stars, dust and molecular clouds but no nebulae to find, no black holes at
their centres and no globular clusters; the night sky showed none of them.

## What Changes

- **The just dead among the stars**: two bins after the dark one hold white dwarfs born
  within 25 kyr and neutron stars and black holes born within 100 kyr, from strata of the
  field tables and from birth events' members at their exact age, in 2 kpc cubes so the few
  there are can be found far off. They are real stars with addresses.
- **Nebulae** (`geoforge_cosmos::nebulae`), found where their makers are: H II regions
  round birth events under 10 Myr with ≥ 10⁴⁷ ionising photons/s (Strömgren radius in
  100 cm⁻³ grown at 10 km/s), reflection nebulae round B stars still in gas, molecular
  clouds over 10⁴ M☉ as dark nebulae, planetary nebulae round the just-dead white dwarfs and
  supernova remnants round the just-dead neutron stars and black holes (free expansion,
  then Sedov–Taylor). Each is named, sized, aged and given an integrated brightness.
- **Nebulae are not an address level.** Each is identified by what makes it (a birth
  event, a cloud, a star) and leads to the star that lights or left it; see the design.
- **Galaxy centres and globular clusters** (`geoforge_cosmos::nuclei`): a black hole from a
  classical bulge's mass or, for pseudobulges and dwarfs, the galaxy's stars, not every dwarf
  holding one; a nucleus feeding at an Eddington ratio drawn by type, named quasar, Seyfert,
  LINER or quiet, some with jets (Sgr A* for the real Milky Way); globular clusters weighing
  2.9×10⁻⁵ of the halo, metal-poor in the halo and metal-rich nearer the bulge.
- **The sky**: `Sky::nebulae` and `Sky::globulars` list what a world sees through the dust;
  `geoforge sky` lists them and draws them on the chart (`SkyChart`, `SkyPatch`), toned by
  surface brightness against the Milky Way's.
- **The viewer**: nebulae glow by kind round the home star in the galaxy view, with its
  globular clusters and centre; the star view lists the nearest round a system and draws
  those reaching into it; each tells about itself and leads to its star; search finds them
  by name and by the fictional materials they hold; active galaxies are ringed in the
  universe view.

## Impact

- `geoforge-cosmos`: `sectors::bins` (`NASCENT_WHITE_DWARFS`, `NASCENT_COLLAPSED`),
  `sectors::tables`, `sectors::groups`, `nebulae` (new), `nuclei` (new), `resolve`
  (`GalaxyAt::{halo_msun, nucleus, globulars}`), `universe::members`, `exotics`.
- `geoforge-almanac`: `Sky::{nebulae, globulars}`.
- `geoforge-store`: `render_sky(width, &SkyChart)`.
- `geoforge-viewer`: `nebulae` (new: `/api/nebulae`), galaxy and universe scenes, search.
- A star view's queries now visit the just-dead bins' cubes; lit and dark stars are drawn as
  before.
