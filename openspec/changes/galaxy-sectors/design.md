# Design — galaxy-sectors

## The zoom

    universe ─ superclusters ─ galaxies        (v2-architecture, later)
                                  │
                               galaxy ─ components (where stars are, how old, how metal-rich)
                                  │
                               sectors ─ nested cubes, one level per brightness
                                  │
                               star ─ a stable identity
                                  │
                               system ─ built from the star (`galactic::build`)
                                  │
                               bodies ─ the surface pipeline

Every level's children are a pure function of the parent's seed and the child's address,
so nothing is generated until it is looked at and nothing generated needs storing: a sky
query generates the cubes it reaches, a visit generates one system. What is stored is what
was asked to be kept (a world's focus system and its layers) and what was edited. The same
scheme, cubes sized by how far their objects are seen, is meant to place a universe's
galaxies later (bright galaxies in big cubes, dwarfs in small ones), so zooming from a
universe to a planet is one mechanism at every scale.

A standalone system is the same path with fewer levels shown: its world's seed makes its
galaxy (the Milky Way model unless the world names another) and draws its place in it, and
its sky comes from that galaxy's sectors around that place.

## Galaxy components

A galaxy is a list of components, each a shape, a share of the stars, a star formation
history and a metallicity relation:

| Component | Shape | Stars | Milky Way (to check against sources in 1.2) |
|---|---|---|---|
| nucleus | Plummer sphere around a black hole | old and young, metal-rich | nuclear cluster, ~10⁷ M☉ within a few pc |
| bulge and bar | triaxial (bar) plus Sérsic spheroid | old, metal-rich | box/peanut bulge, bar half-length ~5 kpc |
| thin disk | exponential in radius; vertical scale height growing with age | all ages, the reference metals and gradient | scale length 2.6 kpc, 50 → 300 pc with age; 0.09 stars pc⁻³ near the Sun |
| spiral arms | logarithmic spirals that raise the young thin disk's density and star formation | young | number of arms and pitch angle per morphology |
| thick disk | exponential, scale height ~0.9 kpc | ~10 Gyr, [Fe/H] ≈ −0.5 | 2.1% of the thin disk's density near the Sun |
| halo | power law in radius, flattened | ~12 Gyr, [Fe/H] ≈ −1.5 | 1.2% of the disk's density near the Sun |
| dust | thin exponential layer (scale height ~100 pc) plus clouds in the arms | — | 1 mag/kpc in the plane near the Sun |

Each component answers: its stellar density at a point, the age and metals of a star born
in it, and its share of stars per brightness bin (a luminosity table, built once per galaxy
by the importance sampling the sky uses today: kinds of star drawn with masses log-uniform
and ages within their lives, each weighted by the IMF and the history; unresolved
companions add their light).

Morphologies are component weights and shapes:

| Morphology | Components |
|---|---|
| Sa, Sb, Sc, Sd (each barred or not) | bulge share falling Sa → Sd; arms tighter in Sa, more open and more of them in Sc–Sd; gas and young stars rising Sa → Sd |
| S0 (lenticular) | bulge and disk, no arms, star formation long over |
| E (elliptical) | one Sérsic spheroid, old, metal-rich, no dust lane |
| dE / dSph | small faint spheroid, old, metal-poor |
| dIrr (Magellanic) | small clumpy disk, metal-poor, still forming stars |
| starburst | a disk with a recent burst |

The existing presets (`milky-way`, `elliptical`, `dwarf`, `starburst`) become morphologies
with components; `Galaxy` stays serde-loadable so a world can carry its own.

## Sectors

Stars are binned by absolute visual magnitude (one magnitude per bin, about −10 to +16;
remnants and brown dwarfs, which the eye never sees, in faint bins). A bin's reach is how
far its brightest stars stay brighter than magnitude 6.5 with no dust (Sun-like 22 pc, an
A star ~100 pc, a B star ~600 pc, a supergiant ~5 kpc), so a dusty sky is never cut short.

Space in the galaxy's frame (pc, origin at the centre, z toward its north pole) is cut into
cubes of side `C₀ / 2^ℓ` (`C₀` covering the galaxy). Bin `b`'s stars live only in cubes of
the level whose side is nearest half its reach, so a query for bin `b` around a point visits
about 5³ cubes. A cube's key is (level, i, j, k).

Per cube and bin, from the seed `galaxy · "stars" · bin · level · i · j · k`:

1. the expected number of stars, the components' densities times their share of the bin,
   integrated over the cube (vertically exactly for exponential layers, at sub-cube points
   otherwise);
2. that many on a Poisson draw;
3. each star's component (by its share in the cube), position (by the component's density
   within the cube, so a kiloparsec cube still puts its stars near the plane), and kind
   (from the component's luminosity table for the bin: mass, age, [Fe/H], companion).

A star's identity is (bin, level, i, j, k, n): stable, independent of who asks, and the seed
of its system. Queries: stars within a radius of a point down to an absolute magnitude, and
stars brighter than an apparent magnitude seen from a point (each bin out to the distance it
can be seen from there). A query generates only the cubes it reaches: zooming in on a
region generates its fainter stars as it goes.

### Young stars in groups

Most stars form in groups that disperse: embedded clusters hold 70–90% of the stars formed in
giant molecular clouds, and fewer than 4–7% of them survive to the Pleiades' age (Lada & Lada
2003). Stars younger than ~300 Myr are not field stars: they belong to birth events,
generated per kiloparsec cube where star formation is (the young thin disk, raised in the
arms, inside dust clouds):

- an event's mass from dN/dM ∝ M⁻² (10 to 10⁵ M☉), its age uniform over the 300 Myr;
- ~5% stay bound as an open cluster, compact and dissolving over 100 Myr to 1 Gyr; the rest
  spread as an association at a few pc per Myr;
- its members drawn from the IMF at its age, generated per bin on demand (seed: event · bin)
  so a far observer sees its O and B stars without generating its dwarfs.

The field's luminosity tables then cover ages over 300 Myr only, so no star is counted twice.

### Dust

Extinction along a ray is the smooth layer's (integrated analytically) plus that of the
clouds the ray passes (Gaussian blobs in the arm cubes, with masses from a cloud mass
function). The same clouds host the birth events and cut dark lanes in the glow.

### As built

- The galaxy-integrated IMF is steeper than the canonical one above a solar mass (α = 2.6;
  Scalo 1986's field slope of 2.5–2.7, Kroupa & Weidner 2003's integrated galactic IMF): with
  the canonical 2.3 the Sun's sky held 11,400 stars and twice the third- and fourth-magnitude
  ones; with 2.6 it holds 8,500–9,200 (9,110 observed), 18–26 brighter than 1.5 (22), and 320–
  360 systems within 10 pc (339).
- Groups under 100 M☉ (a quarter of the young mass under dN/dM ∝ M⁻² from 10 M☉) are young
  field stars, not events: a 10 kpc sky would otherwise draw two million events.
- A disk's plane is split once per cube for all its age slices: the arms' factor is linear in
  the contrast, so each leaf keeps its mass and its arm-weighted mass.
- In crowded places (a bulge) a sky lists stars to a brighter limit, about four Suns' worth,
  and the rest join the glow.
- A world made at a star from a standalone world (`geoforge visit`) gets the home system's
  star added to the galaxy in both worlds, so each sees the other.

## The night sky

The sky from a point: the stars brighter than 6.5 (any limit, for telescopes later) with
their dust, minus the observer's own system and, for a standalone system, generated stars
within 0.5 pc of it; directions turned into the system's frame by its orientation (a
rotation drawn from the system's seed). The glow is the unresolved light of every component
along each line of sight, through the dust, on a 1° grid of directions, cached.

Other galaxies in the sky (Andromeda, satellites) wait for the universe layers; the glow is
built as a sum of sources so they can be added.

## A system's place

`GalacticSite` gains the position in the galaxy (pc) and, for a system made at a star, the
star's identity. A standalone system draws its position from its world's seed (the thin
disk's density, as the radius is drawn today, plus height and azimuth). A system made at a
star takes the star's phenotype (primary mass, age, [Fe/H], companions) and its seed; the
habitable recipe keeps a Sun-like primary at a drawn place. Existing worlds get a position
at their stored radius, so their skies change but their systems do not.

## Edits

`galaxy/stars.json` in the world holds edits applied over the generated stars in every
query:

- **remove** a generated star by identity;
- **change** a generated star: name, position, kind; the whole new record is stored, so a
  changed star survives changes to the generator;
- **add** a star: a name, a position, and either a kind (mass, age, [Fe/H]) or a spectral
  type and magnitude;
- **add a cluster**: a centre, radius, number of members, age and metals, its members
  generated from its own seed (like a birth event) and editable in turn.

Positions are galactic (pc), or given in the focus's sky (right ascension, declination,
distance) and stored galactic. Added stars and cluster members get identities in an edit
namespace. An edit naming a star the generator no longer makes is reported, not an error.

Constellations belong to an observer, since patterns change with the point of view:
`constellations.json` beside the focus system holds named sets of lines between star
identities, which the sky chart draws.

Commands: `geoforge star add|move|change|remove`, `geoforge stars export|import` (the sky's
stars with identities, position, distance and magnitude as CSV or JSON, to edit in a
spreadsheet), `geoforge constellation add|remove`.

## Checks against the real sky

At the Sun's place in the Milky Way model: about 9,110 naked-eye stars (Yale Bright Star
Catalogue; today's smooth disk gives ~6,500), about 22 brighter than 1.5, giants the
commonest naked-eye stars; B stars in a belt of associations as the Gould Belt; the 10 pc
census of 540 stars, brown dwarfs and planets in 339 systems, about 61% M stars (Reylé et
al. 2021; brown dwarfs once the IMF reaches below 0.08 M☉); the band bright toward the bulge
and cut by dust lanes; the farthest naked-eye stars a few kiloparsecs off (V762 Cas,
~16,000 light-years).

## Cost

A naked-eye sky visits about 27 bins × 125 cubes and generates about four times the stars it
keeps (cubes reach past the sphere): ~40,000 stars, under a second. Cubes are pure functions,
so they can be cached and generated in parallel in a fixed order.
