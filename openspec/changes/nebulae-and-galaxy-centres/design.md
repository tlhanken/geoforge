# Design: nebulae-and-galaxy-centres

## Should a nebula have an address of its own?

Everything in a universe has an address, so the question had to be settled: is a nebula a
place (a level or a kind of star reference, `7/h6:…/n:e3,1,0:12`), or only something drawn?

**For an address of its own**
- It can be shared, bookmarked and typed, like any place.
- A fictional material or a note could be attached to it by address.
- Search hits would land on it, not beside it.

**Against**
- A nebula is not a container: nothing is generated "inside" one the way a system holds
  bodies. Its stars, the thing one would visit, already have addresses.
- It is transient and derived: an H II region is a birth event's gas for 10 Myr, a
  remnant a dead star's shell for 100 kyr. Its identity is its maker's; a second identity
  duplicates it and must be kept in step.
- Nebulae overlap and nest (an H II region in its cloud), so "the nebula at a point" is not
  one thing.
- Address parsing, resolution and the store would all grow a level for something with no
  stored state.

**Neither, but between** (chosen): a nebula is selectable and named, with an identity
derived from its maker (`e<cube>:<event>`, `c<cube>:<n>`, `s<star>`), shown and searched
round the home star or a system, and it leads to the star that lights or left it. No level
is added; if nebulae ever gain stored state (edits), the derived identity can become an
address without changing how they are found.

Globular clusters are handled the same way for now: their light and place are drawn, and
their stars are not generated one by one.

## The just dead as bins

Planetary nebulae and supernova remnants need their stars. Such stars are about 10⁻⁵ of a
population's remnants, so a field table's few thousand kinds would almost never hold one, and
finding them kiloparsecs off by drawing every white dwarf is out of reach. Giving them
strata of their own, in bins whose cubes are 2 kpc, keeps them real stars of the field (with
addresses) that can be gathered cheaply. The dead stratum keeps its share of them (some 10⁻⁵
of its weight), so lit and dark stars are drawn exactly as before.
