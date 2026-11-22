# Geoforge v2 realism check

## Earth-like worlds

Six worlds at 2048×1024 (about 18 km cells at the equator), generated with
`geoforge-v2 new … && geoforge-v2 generate … && geoforge-v2 stats …`:

| World | Body |
|---|---|
| `earth42` | `--archetype earth`, seed 42 |
| `s2837` | Tevo IV: 5799 km, 0.89× Earth's starlight, 0.72 bar, day 30.5 h, obliquity 20.6° |
| `s42` | Kibedo IV: 5592 km, 0.88×, 0.84 bar, day 28.7 h, obliquity 18.4° |
| `s7` | Lisano IV: 6308 km, 0.93×, 0.74 bar, day 22.8 h, obliquity 23.3° |
| `s1234` | Vivi IV: 6278 km, 1.02×, 1.19 bar, day 36.8 h, obliquity 28.9° |
| `s99` | Shulor III: 6344 km, 0.91×, 1.46 bar, day 20.5 h, obliquity 29.1° |

A world takes about 37 s with the `iter` build (tectonics 7 s, geology 7 s, elevation with
erosion 17 s, circulation 1 s, climate 5 s, hydrology 2 s, biomes <0.1 s).

| Quantity | earth42 | s2837 | s42 | s7 | s1234 | s99 | Earth |
|---|---|---|---|---|---|---|---|
| Largest plate, of surface | 0.25 | 0.22 | 0.18 | 0.24 | 0.25 | 0.24 | 0.20 |
| Seven largest plates | 0.81 | 0.81 | 0.82 | 0.85 | 0.83 | 0.81 | 0.81 |
| Land fraction | 0.31 | 0.23 | 0.28 | 0.29 | 0.30 | 0.31 | 0.29 |
| Mean land elevation, m | 804 | 1482 | 754 | 1137 | 795 | 954 | 840 |
| Land above 1000 m | 0.29 | 0.76 | 0.15 | 0.65 | 0.15 | 0.23 | 0.27 |
| Highest point, m | 8171 | 10682 | 9549 | 9629 | 10472 | 10501 | 8849 |
| Mean ocean depth, m | 3773 | 3465 | 3811 | 3357 | 3684 | 4010 | 3680 |
| Surface rock: Cenozoic / Precambrian | 0.24 / 0.38 | 0.28 / 0.34 | 0.23 / 0.38 | 0.25 / 0.36 | 0.27 / 0.39 | 0.29 / 0.38 | 0.28 / 0.32 |
| Global mean temperature, °C | 14.8 | 9.2 | 10.0 | 10.9 | 16.8 | 13.5 | 14.5 |
| Equator (5°S–5°N), °C | 26.0 | 21.8 | 22.0 | 22.2 | 25.3 | 21.9 | 26.5 |
| 40–50°N, °C | 10.2 | 6.7 | 5.8 | 7.5 | 12.4 | 10.5 | 10 |
| 55–65°N, °C | −1.8 | −8.0 | −10.3 | −9.0 | 6.0 | 3.6 | 0.5 |
| Annual range, land 40–60°, °C | 30.3 | 34.9 | 23.8 | 32.2 | 34.7 | 18.8 | ~28 |
| Annual range, ocean 40–60°, °C | 5.8 | 5.6 | 5.0 | 6.9 | 6.6 | 3.5 | ~7 |
| Global precipitation, mm/yr | 1024 | 927 | 897 | 929 | 1089 | 980 | 1000 |
| Land precipitation, mm/yr | 851 | 743 | 628 | 896 | 1066 | 631 | 780 |
| Land drier than 250 mm/yr | 0.39 | 0.49 | 0.50 | 0.36 | 0.33 | 0.45 | ~0.2 |
| Precipitation 5°S–10°N, mm/yr | 2264 | 2321 | 2607 | 1749 | 2241 | 975 | ~2000 |
| Precipitation 40–55°N, mm/yr | 1023 | 943 | 722 | 1181 | 670 | 1644 | ~1000 |
| River discharge to the sea, km³/yr | 54 473 | 40 074 | 30 901 | 62 961 | 63 756 | 41 064 | 37 000 |
| Rivers over 10 000 m³/s | 24 | 20 | 17 | 31 | 24 | 26 | 8 |
| Lakes over 10 000 km² | 5 | 4 | 5 | 3 | 6 | 16 | 16 |
| Endorheic land | 0.08 | 0.01 | 0.06 | 0.03 | 0.01 | 0.14 | 0.18 |
| Mediterranean biome, of land | 0.046 | 0.022 | 0.041 | 0.039 | 0.028 | 0.055 | 0.02 |
| Hot desert, of land | 0.20 | 0.05 | 0.17 | 0.16 | 0.20 | 0.34 | 0.14 |
| Temperate grassland, of land | 0.14 | 0.15 | 0.11 | 0.12 | 0.11 | 0.13 | 0.08 |

Figures are for the code as of the real-Earth check below (climate v3, circulation v3,
biomes v3).

The generated worlds are not Earth: most get 0.88–0.93× Earth's starlight and differ in
size, air pressure, day length and tilt. The carbonate–silicate thermostat (set in the solar
system step) puts CO₂ where the world weathers as fast as Earth with Earth's outgassing
(Walker, Hays & Kasting 1981: rate ∝ pCO₂^0.3·e^(ΔT/13.7 K)), so dim worlds build up CO₂
and run a few degrees cooler, bright ones draw it down and run warmer. At Earth's spin:
0.8× Earth's light, 7300 ppm and a 5 °C mean; 0.88×, 2800 ppm and 8 °C; 1.1×, 50 ppm and
24 °C. It balances weathering, which runs on warm unfrozen ground, not the mean
temperature, using the same latitude-band energy balance as the climate (heat transport ∝
day^(2/3), Kaspi & Showman 2015): a 6-hour day keeps about Earth's CO₂ (350 ppm) but runs
8 °C colder with poles at −39 °C; a tidally locked Earth, weathering only by day (Edson et
al. 2012), builds up ~2700 ppm with 32 °C under the star. Where weathering jumps at the
ice-albedo tipping point the thermostat settles on the warm side. (Until this check it
cancelled a fixed 75% of the weathering difference from Earth, which fails across the
tipping point: seed 42's world, at 0.88× Earth's light, sat on it and froze over, a −41 °C
snowball; it now runs at 7 °C. The balance it solves also used 10° bands, which hold partly
frozen states that the climate's finer bands do not; it now uses 5°.) Thick air damps the seasons (s99,
1.46 bar: 20 °C land range at 40–60°); thin air widens them.

Drainage networks: on the Earth-like world at 1024×512, the exponent of the drainage-area
distribution (P(A ≥ a) ∝ a^−β) went from 0.66 without erosion to 0.54 with it; Earth's
networks give about 0.43–0.45.

## Other worlds

Presets (`--archetype`, real sizes, orbits and air) at 2048×1024, seed 42. Elevation is
relative to the mean radius.

| Quantity | Mars | Venus | Moon | Mercury |
|---|---|---|---|---|
| Generators | lid ×3, circulation, climate | lid ×3, circulation, climate | lid ×3 | lid ×3 |
| Lowest / highest, km | −9.4 / +20.1 | −3.6 / +12.0 | −10.2 / +5.4 | −7.6 / +5.5 |
| Real, km | −8.2 / +21.2 | about −3 / +11 | 19.9 km range | −5.4 / +4.5 |
| Surface within 1 km of the mean | 0.13 | 0.80 (real: over 0.80) | 0.27 | 0.17 |
| Lowland (thin crust) share | 0.35 | 0.80 | 0.18 (maria ≈ 0.17) | 0.27 |
| Global mean temperature, °C | −65 (real ≈ −63) | 466 (real 464) | — | — |
| Equator / north pole, °C | −48 / −115 | 467 / 459 | — | — |
| Seasonal range, 40–60°, °C | 71 | 2 | — | — |

The climate model's greenhouse is fitted to Earth, Mars, Venus and Titan presets
(`cargo run -p geoforge-core --example radiation_check`, and
`-p geoforge-surface --example zonal_check` for zonal profiles).

Generated systems also give Venus-like inner worlds (runaway greenhouse, 35–360 bar CO₂),
cold super-Earths with thin CO₂, ice-shell ocean moons around giants, frozen cratered moons
and belt members; `geoforge-v2 show` lists each body with its regime, air and water.

## Real Earth

`scripts/real-earth/run.sh` puts Earth's actual topography (ETOPO1, as authored elevation)
under an `--archetype earth` world. It generates circulation, climate, hydrology and biomes,
then compares each layer with the observed Earth: ERA-Interim winds, ERSST v5 sea
temperatures, the Köppen–Geiger map and 31 station normals (see the script's README). With
the same land and sea as Earth, any difference comes from the model, not from the map.

The first run answered the grassland question. On Earth's own map the model gave 2% hot
desert and 20% temperate grassland. The cause was the climate, not where the land lay: the
Sahara, Arabia, Australia and the US Southwest came out wet, and the Indian monsoon was
missing. The changes since:

- **The ITCZ moves with the continents.** It follows the sun over the oceans but moves into
  the summer hemisphere over large heated land. The trades that cross the equator turn into
  monsoon westerlies. Rain falls on the equatorward side of that trough.
- **Monsoons flow into subtropical heat lows**, turned by the Coriolis force and stopped by
  high ranges (the Andes keep the Amazon from pulling the Pacific ashore at Peru).
- **Sinking air dries:**
  - a broad Hadley sinking belt, strongest in the winter hemisphere;
  - Sverdrup balance on the flanks of the subtropical highs: air turned equatorward sinks,
    and that sinking reaches inland with the wind;
  - northerlies on the poleward side of a heat low;
  - stable air off cold water.
- **Rain needs more humid air where it sinks** (30% to 70% of saturation).
- **Summer land under sinking air is hotter** (up to 8 °C in the moisture model).
- **Other moisture and heat fixes:**
  - air over colder ground condenses gradually instead of all at once (this had dumped the
    Pacific on Southern California);
  - storm-track eddies carry vapour deeper inland;
  - land recycling scales with the season's warmth, but land gives back little more than open
    water would evaporate into the same air (wet regions had fed on their own rain);
  - lifting over terrain wrings out the vapour within a scale height or two (~1.5 km), so a
    4 km plateau rim rains about twice what a 1 km rise does, not four times;
  - broad plateaus warm the air above them (about 2 °C per km);
  - maritime air reaches about 2500 km inland instead of 1200;
  - deserts get +5 °C mean and +7 °C range instead of +3.5 and +5.
- **Ocean:** the western boundary current crosses its basin where the gyres meet and runs up
  the far shore (a North Atlantic Drift), but only in a basin of a gyre's width and only up a
  coast that runs north–south. The SST anomaly remembers for months, and upwelling is
  stronger. Only flow that differs from its latitude's mean makes an anomaly, the gradient it
  crosses stops at the freezing point of seawater, and the zonal mean taken off is smoothed
  over 10° of latitude (each had drawn warm and cold bands round the globe).
- **Biomes:** ground that never thaws is ice sheet, not alpine, even on a high plateau.

| Real Earth, 2048×1024 | Before | After | Observed |
|---|---|---|---|
| Köppen major group, cell agreement | 0.49 | 0.58 | — |
| Desert (BW), of land | 0.09 | 0.21 | 0.18 |
| Hot desert (BWh), of land | 0.02 | 0.16 | 0.14 |
| Steppe (BS), of land | 0.18 | 0.12 | 0.12 |
| Temperate (C) / continental (D), of land | 0.25 / 0.13 | 0.17 / 0.14 | 0.15 / 0.22 |
| Hot desert biome / temperate grassland | 0.02 / 0.20 | 0.16 / 0.09 | 0.14 / 0.08 |
| Ice sheet / alpine biome | 0.04 / 0.07 | 0.10 / 0.01 | 0.10 / — |
| Station RMSE: mean temperature, °C | 5.4 | 4.3 | — |
| Station RMSE: log precipitation | 1.89 | 1.45 | — |
| Wind direction vs ERA-Interim, Jan / Jul | 0.75 / 0.70 | 0.77 / 0.65 | — |
| July wind, northern Indian Ocean (u, v), m/s | (−3.8, 3.3) | (4.2, 4.4) | (12.6, 2.9) |
| SST departure from zonal mean, correlation | 0.27 | 0.38 | — |
| NE Atlantic SST departure, °C | +0.1 | +1.2 | +4.4 |
| Zonal temperature 25–35°N, °C | 17.1 | 18.4 | 20 |

Stations now in the right Köppen group include:
- deserts: Cairo, In Salah and Alice Springs;
- the Mediterranean coasts: Los Angeles and Perth;
- Mumbai, Manaus and Singapore (tropical).

Known misses on Earth's map:

1. **Peru's coast is wet** (1600 mm vs 15). The cold current is a third as cold as observed,
   and the southern summer's rain belt still reaches it.
2. **The East African highlands are drenched** (Nairobi 5100 mm vs 870). There is no Walker
   circulation, and nothing makes the Somali coast's flow diverge.
3. **The Indian monsoon is a third as strong as observed.** Mumbai gets 1040 mm vs 2200.
4. **Eastern continental interiors at 35–45° are dry**: Chicago 430 mm vs 950, Buenos
   Aires 375 vs 1200, Beijing 350 vs 570. Moist poleward flow from the Gulf of Mexico, the
   South Atlantic and the South China Sea is too weak. Siberia is dry too (Yakutsk 51 mm vs
   240).
5. **The North Atlantic Drift brings a quarter of the observed warmth.** London is 7.7 °C
   vs 11.3, Reykjavik −1.4 vs 4.7.
6. **The Antarctic plateau is 15 °C too warm.** The zonal energy balance does not see the
   ice sheet's albedo, and there is no winter inversion. Yakutsk's winters are too mild for
   the same reason.
7. **Subtropical oceans get half their rain** (20–30°: 360 mm vs ~700).
8. **Rivers can't be checked on this map.** ETOPO1 at 20 km blocks narrow gorges, so the
   Congo ponds behind its rapids, and the upper Mississippi and middle Yangtze end in closed
   basins. The known cases below burn big rivers into their map and check rivers there.

## Known cases

Three `cargo test` suites give each layer real inputs and check its output against what the
real body does. The inputs are small fixtures in
`crates/libs/geoforge-pipeline/tests/fixtures`. Cases the model is known to miss are
recorded as gaps, each with its reason. Gaps never fail the main tests, and each suite's
`#[ignore]`d `…_known_gaps` test asserts them:

```
cargo test -p geoforge-pipeline known_cases -- --nocapture   # the three scorecards
cargo test -p geoforge-pipeline --test known_planets          # one suite
cargo test -p geoforge-pipeline known_gaps -- --ignored       # the open gaps
```

Each run writes its scorecard to `target/known-cases/<suite>.txt`.

| Suite | Inputs | Checked against | Cases | Gaps |
|---|---|---|---|---|
| `known_earth_tectonics` | PB2002 plates and Euler poles over the hotspots, half-degree; crust from ETOPO1; then circulation, climate and erosion on the model's own heights | PB2002 boundary classes; crust ages; provinces at named places; ETOPO1 heights; earthquake rates and corner magnitudes | 58 | 7 |
| `known_earth_climate` | ETOPO1, one degree, as authored elevation | winds and currents; SST departures; 31 station normals; Köppen groups and full classes; closed basins and the Nile; biome shares and 20 places; the biome and Köppen classifiers fed station normals directly; cyclone tracks; sea-level pressure | 223 | 51 |
| `known_planets` | generated Mars and Moon, median of seeds 1–3 | MOLA, LOLA, InSight and the CO₂ frost point | 15 | 4 |

Biggest open gaps, by layer:

- **Elevation from plates:**
  - The climate on the model's own heights misplaces the rain (the Himalayan front gets
    170 mm, the central Andes 1900 mm), so rain-aware erosion leaves the Himalaya too high
    and wears the wet-in-the-model Andes and Cascades too low.
  - No flat-slab segments: the Andean plateau runs on south of the Altiplano over the
    Pampean ranges (2.7 km, ETOPO1 0.7), and Peru's cordillera is half as wide as it should be.
  - No plume uplift (Ethiopia), no intraplate ranges (the Rockies); the Amazon lowland
    beyond the Andean foreland sits at the continental base (~300 m).
  - Ridge–transform staircases are smoothed into ridges.
- **Rain:**
  - Winter storms still drench the interior Southwest (Phoenix 930 mm, and the Great Basin
    drains to the sea), London and Kazakhstan.
  - The West African monsoon reaches too far north (Timbuktu 810 mm); Nairobi gets 2800 mm.
  - North-east Siberia is too dry (Yakutsk 15 mm).
- **Temperature:** the subtropical interiors are ~5 °C too cold (Phoenix, Riyadh, Timbuktu),
  Europe turns continental too fast, Siberian winters are too mild, Lima's coastal cold
  does not reach the air, and the Antarctic plateau is ~17 °C too warm.
- **Biomes:**
  - At one degree, every lowland cell on a big river turns wetland.
  - The classifier gets four of 31 stations wrong: Tehran (savanna), Nairobi and Reykjavik
    (temperate forest), and Vancouver (just short of rainforest).
- **Mars and the Moon:**
  - Polar winters fall below the CO₂ frost point.
  - Mars's relief overshoots in two of three seeds.
  - The Moon has no far-side highland bulge and no basin the size of South Pole–Aitken.

### Continental collisions

Collisions are built asymmetrically (`geology::provinces::Collision`):
- **Which plate goes under:** the smaller continent. India and Arabia went under Asia
  because they rode in behind ocean floor that was already subducting beneath it. Once their
  buoyant upper crust is scraped off, the rest is dense enough to follow. The exception is a
  plate that already rides over an ocean slab within ~600 km: its crust is hot and weak, so
  the other drives under it (in Taiwan, China's margin goes under the Luzon arc).
- **What forms:** a foreland basin on the lower plate, a steep front range across the
  boundary, and a plateau behind it on the upper plate. Rivers do not cross the plateau's
  interior to the sea, so it stays flat and dotted with lakes.
- **Speed:** sets the range's height and the plateau's width and height. Collisions
  slower than 0.5 cm/yr build nothing; a plateau needs 2 cm/yr or more.
- **Gravity:** heights scale as 1/g (`relief::relief_scale`), the load a crust can bear.

On PB2002's plates this gives Tibet 4.95 km (ETOPO1 4.95), the Himalaya 3.1 km (3.1), the
Ganges plain 230 m (200), Mesopotamia 150 m (185), and the Iranian plateau 1.9 km (1.3).
Iran is too high because rate is all the model knows: Arabia converges almost as fast as
India but has collided half as long.

### Ocean–continent margins

Where ocean floor goes under a continent (`geology::andean_margin`), the margin is built as
the Andes, Cascadia and Sumatra are, not as an island arc:
- **Across the margin:** a narrow wedge under a steep slope with a shelf a few tens of km
  wide, so the coast lies ~100 km from the trench and runs straight along it; an onshore
  forearc (the coastal range and central valley); then a volcanic cordillera 150–250 km
  wide, ~1.5 km high on slow margins and ~4 km on fast ones. Behind it the upper plate is
  shortened, so a retroarc foreland basin replaces the extensional backarc.
- **The Andean plateau:** where convergence is fast (4 cm/yr or more) and the slab entering
  the nearest trench is 50–100 Ma old, a plateau three quarters of Tibet's height stands
  behind the arc. Capitanio et al. (2011) tie the Altiplano to the Nazca slab being oldest
  there; slabs much older than ~100 Ma roll back and stretch the upper plate instead
  (Molnar & Atwater 1978).
- **Rifts** do not cut across belts already built at convergent margins.

On PB2002's plates the central Andes (15–25°S, 70–65°W) average 2.5 km (ETOPO1 3.1), and the
median arc–trench gap is 260 km (Earth 150–300).

## What matches

- Plate sizes (largest ~0.2–0.25, seven largest ~0.81–0.85); continents inside plates with
  passive margins.
- Ocean depth by crust age, land fraction and hypsometry for Earth-like water; surface rock
  ages close to Earth's mix.
- Rivers carve dendritic valleys; mountain belts stay high.
- Western boundary currents run warm and eastern margins cold with upwelling; zonal
  temperatures within ~2 °C of Earth from the equator to 50° on the Earth preset; land/ocean
  seasonal ranges; global rainfall and the ITCZ; summer-dry west coasts in the subtropics
  (Mediterranean 2–6% of land).
- On Earth's own map: the deserts of North Africa, Arabia, Australia and Central Asia; the
  Mediterranean climates of California, Chile, the Cape and south-west Australia; the
  tropical rainforest belt; the ice sheets.
- Mars, Venus, the Moon and Mercury: relief ranges, Venus's narrow hypsometry, Mars's
  dichotomy and volcanic rises, mean temperatures of Mars and Venus.

### Köppen view

The Köppen–Geiger view (a derived view, `render <dir> koppen`) classifies the stored
climate by Peel et al. (2007)'s rules. On the one-degree Earth run, land outside
Antarctica splits A 19%, B 36%, C 20%, D 15%, E 11%; Peel et al.'s map gives A 19.0, B
30.2, C 13.4, D 24.6, E 12.8. Hot desert (BWh) is 18% (observed 14.2%), savanna (Aw) 10%
(11.5%). D is short because Siberia and Canada come out dry and partly winter-dry: Dw is 9%
of land, where Earth's Dw is mostly East Asia. The Siberian dryness is a listed climate gap.
The model puts 13 of 31 stations in their full class.

Fed station normals directly, the classifier gets 18 of 31. Of the 13 misses, 9 sit on a
class line with rounded normals (Lhasa's 430 mm against a 440 mm aridity line, Mumbai's
dry season against Am's 12 mm). The other 4 show what the rebuilt months cannot hold. The
climate stores means, ranges and half-year rain shares, so months are sinusoids in phase.
That cannot show two rainy seasons (Nairobi), a short tropical dry season (Kinshasa), or a
January colder than the mean less half the range (Beijing).

### Pressure

Sea-level pressure recovered from the modelled winds on the ETOPO-forced Earth run
(`known_earth_climate`), against ERA5's 1991–2020 seasonal means: within 8 hPa at 10 of 12
centres: in January the Aleutian (997 hPa) and Icelandic (994) lows, the Azores (1018) and
South Pacific (1022) highs and the Australian heat low (997); in July the Azores–Bermuda
(1029) and North Pacific (1028) highs, the Asian heat low (1005), and the South Atlantic and
Australian highs. Gaps: the Siberian high comes out 1020 hPa (observed ~1036), because the
modelled Siberian winter is 5–15 K too warm; the Antarctic trough is 1000 hPa (observed
~987), because the modelled westerlies peak at 40–45°S, not 50–55°S.

### Hazards

On Earth's real plates (`known_earth_tectonics`): 106 magnitude 6+ quakes a year worldwide
(USGS ~140); 0.03 a year per 10,000 km² round Japan (JMA's ~17 a year over 2.6×10⁶ km² is
0.065), 0.02 in central Chile, under 0.002 in central Australia, the Siberian craton and the
Brazilian shield; corner magnitudes 9.5 in Chile, 8.5 in the Himalaya and on the San Andreas
(Fort Tejon 1857: 7.9). On Earth's modelled climate (`known_earth_climate`): cyclones pass
0.7 a year per 10,000 km² in the Philippine Sea, 0.6 in the Gulf of Mexico, 0.5 in the Bay
of Bengal, none off Peru or on the equator. Gap: 0.25 in the South Atlantic, where the
modelled summer sea is too warm and wet and wind shear (which stops real storms there) is
not modelled. Resource rates have no Earth check yet (see the roadmap).

### Magnetic fields

`geoforge_core::magnetism`'s tests: Earth 30 µT, magnetopause 9.5 radii (observed ~10),
auroral oval 71° (observed 65–75°); Mars's dynamo stopped ~4.1 Gyr ago, Venus, the Moon and
Mercury have none; Jupiter 420 µT with a 40-radius magnetopause (a pure dipole; Io's plasma
inflates the real one to 50–100); ice giants multipolar and steeply tilted. Missed: Mercury's
weak field, Ganymede's dynamo, Saturn's weakness for its mass, and moons inside a giant's
magnetosphere.

### Star systems

`geoforge-cosmos`'s unit tests check systems drawn from the Milky Way model against the
Galaxy's census (`census` prints the same numbers for any galaxy and range of seeds):

| Check | Model | Milky Way |
|---|---|---|
| Red dwarfs among stars | 75% | ~75% of stars (61% of stars and brown dwarfs within 10 pc, Reylé et al. 2021) |
| White dwarfs | 6.1% | 4.8×10⁻³ pc⁻³ of ~0.08 stars pc⁻³ (Holberg et al. 2016) |
| Neutron stars, black holes | 0.55%, 0.17% | 10⁸–10⁹ of each among 1–4×10¹¹ stars |
| Multiples | 25% of red dwarfs, 46% of Sun-like stars | 23–27% (Winters et al. 2019), 46% (Raghavan et al. 2010) |
| Metallicity near the Sun | mean −0.1, 0.2 dex wide | peak −0.15 to 0 (Gilmore 2007; Haywood 2001) |
| Small planets per red dwarf within 200 days | 1.5–3.5 | 2.5 (Dressing & Charbonneau 2015) |
| Sun-like stars with giant planets; hot Jupiters | 10–18%; 0.3–2% | ~14%; ~1% |
| Pulsars with planets | 0.5% | under 0.5% (Nițu et al. 2022) |
| Neighbouring planets, mutual Hill radii | at least 8 | Jupiter–Saturn ~8; long-lived multi-planet systems 8–10 (Chambers et al. 1996) |
| Asteroid and Kuiper belts | the first giant's 4:1–2:1; the outermost giant's 3:2–2:1 | Jupiter's 2.06–3.28 AU; Neptune's 39.4–47.8 AU |

Their gaps:

- No brown dwarfs, and no binary evolution: mass transfer, mergers and X-ray binaries are
  missing, and a supernova that blows off half a pair's mass always unbinds it (in fact
  ~8% of black holes keep a companion).
- Only the primary (or a close pair) has planets; companions' own planets are not made.
- A giant's light turns its rocky planets into Venuses, but icy moons do not thaw: icy
  interiors ignore starlight.
- Red dwarf worlds get the same air and water rules as Sun-like ones: no flare stripping
  and no water lost while the young star was bright.
- A pulsar's wind does not heat its planets.

### Other worlds' circulation

- Cells beyond Earth's three are a pattern, not a simulation: very fast rotators (5–7 cells)
  show thin row-by-row streaks in ocean rain, and the ocean gyres still follow one jet.
- Tidally locked worlds: the substellar point is always at 0°E; ocean currents do not carry
  heat to the night side (the SST anomaly is left out).
- The thermostat assumes land is spread evenly; a continent under the star would weather
  more and hold less CO₂ (Edson et al. 2012), and Kite et al. (2011)'s runaway substellar
  weathering is not modelled.
- Giants' bands are drawn from the wind profile, not simulated: no polar hexagon or
  cyclone clusters, storms do not drift, and colour is a palette by class (no chemistry).

## Known gaps

1. **Generated worlds now run dry.** Hot desert is 0.05–0.34 of land vs 0.14, land drier
   than 250 mm is 0.33–0.50, and the 20–30° belt gets 250–630 mm vs ~800. Their continents
   are larger and less broken up than Earth's, which leaves more interior. On Earth's own
   map the desert share is close (0.16 vs 0.14). Grassland is 0.11–0.15 vs 0.08.
2. **Too many large rivers** (8–38 over 10 000 m³/s vs 8): continents here are bigger and
   less broken up than Earth's, and runoff is somewhat high.
3. **Drainage networks branch less than Earth's** (0.54 vs ~0.44); longer erosion flattens
   the land without improving it, so the remaining gap is likely resolution.
4. **Some worlds are very high** (s2837: 76% of land above 1000 m): a broad plateau term can
   lift a whole large continent. `plateau_m` controls it.
5. **Ocean currents are a diagnostic gyre model**, not a simulation.
6. **Closed basins depend on terrain** (endorheic land 0–15% vs Earth's 18%).
7. **The gray greenhouse warms the tropics more than the poles** for added CO₂ (real CO₂
   forcing is flatter), so dim, CO₂-rich worlds have cold high latitudes.
8. **Mercury's basins come out deep** (−7.6 km vs −5.4 km), and the Moon's highest highlands
   low (+5.4 vs +10.8 km): basin depth does not yet depend on the lithosphere's heat.
9. **No generator yet for ice-shell surfaces** (Europa, Titan) or giants; their layers report
   why they are skipped. Titan's methane cycle is out of scope.
