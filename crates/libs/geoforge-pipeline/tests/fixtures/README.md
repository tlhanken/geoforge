# Known-case fixtures

Real-world inputs for the known-case tests (`tests/known_*.rs`). Built by
`scripts/known-cases/make_fixtures.py`; all grids are equirectangular with pixel-centred
cells, row 0 at the north pole.

| File | What | Source |
|---|---|---|
| `earth_etopo1_360x180.png` | Earth's topography, one-degree block means, 16-bit, metres + 12000 | ETOPO1 (NOAA NCEI, Amante & Eakins 2009) at 10′, Fatiando a Terra mirror; public domain |
| `pb2002_plates_720x360.png` | Plate index per half-degree cell, 8-bit, row order of `pb2002_plates.csv` | PB2002 (Bird 2003), GeoJSON by Hugo Ahlenius (github.com/fraxen/tectonicplates); ODC-By 1.0 |
| `pb2002_plates.csv` | Plate code, name, Euler pole relative to the Pacific (°N, °E, °/Myr) | PB2002 |
| `pb2002_steps.csv` | Boundary steps: midpoint, PB2002 class (SUB, OSR, OTF, …), plate pair, length | PB2002 |

## How the topography is conditioned

- **Block means** of a 6×6 sample of each one-degree cell.
- **Big rivers are burned in**: a block that a river draining ≥ 500 000 km² crosses takes the
  river's level, where that is lower than the block mean. Without this, block means dam
  valleys narrower than a cell. The Nile below Aswan, for one, would end in a closed basin.
- **Cut-off seas are raised**: below-sea land the ocean does not reach is raised to +1 m,
  unless it is a sea of 100 000 km² or more, like the Caspian.

The 10′ source still dams the Congo where it cuts the Crystal Mountains. The cuvette fills
as a lake, so no test asks whether the Congo reaches the sea.

## Attribution

Bird, P. (2003), An updated digital model of plate boundaries, *Geochem. Geophys. Geosyst.*
4(3), 1027, doi:10.1029/2001GC000252. Amante, C. & Eakins, B. W. (2009), ETOPO1 1 Arc-Minute
Global Relief Model, NOAA Technical Memorandum NESDIS NGDC-24.
