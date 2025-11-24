# Subduction Zone Geological Provinces - Color Reference

## Base Layer Colors

### Oceanic Plates
- **Ocean Blue `[70, 130, 180]`** = Abyssal Plain (oceanic crust base layer, ~5000m depth)

### Continental Plates
- **Bright Yellow `[255, 255, 0]`** = Platform (sedimentary-covered continental basement)
- **Orange `[255, 140, 60]`** = Craton/Shield (exposed Precambrian basement)

## Boundary Line Colors
- **Red** = Convergent boundaries (plates colliding)
- **Blue/Cyan** = Divergent boundaries (plates spreading apart)
- **Green** = Transform boundaries (plates sliding past)

---

## OCEANIC-OCEANIC SUBDUCTION (Island Arc System)

### Subducting Plate (oceanic) - Going DOWN
1. **Dark Navy Blue `[0, 50, 100]`** = Ocean Trench
   - Location: Right at the convergent boundary
   - Width: ~50 km (very narrow, just boundary pixels)
   - Depth: 7-11 km below sea level
   - Example: Mariana Trench

### Overriding Plate (oceanic) - Staying UP
Sequence from boundary moving inland:

1. **Tan/Beige `[255, 200, 100]`** = Accretionary Wedge
   - Location: Adjacent to boundary on overriding side
   - Width: 100-200 km
   - Description: Scraped sediments piled up from subducting plate

2. **Medium Gray `[160, 160, 160]`** = Forearc Basin
   - Location: Behind accretionary wedge
   - Width: 150-300 km
   - Description: Subsided basin between wedge and arc

3. **Bright Red `[255, 50, 50]`** = Volcanic Arc
   - Location: 300-600 km from trench
   - Width: 200-400 km
   - Description: Island arc volcanoes (e.g., Aleutian Islands, Mariana Islands)

4. **Light Gray `[180, 180, 180]`** = Backarc Basin
   - Location: Behind volcanic arc
   - Width: 200-400 km
   - Description: Extensional basin with active spreading
   - Example: Sea of Japan, Mariana Trough

5. **Ocean Blue `[70, 130, 180]`** = Abyssal Plain (base layer continues)

---

## OCEANIC-CONTINENTAL SUBDUCTION (Andean-Type)

### Subducting Plate (oceanic) - Going DOWN
1. **Dark Navy Blue `[0, 50, 100]`** = Ocean Trench
   - Location: Right at the convergent boundary
   - Width: ~50 km
   - Example: Peru-Chile Trench

### Overriding Plate (continental) - Staying UP
Sequence from boundary moving inland:

1. **Tan/Beige `[255, 200, 100]`** = Accretionary Wedge
   - Location: At continental margin, adjacent to boundary
   - Width: 100-200 km
   - Description: Scraped sediments at plate edge

2. **Medium Gray `[160, 160, 160]`** = Forearc Basin
   - Location: Between wedge and arc
   - Width: 150-300 km
   - Example: Central Valley of Chile

3. **Bright Red `[255, 50, 50]`** = Volcanic Arc
   - Location: 300-600 km from trench
   - Width: 200-400 km
   - Description: Continental arc volcanoes
   - Example: Andes Mountains

4. **Light Gray `[180, 180, 180]`** = Backarc Basin
   - Location: Behind volcanic arc (if present)
   - Width: 200-400 km
   - Description: Only forms if extension occurs behind arc
   - Note: Not always present in continental settings

5. **Bright Yellow `[255, 255, 0]`** = Platform
   - Location: Continental interior beyond arc system
   - Description: Stable continental basement

---

## Other Divergent Boundary Features

### Mid-Ocean Ridges (Divergent Boundaries)
- **Cyan-Green `[100, 200, 200]`** = Mid-Ocean Ridge
  - Location: At divergent boundaries (blue lines)
  - Width: Variable by spreading rate (60-150 km)
  - Description: Active seafloor spreading center

### Large Igneous Provinces
- **Purple `[150, 50, 150]`** = Continental Flood Basalts
- **Light Purple `[180, 100, 180]`** = Oceanic Plateaus
- **Pale Purple `[200, 120, 200]`** = Hotspot Tracks

### Continental Margins
- **Pale Yellow `[255, 240, 120]`** = Extended Crust / Passive Margins
  - Location: Edges of continental plates (not at subduction zones)
  - Width: 130-350 km

---

## Current Width Settings (as of latest commit)

| Feature | Min Width (km) | Max Width (km) | Notes |
|---------|---------------|----------------|-------|
| Ocean Trench | 50 | 50 | Just boundary pixels |
| Accretionary Wedge | 100 | 200 | |
| Forearc Basin | 150 | 300 | |
| Volcanic Arc | 200 | 400 | Plus 300-600 km offset from trench |
| Backarc Basin | 200 | 400 | |

**Total Subduction Zone Width:** ~750-1350 km from trench to end of backarc

---

## Notes on Polar Regions

At high latitudes (>60° N/S), features that span the full image width (x=0-1799) are NOT wrapping around the entire planet. Due to equirectangular projection compression:

- At 80° latitude: Full image width = only ~6,950 km actual circumference
- At the poles (90°): All x-coordinates represent the same point

This means large backarc basins at polar regions can legitimately appear to cover "the entire width" while actually being appropriately sized in 3D space.
