# Geoforge Visualization Guide

This document describes the visualization outputs from the complete Stage 1 pipeline.

## Generated PNG Visualizations

When you run `cargo run --release --features export-png`, three visualization PNGs are created in the `outputs/` directory:

### 1. **tectonics.png** - Plate Boundaries
**What it shows:** Basic tectonic plate layout with each plate colored randomly.

**Use case:** Overview of plate sizes and distribution

**Details:**
- Each plate has a unique color (seeded random colors for consistency)
- Shows the final result after all Stage 1 processing
- Good for seeing relative plate sizes at a glance

---

### 2. **boundaries.png** - Boundary Classification ✨ NEW
**What it shows:** Plate boundaries color-coded by interaction type.

**Color Key:**
- 🔴 **Red** = **Convergent** boundaries (plates colliding)
- 🔵 **Blue** = **Divergent** boundaries (plates spreading apart)
- 🟢 **Green** = **Transform** boundaries (plates sliding past each other)

**Use case:** Understanding plate tectonics and future geological features

**Details:**
- Background: Grayscale plates for context
- Foreground: Color-coded boundary pixels
- Critical for Stage 2 (Geologic Provinces) - convergent boundaries become mountains!

**Example from latest run:**
- 15 convergent boundaries → Future mountain ranges
- 14 divergent boundaries → Future rift zones and mid-ocean ridges
- 25 transform boundaries → Future fault lines

---

### 3. **plate_motion.png** - Motion Vectors ✨ NEW
**What it shows:** Each plate colored by its motion direction and speed.

**Color Encoding (HSV Color Space):**
- **Hue (Color)** = Motion direction (0-360°)
  - 🔴 **Red** → Eastward (0°)
  - 🟡 **Yellow** → Northward (90°)
  - 🔵 **Cyan** → Westward (180°)
  - 🟣 **Blue/Magenta** → Southward (270°)

- **Saturation (Brightness)** = Motion speed
  - **Vivid/Bright** → Fast plates (8-10 cm/year)
  - **Grayish/Dull** → Slow plates (1-3 cm/year)

- **Value** = Fixed at 90% for visibility

**Use case:** Visualizing global tectonic flow patterns

**Details:**
- Each plate is uniformly colored by its motion vector
- Intuitive at-a-glance understanding of motion patterns
- Speed range: 1-10 cm/year (Earth-realistic values)

---

## Interpreting the Motion Visualization

### Example Scenarios

**Bright Red Plate:**
- Moving eastward at high speed (~8-10 cm/year)
- Example: Pacific Plate-like behavior

**Dull Yellow Plate:**
- Moving northward at slow speed (~2-3 cm/year)
- Example: Stable continental drift

**Cyan Plate:**
- Moving westward at medium speed
- Creates interesting collision zones with eastward plates

**Adjacent Complementary Colors (e.g., Red and Cyan):**
- Likely to form **convergent boundaries** (head-on collision)
- Will generate mountain ranges in Stage 2

**Adjacent Similar Colors:**
- Likely to form **transform boundaries** (sliding past)
- Will generate fault lines in Stage 2

---

## Real-World Comparison

Our generated world (latest run):
- **20 plates** (Earth has ~15 major plates + many smaller ones)
- **54 boundaries total**
- **Average velocity: 8.26 cm/year** (Earth's plates: 1-10 cm/year)
- **Plate size ratio:** 22.9% superplate to <1% microplates (realistic Earth-like distribution)

---

## Technical Details

### Resolution
- Default: 1800×900 pixels (0.2° per pixel)
- Represents a full planetary sphere in equirectangular projection

### File Sizes
- tectonics.png: ~56 KB
- boundaries.png: ~82 KB (slightly larger due to edge detail)
- plate_motion.png: ~53 KB
- world.map: ~3.1 MB (binary data with full metadata)

### Color Space
- Plate colors: RGB (random seeded)
- Boundary colors: RGB (categorical)
- Motion colors: HSV → RGB conversion
  - Ensures intuitive directional mapping
  - Natural perceptual relationship between hue and direction

---

## Using the Visualizations

### For Development
- **boundaries.png** - Verify classification algorithm works correctly
- **plate_motion.png** - Debug motion assignment and relative velocities
- **tectonics.png** - Quick sanity check on plate generation

### For Presentations
- **plate_motion.png** - Most visually striking, shows global dynamics
- **boundaries.png** - Clear demonstration of boundary classification
- All three together show complete tectonic system

### For Stage 2 Planning
- **boundaries.png** identifies where to place:
  - Mountain ranges (red convergent boundaries)
  - Rift valleys (blue divergent boundaries)
  - Fault systems (green transform boundaries)
- **plate_motion.png** shows relative velocities for mountain height estimation

---

## Example Output (Latest Run)

```
🎉 STAGE 1: TECTONIC FOUNDATION COMPLETE!

Statistics:
  • 20 plates generated
  • 54 boundaries identified
  • 15 convergent (mountain-building zones)
  • 14 divergent (rift/ridge zones)
  • 25 transform (fault zones)
  • Total boundary length: 538,499 km
  • Average relative velocity: 8.26 cm/year

Plate distribution:
  • 5 continental (large, stable)
  • 6 oceanic (denser, faster)
  • 9 mixed (composite)

Files:
  ✅ outputs/tectonics.png
  ✅ outputs/boundaries.png
  ✅ outputs/plate_motion.png
  ✅ outputs/world.map
```

---

## Next Steps

These visualizations provide the foundation for **Stage 2: Geologic Provinces**:

1. Use **convergent boundaries** (red) to place orogenic belts
2. Distinguish oceanic-continental vs continental-continental collision
3. Use **divergent boundaries** (blue) for rift zones and mid-ocean ridges
4. Use **transform boundaries** (green) for strike-slip fault systems
5. Use motion vectors to calculate mountain heights and volcanic activity

**Ready to build realistic geology on top of this tectonic foundation!** 🏔️
