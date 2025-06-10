# Spherical Plate Generation Improvements

## Overview

The plate generation algorithm has been completely rewritten to work on a proper sphere, eliminating the pinching artifacts and creating more realistic tectonic plates.

## Key Changes

### 1. **3D Spherical Coordinate System**
- Added `Point3D` struct to represent points on a unit sphere
- All distance calculations now use geodesic (great circle) distances
- Seed placement and region growing happen in 3D space
- Only project to 2D map at the final step

### 2. **Geodesic Distance Calculations**
```rust
fn geodesic_distance(&self, other: &Point3D) -> f64 {
    // Dot product gives cos of angle between unit vectors
    let dot = self.x * other.x + self.y * other.y + self.z * other.z;
    let dot = dot.max(-1.0).min(1.0);
    dot.acos()
}
```

### 3. **Variable Plate Sizes**
- Implements power-law distribution similar to Earth's plates
- Larger plates are much bigger than smaller ones
- Controlled growth based on target sizes

### 4. **Spherical Region Growing**
- Uses priority queue with geodesic distances
- Growth happens on sphere surface, not distorted 2D projection
- Eliminates pinching at poles

### 5. **Realistic Boundary Complexity**
- Added fractal noise to boundaries after generation
- Multiple octaves of noise for natural-looking edges
- Geodesic-aware smoothing preserves spherical geometry

## Technical Details

### Uniform Sphere Sampling
Seeds are placed using uniform sphere distribution:
```rust
let u = rng.gen::<f64>();
let v = rng.gen::<f64>();
let theta = 2.0 * PI * u;
let phi = acos(2.0 * v - 1.0);
```

### Priority Queue Growth
Region growing uses a min-heap priority queue ordered by geodesic distance:
```rust
struct GrowthPoint {
    point: Point3D,
    plate_id: u16,
    distance: f64,
}
```

### Size-Controlled Growth
Growth probability adjusts based on current vs target size:
```rust
let growth_probability = if current_fraction < target_fraction {
    0.95 // High probability if under target
} else {
    0.95 * (target_fraction / current_fraction).min(1.0)
};
```

## Results

With these changes, you should see:
- **No pinching at poles** - Plates maintain blob-like shapes
- **Realistic size distribution** - Mix of large and small plates like Earth
- **Natural boundaries** - Fractal, complex edges instead of smooth curves
- **Proper polar behavior** - Plates at poles behave like Antarctic plate

## Usage

The API remains unchanged:
```rust
let mut generator = TectonicPlateGenerator::new(3600, 1800, 15)?;
generator.generate("region_growing", true)?;
```

But the results are now spherically correct!
