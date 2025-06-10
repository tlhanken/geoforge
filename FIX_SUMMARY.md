# Geoforge - Fixed Spherical Tectonic Plate Generation

## The Fix

The issue with the messy output was caused by a fundamental problem in the `generate_plates_region_growing` method:

### The Problem
1. The frontier queue was storing 3D points (Point3D) rather than pixel coordinates
2. When popping from the frontier, the code was converting 3D points back to pixel coordinates
3. This conversion is not one-to-one - multiple 3D points can map to the same pixel
4. This caused pixels to be processed multiple times, creating the artifacts you saw

### The Solution
I've fixed the implementation by:

1. **Pixel-based Frontier**: Changed `GrowthPoint` to `GrowthPixel` which stores pixel coordinates (x, y) instead of 3D points
2. **Direct Pixel Processing**: The frontier now works entirely in pixel space, avoiding conversion issues
3. **Proper Distance Tracking**: Still uses geodesic distances for growth prioritization, but tracks pixels not 3D points
4. **Consistent Neighbor Addition**: The `add_pixel_neighbors_to_frontier` method now properly manages pixel-based growth

## Running the Fixed Version

To run with PNG output visualization:

```bash
# Compile with PNG export feature
cargo build --features export-png

# Run the program
cargo run --features export-png
```

The output files will be in the `outputs/` directory:
- `test.bin` - Raw binary data
- `test_world.bin` - World data
- `test_world.png` - Visual representation (if compiled with export-png feature)

## Key Changes in the Code

### Before (Broken):
```rust
struct GrowthPoint {
    point: Point3D,      // 3D point on sphere
    plate_id: u16,
    distance: f64,
}

// When processing:
let (lat, lon) = growth_point.point.to_lat_lon();
let x = ((lon + 180.0) / 360.0 * self.width as f64) as usize;
let y = ((90.0 - lat) / 180.0 * self.height as f64) as usize;
// Multiple 3D points could map to same pixel!
```

### After (Fixed):
```rust
struct GrowthPixel {
    x: usize,           // Direct pixel coordinates
    y: usize,
    plate_id: u16,
    distance: f64,
}

// When processing:
let idx = self.get_pixel_index(growth_pixel.x, growth_pixel.y);
// Each pixel processed exactly once
```

## Verification

The fixed version should produce:
- Clean blob-like tectonic plates without horizontal line artifacts
- Proper handling of poles without pinching
- Realistic plate size distribution (power law)
- Natural-looking fractal boundaries

## Notes

- I've backed up the broken version to `lib_spherical_broken.rs` for reference
- The spherical geometry calculations remain unchanged - they were correct
- The fix maintains all the improvements from SPHERICAL_IMPROVEMENTS.md
- Performance should be similar or slightly better due to reduced conversions
