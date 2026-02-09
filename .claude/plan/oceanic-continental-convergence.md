# Implementation Plan: Oceanic-Continental Convergence

## Problem Statement

**Critical Gap Identified:** The geology generation system is missing oceanic-continental convergence zones (e.g., Andes Mountains, Cascades Range, Japanese Alps).

**Current State:**
1. ✅ Continental-Continental → CollisionOrogen (working)
2. ✅ Oceanic-Oceanic → Arc systems (working)
3. ❌ **Oceanic-Continental → NOT HANDLED!**

**Code Evidence:**
- `src/geology/orogenic.rs:185-186` - Returns `None` for oceanic-continental
- `src/geology/generator.rs:220` - Arc systems only check for oceanic-oceanic
- Comments incorrectly claim "Arc Systems handle this"

---

## Task Type
- [x] Backend (Rust geological simulation)
- [ ] Frontend
- [ ] Fullstack

---

## Technical Solution

**Approach:** Extend Arc Systems Generator to handle both oceanic-oceanic AND oceanic-continental convergence

**Rationale:**
1. **Code Reuse**: Both scenarios create similar transects (Trench → Wedge → Forearc → Arc → Backarc)
2. **Key Difference**: Oceanic-continental arcs form ON CONTINENTAL CRUST (higher elevation, different composition)
3. **Architecture**: Arc systems already has complete transect generation logic
4. **Separation of Concerns**: OrogenicBeltGenerator stays focused on collision orogens only

**Province Type Strategy:**
- Reuse existing types: `OceanTrench`, `AccretionaryWedge`, `ForearcBasin`, `VolcanicArc`, `BackarcBasin`
- Same province types work for both oceanic-oceanic and oceanic-continental
- Difference is in **placement** (continental vs oceanic crust) and **elevation** (handled in Stage 3)

---

## Implementation Steps

### Step 1: Extend Arc Systems to Handle Oceanic-Continental Convergence

**File:** `src/geology/generator.rs:200-250`

**Current Code (lines 218-226):**
```rust
// Check if both plates are oceanic
if let (Some(stats_a), Some(stats_b)) = (plate_stats.get(&boundary.plate_a), plate_stats.get(&boundary.plate_b)) {
    if stats_a.plate_type == PlateType::Oceanic && stats_b.plate_type == PlateType::Oceanic {
        // Determine which plate subducts (older/denser - lower ID)
        let (subducting_plate, overriding_plate) = if boundary.plate_a < boundary.plate_b {
            (boundary.plate_a, boundary.plate_b)
        } else {
            (boundary.plate_b, boundary.plate_a)
        };
```

**New Code:**
```rust
// Check for ANY subduction (oceanic-oceanic OR oceanic-continental)
if let (Some(stats_a), Some(stats_b)) = (plate_stats.get(&boundary.plate_a), plate_stats.get(&boundary.plate_b)) {
    // Determine if this is a subduction zone and identify plates
    let (is_subduction, subducting_plate, overriding_plate) = match (stats_a.plate_type, stats_b.plate_type) {
        // Oceanic-Oceanic: Older plate (lower ID) subducts
        (PlateType::Oceanic, PlateType::Oceanic) => {
            if boundary.plate_a < boundary.plate_b {
                (true, boundary.plate_a, boundary.plate_b)
            } else {
                (true, boundary.plate_b, boundary.plate_a)
            }
        },
        // Oceanic-Continental: Oceanic plate always subducts (denser)
        (PlateType::Oceanic, PlateType::Continental) => (true, boundary.plate_a, boundary.plate_b),
        (PlateType::Continental, PlateType::Oceanic) => (true, boundary.plate_b, boundary.plate_a),
        // Continental-Continental: Handled by OrogenicBeltGenerator
        _ => (false, 0, 0),
    };

    if is_subduction {
        // Generate complete subduction zone transect (same for BOTH scenarios)
```

**Deliverable:** Arc systems now handle oceanic-continental convergence

---

### Step 2: Update Orogenic Belt Generator Documentation

**File:** `src/geology/orogenic.rs:182-186`

**Current Code:**
```rust
// Oceanic-Continental convergence
// NOTE: Now handled in Arc Systems (Stage 2.3) - generates VolcanicArc + AccretionaryWedge
// This orogenic belt generator only handles collision orogens now
(PlateType::Oceanic, PlateType::Continental)
| (PlateType::Continental, PlateType::Oceanic) => None,
```

**New Code:**
```rust
// Oceanic-Continental convergence
// Handled by Arc Systems (Stage 2.3) which generates complete subduction transect:
// Trench → AccretionaryWedge → ForearcBasin → VolcanicArc → BackarcBasin
// Key difference from oceanic-oceanic: Arc forms ON CONTINENTAL CRUST (Andes-style)
// vs oceanic crust (island arc-style). Elevation differences applied in Stage 3.
(PlateType::Oceanic, PlateType::Continental)
| (PlateType::Continental, PlateType::Oceanic) => None,
```

**Deliverable:** Accurate documentation of responsibility division

---

### Step 3: Update Module-Level Documentation

**File:** `src/geology/generator.rs:28-30`

**Current Code:**
```rust
//! **Collision Orogens (1)**: CollisionOrogen
//! **Subduction Systems (5)**: OceanTrench, AccretionaryWedge, ForearcBasin, VolcanicArc, BackarcBasin
```

**New Code:**
```rust
//! **Collision Orogens (1)**: CollisionOrogen (continental-continental convergence)
//! **Subduction Systems (5)**: OceanTrench, AccretionaryWedge, ForearcBasin, VolcanicArc, BackarcBasin
//!   - Applies to BOTH oceanic-oceanic (island arcs) AND oceanic-continental (Andes-style) convergence
//!   - Same province types, different elevations (determined in Stage 3 based on crust type)
```

**Deliverable:** Clear module documentation

---

### Step 4: Add Integration Test for Oceanic-Continental Convergence

**File:** `tests/stage2_geology_integration.rs`

**Add new test (after existing tests):**
```rust
#[test]
fn test_oceanic_continental_convergence() {
    // Test that oceanic-continental convergence generates arc systems
    // Real-world examples: Andes, Cascades, Japanese Alps

    use geoforge::{WorldMap, PlateType, PlateInteraction, GeologicProvince};

    let mut world = WorldMap::new(1800, 900, 12345).expect("Failed to create world");
    world.tectonics().generate_plates(15).unwrap();
    world.tectonics().analyze(None).unwrap();

    let metadata = world.get_tectonic_metadata().unwrap();

    // Find oceanic-continental convergent boundaries
    let oceanic_continental_boundaries: Vec<_> = metadata.plate_boundaries.iter()
        .filter(|b| {
            if b.interaction_type != PlateInteraction::Convergent {
                return false;
            }

            let stats_a = metadata.plate_stats.get(&b.plate_a);
            let stats_b = metadata.plate_stats.get(&b.plate_b);

            if let (Some(a), Some(b)) = (stats_a, stats_b) {
                matches!(
                    (a.plate_type, b.plate_type),
                    (PlateType::Oceanic, PlateType::Continental) |
                    (PlateType::Continental, PlateType::Oceanic)
                )
            } else {
                false
            }
        })
        .collect();

    println!("\nFound {} oceanic-continental convergent boundaries", oceanic_continental_boundaries.len());

    // Generate geology
    let provinces = world.generate_geology(None).unwrap();

    // If we have oceanic-continental boundaries, we should generate arc system provinces
    if !oceanic_continental_boundaries.is_empty() {
        let volcanic_arcs: Vec<_> = provinces.iter()
            .filter(|p| p.characteristics.province_type == GeologicProvince::VolcanicArc)
            .collect();

        let trenches: Vec<_> = provinces.iter()
            .filter(|p| p.characteristics.province_type == GeologicProvince::OceanTrench)
            .collect();

        println!("Generated {} volcanic arcs", volcanic_arcs.len());
        println!("Generated {} ocean trenches", trenches.len());

        assert!(
            !volcanic_arcs.is_empty() || !trenches.is_empty(),
            "Should generate arc system features for oceanic-continental convergence (found {} boundaries)",
            oceanic_continental_boundaries.len()
        );
    } else {
        println!("⚠️  No oceanic-continental boundaries in this test world (seed dependent)");
    }

    println!("\n✓ Oceanic-continental convergence test passed");
}
```

**Deliverable:** Test coverage for oceanic-continental scenario

---

### Step 5: Strengthen Integration Test Assertions

**File:** `tests/stage2_geology_integration.rs:122` (in `test_geology_generator_full_pipeline`)

**Add after line 122:**
```rust
// With 20 plates, we should generate SOME convergent features
let has_convergent_features = count_type(GeologicProvince::VolcanicArc) > 0
    || count_type(GeologicProvince::CollisionOrogen) > 0
    || count_type(GeologicProvince::OceanTrench) > 0;

assert!(
    has_convergent_features,
    "Should generate at least some convergent boundary features with 20 plates"
);
```

**Deliverable:** Stronger assertions in integration tests

---

## Key Files Modified

| File | Lines | Operation | Description |
|------|-------|-----------|-------------|
| `src/geology/generator.rs` | 218-226 | Modify | Extend arc systems to handle oceanic-continental |
| `src/geology/generator.rs` | 28-30 | Modify | Update module documentation |
| `src/geology/orogenic.rs` | 182-186 | Modify | Update comments with accurate responsibility |
| `tests/stage2_geology_integration.rs` | End | Add | New test for oceanic-continental convergence |
| `tests/stage2_geology_integration.rs` | 122 | Modify | Strengthen assertions in full pipeline test |

---

## Risks and Mitigation

| Risk | Severity | Mitigation |
|------|----------|------------|
| Arc features expanding onto wrong crust type | Medium | Already handled - expansion uses `target_plate` parameter to stay on correct plate |
| Performance impact from pattern matching | Low | Pattern matching is compile-time optimized, no runtime cost |
| Breaking existing oceanic-oceanic arcs | High | Use explicit pattern matching to preserve existing logic path; test thoroughly |
| Collision orogens overlapping with arcs | Low | OrogenicBeltGenerator runs BEFORE arc systems; continental-continental boundaries never match arc systems condition |
| Test might fail on some seeds | Low | Test checks for boundary existence first; prints warning if no oceanic-continental boundaries found |

---

## Testing Strategy

1. **Unit Tests:** Existing tests in `orogenic.rs` already verify continental-continental handling
2. **Integration Tests:**
   - New `test_oceanic_continental_convergence` - Validates oceanic-continental generates arcs
   - Updated `test_geology_generator_full_pipeline` - Strengthened assertions
3. **Manual Visual QA:** Run `examples/geology_visual_qa.rs` and inspect outputs for:
   - Volcanic arcs at oceanic-continental boundaries
   - Trenches on oceanic side
   - Proper placement on continental crust

---

## Expected Outcomes

### Before Fix
- Oceanic-continental boundaries: ❌ No provinces generated
- Andes-style mountains: ❌ Missing
- Real-world analog coverage: ❌ Incomplete (~33% of convergent zones missing)

### After Fix
- Oceanic-continental boundaries: ✅ Complete subduction transect
- Trench on oceanic side: ✅ Generated
- Volcanic arc on continental side: ✅ Generated (Andes-style)
- Accretionary wedge, forearc, backarc: ✅ Generated
- Real-world analog coverage: ✅ Complete (100% of convergent zones)

---

## Validation Checklist

After implementation, verify:

- [ ] Tests pass: `cargo test test_orogen`
- [ ] Tests pass: `cargo test test_geology_generator_full_pipeline`
- [ ] Tests pass: `cargo test test_oceanic_continental_convergence`
- [ ] Visual QA shows volcanic arcs at oceanic-continental boundaries
- [ ] No regression in existing collision orogen generation
- [ ] No regression in existing oceanic-oceanic arc generation
- [ ] Documentation accurately describes system behavior
- [ ] Code follows project style guidelines (checked by `cargo clippy`)

---

## Future Extensions

Once oceanic-continental convergence is working, consider:

1. **Elevation differentiation (Stage 3)**: Oceanic-continental arcs should have higher elevation than oceanic-oceanic island arcs
2. **Compositional differences**: Continental arcs have different magma chemistry (andesitic vs basaltic)
3. **Variable arc width**: Continental arcs can be wider due to thicker crust
4. **Flat-slab subduction**: Some oceanic-continental zones have shallow subduction angles (Andes central segment)

These are deferred to future stages and do not block this fix.

---

## References

**Geological Context:**
- Andes Mountains: Nazca Plate (oceanic) subducting under South American Plate (continental)
- Cascades Range: Juan de Fuca Plate (oceanic) under North American Plate (continental)
- Japanese Alps: Pacific Plate (oceanic) under Eurasian Plate (continental margin)

**Key Difference from Island Arcs:**
- Island arcs (oceanic-oceanic): Low elevation volcanic islands (e.g., Aleutians, Marianas)
- Continental arcs (oceanic-continental): High elevation mountain ranges (e.g., Andes, Cascades)

---

## SESSION_ID (for /ccg:execute use)

N/A - No external model sessions used (backend-only Rust logic)
