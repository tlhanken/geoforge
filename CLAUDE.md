# Geoforge Project

**Procedural worlds you can edit at any layer**

## Project Overview
Geoforge makes a solar system and generates its focus world layer by layer — tectonics,
geology, elevation, circulation, climate, erosion, hydrology, biomes, resources — with every
layer exportable to images and GIS files, editable or hand-drawn, and re-importable;
downstream layers regenerate from the edits. Derived views (Köppen climates, soils and
deposits, hazards, oceans, pressure) and time (orbits, calendar, eclipses, night sky, storms,
daily weather) are worked out from the stored world. Everything is deterministic by seed.

## Current Status
- Branch `v2-implementation`: the v2 workspace below; v1 has been removed.
- Each layer is checked against real Earth (and Mars, the Moon) in the `known_earth_*` and
  `known_planets` tests; scorecards and open gaps are in `docs/V2_REALISM.md`.
- What comes next: `docs/TODO.md` (the ordered queue) and `docs/V2_ROADMAP.md` (the detail).

## Development Workflow

**IMPORTANT:** When implementing new features, always follow this workflow:

### 1. API Design First
- **Start with interface design** - Define public APIs and function signatures first
- **Review and proposal** - Present the API design for review before implementation
- **Consider ergonomics** - Ensure APIs are intuitive and follow Rust best practices
- **Document early** - Write doc comments during API design phase

### 2. Test-Driven Development
- **Write tests second** - After API approval, write comprehensive test cases
- **Cover edge cases** - Include boundary conditions, error cases, and integration scenarios
- **Test organization:**
  - Unit tests: `#[cfg(test)] mod tests` in source files
  - Integration tests: `/tests/stageN_*_integration.rs` for pipeline stages
  - Examples: `/examples/` for user-facing demonstrations
- **Test before implementation** - Tests should fail initially (Red-Green-Refactor)

### 3. Implementation
- **Implement to tests** - Write code to satisfy the test requirements
- **Maintain test coverage** - Ensure all tests pass before marking complete
- **Refactor iteratively** - Improve code quality while keeping tests green

### 4. Code Review
- **Ensure all code is sensible** - Make sure code is sensible and not bloated.  Implement what is needed efficiently, nothing more.
- **Meets style guidelines** - Ensure all code uses best practices for rust
- **Remove unnessary code** - Simplify down to just needed code, and clean it up

### 5. Workflow Summary
```
1. Design API → Review → Approve
2. Write Tests → Review
3. Implement → Tests Pass → Done
4. Code Review
```

### 6. Test Coverage Standards
- **Unit tests** - All public functions and methods
- **Integration tests** - Full pipeline workflows for each major stage
- **Edge cases** - Boundary conditions, empty inputs, extreme values
- **Error handling** - All error paths and failure modes
- **Determinism** - Verify the same seed reproduces the same output every time.

### 7. Ideas and Proposals Convention

Unscheduled ideas and committed design proposals live in `docs/` under a two-stage lifecycle:

```
docs/ideas/              # Seeds, brainstorming, future possibilities
                         # No template required — just capture the thought
                         # Descriptive kebab-case names (e.g., geologic-cartography-as-narrative-memory.md)

docs/proposals/          # Committed design proposals
                         # Numbered: 0001-name.md, 0002-name.md, ...
                         # Uses docs/proposals/0000-template.md
                         # Status field: draft | accepted | implemented | rejected | deferred
```

Lifecycle:
```
Idea (docs/ideas/)           → write freely, no structure required
  ↓ you decide to pursue it
Proposal (docs/proposals/)   → move file, add structure, assign next number
  ↓ design is approved
Implementation               → reference proposal # in commits, CLAUDE.md, tests
```

- Most ideas will never become proposals — that's fine. They stay in `docs/ideas/` as a record of thinking.
- Proposal numbers are referenced in commit messages (e.g., "Implements proposal 0003") and can be linked from the Pipeline Roadmap.
- When an idea graduates to a proposal, move the file from `docs/ideas/` to `docs/proposals/`, rename with the next number, and expand using the template.

## Pipeline Roadmap

See `docs/V2_ROADMAP.md`. Layers and products: `docs/V2_PIPELINE.md`.

## Development Commands
```bash
just build        # cargo build --release -p geoforge-cli → target/release/geoforge
just test         # cargo test --release --workspace
just check        # clippy + rustfmt --check
just world 42     # make worlds/world-42 and generate it with previews
nix flake check --keep-going   # the CI checks
```
Getting started with the command: `docs/GETTING_STARTED.md`. Formatting follows
`rustfmt.toml` (long lines); run `cargo fmt --all`.

## Architecture
| Crate | Role |
|---|---|
| `geoforge-core` | data and the body's own physics: keys, seeds, geometry, grids, the five shapes, body and system types, products, the `Generator` trait, radiation, magnetism |
| `geoforge-names` | pronounceable names from a seed (no dependencies) |
| `geoforge-cosmos` | galaxies, stellar evolution, systems and their bodies, the stellar neighbourhood |
| `geoforge-surface` | the surface layers as `Generator`s |
| `geoforge-analysis` | derived readings of the surface: Köppen, soils and deposits (and the giant-deposit generator), hazards, oceans, pressure |
| `geoforge-store` | the world directory: manifest, products, edit views and patches, renders, derived cache |
| `geoforge-almanac` | the world in time: ephemeris and calendar, eclipses, sky, storms, weather, hazards near a place |
| `geoforge-pipeline` | the generator registry, plan, provenance, generation, export/import, views, place reports |
| `geoforge-cli` | the `geoforge` command |

Rules: `geoforge-core` has no algorithms that generate and no I/O; generators never branch
on the kind of body (another kind gets another generator); surface layers see a
`BodyContext` and earlier products, never cosmology; algorithms ask a `Grid` for positions
and neighbours; only `geoforge-store` reads or writes files.
