# Adding a layer, a view or an almanac feature

A checklist for anyone adding to geoforge, so a new piece reaches every place it should:
the pipeline, the renders, the report, the Earth A/B page and the tests. Follow the
workflow in `CLAUDE.md` (API first, then tests, then code).

## 1. Decide what it is

| It is… | When | Lives in |
|---|---|---|
| A **product** of a layer | stored, editable, regenerated when its inputs change (heights, winds, rain, biomes) | a `Generator` in `geoforge-surface` |
| A **derived view** | worked out from stored products whenever asked, never edited (Köppen, soils, hazards, sea ice) | `geoforge-pipeline/src/views.rs`, its rules in `geoforge-analysis` |
| An **almanac feature** | depends on the date or a place over time (weather, droughts, sky, tides) | `geoforge-almanac` |

If it changes with the date, it is not a product.

## 2. A new product

1. **Declare it** in `geoforge-core/src/product.rs`: a `ProductId` variant, its entry in
   `ProductId::ALL`, `name()` and `from_name()`; a `ProductData` variant (the shape: `Field<f32>`
   for a quantity, `Field<[f32; 2]>` for a flow, a `Partition` or catalog for classes); a
   typed accessor on `Products`.
2. **Generate it** in the layer's `Generator` (`geoforge-surface`): list it among the layer's
   products; read only what `needs()`/`wants()` declare; respect products already in the
   inputs (authored or patched); put every tunable number in `Params` with an Earth default.
3. **Version it**: bump the generator's `VERSION` whenever its output changes for the same
   inputs (this is how stored worlds know they are stale).
4. **Register a new generator** in `geoforge-pipeline/src/registry.rs` (`Implementation`
   and `dispatch!`). A different kind of body gets a different generator, never a branch.
5. **Render it** in `geoforge-store/src/render/mod.rs`: a ramp or palette and a `legend`,
   with units in the legend.
6. **Report it**: add its name to `GROUPS` in `geoforge-pipeline/src/report.rs`, in
   pipeline order inside its layer's group, and a caption in `caption()`.
7. **Compare it**: add the variant to `Values::of` in `geoforge-pipeline/src/compare.rs`
   (and to `over()` if it belongs to land or sea only), so `geoforge compare` sets it beside
   an observed map.
8. **Document it** in `docs/V2_PIPELINE.md` (the products table) and, if the export or
   import format is new, the editing section.

## 3. A new derived view

1. Add a `View` variant and its `compute` arm in `views.rs`; put the rule itself in
   `geoforge-analysis`, as a function of products, so it can be unit-tested alone.
2. Report it (`GROUPS`, after the layer it is read from) and, when it is a class map,
   compare it (as `compare.rs` does for Köppen).

## 4. A new almanac feature

- Seed every draw from the world (`root_seed.child(b"…")`) and the step or month index,
  never from the query, so the same date gives the same answer however it is asked.
- Time is `geoforge_core::time` (`Instant`, `Span`, a world's `Clock`); the older APIs
  still take Earth days from the present (the focus's northern spring equinox) until the
  `time-and-calendars` change's tasks 3.x land. Show dates through a `Calendar`
  (`geoforge_almanac::calendar::of_world`), and keep statistics over months that follow
  the seasons on twelfths of the year (see the `time` spec).
- Reach it from a command in `geoforge-cli/src/almanac.rs` and from `place` when it is a
  fact about a place.

## 5. Check it against Earth

1. **A baseline.** Find the observed data (`docs/REFERENCE_DATA.md`); convert it once into
   a product of the real-Earth world (`scripts/real-earth/`), so `geoforge compare` shows
   the observed map first.
2. **Known cases.** Add rows to the `known_earth_*` suite for the layer: named places with
   observed values and honest ranges. A case the model misses becomes a gap whose text says
   what was measured and why, never a widened range.
3. **Before and after.** Run the suite and the A/B page before and after; report both, and
   what was lost as well as gained. On the A/B page each layer runs on observed inputs, so
   a miss there belongs to that layer.

## 6. Tests and finish

- Unit tests for the new functions; a determinism test (same seed, same output); an
  integration test that generates a small world and reads the product.
- Physical constants carry units in their names (`_km`, `_m_s`, `_days`, `_ma`) and a doc
  comment with what they mean and where the value comes from.
- `cargo fmt --all`, `cargo clippy --workspace --all-targets -- -D warnings`,
  `cargo test --release --workspace`.
- Update `docs/TODO.md` (done and still needed) and `docs/V2_ROADMAP.md`.
