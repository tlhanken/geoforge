# Agents working on geoforge

Read these first, whatever agent you are:

- `openspec/` — specs and changes (OpenSpec, `schema: spec-driven`). A change has a
  `proposal.md`, a `design.md`, `tasks.md` and delta specs under `specs/<capability>/`.
  New features start as a change; tick its tasks as they land.
- `docs/ADDING_A_LAYER.md` — where a new product, derived view or almanac feature must be
  wired (renderer, report, Earth A/B page, tests), with file paths.
- `docs/TODO.md` — the ordered queue; `docs/V2_ROADMAP.md` — the detail and findings;
  `docs/V2_REALISM.md` — the scorecards against Earth.
- `CLAUDE.md` — the development workflow (API first, then tests, then code), commands and
  the crate map; it applies to every agent.

Before finishing: `cargo fmt --all`, `cargo clippy --workspace --all-targets -- -D warnings`,
`cargo test --release --workspace`, and the `known_earth_*` scorecards for layers you
touched.
