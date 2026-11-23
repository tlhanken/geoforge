# layer-integration

## Why

Each new layer, view or almanac feature has to reach several places — the renderer, the
report, the Earth A/B page, the tests — and a contributor (a person or a smaller model)
working on one piece does not know them. The checklist lives in `docs/ADDING_A_LAYER.md`;
this change makes its rules requirements.

## What Changes

- Requirements for where a new product, view or almanac feature must be wired, and how it
  must be checked against Earth.

## Capabilities

- `layer-integration` — the contract every new layer, view and almanac feature meets

## Impact

Documentation and review only; `docs/ADDING_A_LAYER.md` is the how-to with file paths.
