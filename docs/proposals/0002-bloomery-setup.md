# Geoforge on Bloomery — Build Migration

- Proposal #: 0002
- Start Date: 2026-10-11
- Status: draft
- Origin: emerged in discussion (deploy geoforge using bloomery, 2026-10-11)

## Assumptions (unconfirmed — strike in review)

- B1 **Target scope.** "Deploy" means: replace the current flake outputs with a
  `bloomery.mkFlake` workspace and get `nix build` + `nix flake check` green for the
  whole repo. Not CI, not publishing, not releasing.
- B2 **v2 workspace location.** Depends on Proposal 0001 assumption A3 (`v2/` is its own
  Cargo workspace). If v2 is instead folded into the root workspace, §4 changes.
- B3 **Toolchain pinning.** We accept nixpkgs' `rustc` and stop pinning via
  `rust-toolchain.toml`, *unless* review demands otherwise (see Open Questions).

## Summary

Replace Geoforge's hand-rolled `blueprint + crane + rust-overlay` flake with a single
`bloomery.mkFlake` workspace driven by `.bloomery/config.toml`, plus a generated
`bloomery.lock`. Bloomery compiles each crate in the lockfile as an independent Nix
derivation via `rustc`, so leaf edits stop invalidating the dependency tree and the
`test` / `clippy` / `doc` / `doctest` checks run as parallel derivations without Cargo
lock contention.

## Motivation

The current flake expresses the build with crane's `buildDepsOnly` / `cargoBuild` layering,
which compiles the workspace as a Cargo invocation inside a derivation. That has two costs
we are already paying: monolithic invalidation (touching one leaf re-runs a large Cargo
build), and Cargo's target-dir lock serializing the check suite. Bloomery's per-crate
derivation graph removes both. The second motivation is fit: Geoforge's next phase
(Proposal 0001) adds many small crates; a builder that treats each crate as its own
derivation and each check as its own derivation scales with that shape rather than
fighting it.

This is a build-system migration only. No source changes, no behavior changes.

## Design

### 1. Inputs

Remove `crane`, `blueprint`, `rust-overlay`, `advisory-db`, `treefmt-nix`, `systems`,
`nixpkgs-unstable`, `nixpkgs-lib`. Bloomery needs only `nixpkgs` and `bloomery`.

```nix
{
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-25.11";
    # During local iteration, point at the working copy instead of the remote:
    # bloomery.url = "path:/home/tlhanken/workspace/bloomery";
    bloomery.url = "github:actionable-work/bloomery";
  };

  outputs = { nixpkgs, bloomery, ... }:
    bloomery.mkFlake {
      inherit nixpkgs;
      root = ./.;
      systems = [ "x86_64-linux" ];
      overrides = {
        # See §3 — gdal is a sys crate behind the optional export-tiff feature.
        gdal = {
          nativeBuildInputs = [ nixpkgs.legacyPackages.x86_64-linux.pkg-config
                                nixpkgs.legacyPackages.x86_64-linux.clang ];
          buildInputs = [ nixpkgs.legacyPackages.x86_64-linux.gdal ];
        };
      };
    };
}
```

Note: `systems` is passed explicitly because the current flake pins `x86_64-linux` only;
bloomery's default is the three `nix-systems` defaults.

### 2. `.bloomery/config.toml` (workspace root)

Mapped directly from the existing `Cargo.toml` `[profile.release]` and the current flake's
toolchain. Only keys Bloomery recognises are allowed — a typo is an evaluation error, by
design.

```toml
[build]
profileName = "release"
libPackages = false
devPackages = true
unify = false

[toolchain]
linker = "lld"

[profile.release]
optLevel = 3
lto = "fat"          # matches Cargo.toml [profile.release] lto = "fat"
codegenUnits = 1     # matches Cargo.toml codegen-units = 1

[flags]
doc = ["-Dwarnings"]

[devShell]
enable = true

[checks]
enable = true
includePackageChecks = true
noDefaultFeatures = false
workspaceDependencies = false
```

`[checks].noDefaultFeatures` and `workspaceDependencies` default to `true` in Bloomery;
they are set to `false` here because they would hard-fail on the current manifests, which
predate those policies. Tightening them is a follow-up, not part of the migration.

### 3. Sys-crate overrides

The v1 `geoforge` crate has two optional deps behind features:

- `image` — pure Rust, no override needed.
- `gdal` — a sys crate binding to native GDAL. Requires `pkg-config`, `clang` (bindgen),
  and `pkgs.gdal`. Bloomery ships built-in overrides for several sys crates but **not
  gdal**, so it is supplied explicitly in `overrides` (§1).

The `export-png` / `export-tiff` / `export-full` Cargo features are selected per-crate
through an override (`features = [...]`), not through `config.toml` — Bloomery's config
has no global feature table. If the release binary must ship with `export-full`, the
`geoforge` crate gets a colocated `overrides.nix` setting its features.

### 4. v2 workspace composition (depends on B2)

Two supported shapes, both first-class in Bloomery:

- **Single workspace.** v2 crates are added to the root `Cargo.toml` members and the whole
  repo is one Bloomery workspace. Simplest; one `bloomery.lock`.
- **Sub-workspace (`v2/`).** `v2/` is its own Cargo workspace with its own flake, composed
  under the root via Bloomery's sub-flake table:

  ```toml
  [flakes.v2]
  path = "v2"
  checks = "aggregate"
  ```

  `mkFlake` then loads `v2/flake.nix` (itself a `bloomery.mkFlake`) and elevates its
  outputs under `v2:`-prefixed names, with an aggregate check. This is the shape Proposal
  0001's A3 implies, and Bloomery supports it directly — so a separate v2 tree does not
  force a second, disconnected build system.

### 5. Lock and check workflow

1. Populate `config.toml` (§2) and swap `flake.nix` (§1).
2. `nix run github:actionable-work/bloomery#bloomery-wrapped -- sync` (or, using the local
   copy, `nix run ../bloomery#bloomery-wrapped -- sync`) to generate `bloomery.lock` from
   the current `Cargo.lock`. Requires `cargo` on the invoking machine; requires `nix` only
   when a flake lockfile is updated.
3. Commit `Cargo.lock`, `bloomery.lock`, and `flake.lock` (per Bloomery's own guidance).
4. `nix flake check` runs the independent check derivations
   (`<crate>:test|clippy|doc|doctest|bin|lib`, `workspace:lock`).
5. `bloomery check` is the merge gate (per Bloomery's AGENTS.md): static spec validation
   plus selected Nix checks, with retained logs.

### 6. Toolchain pinning (Open Question)

Bloomery's `[toolchain]` table resolves `rustc`, `clippy`, `cargo` as **nixpkgs attribute
paths** (`build-config.nix` → `resolvePackage`). There is no rust-overlay input. So a
specific Rust version is pinned one of two ways:

- **Accept nixpkgs' rustc** (from `nixpkgs-25.11`), and drop the pinned
  `rust-toolchain.toml`. Least moving parts; loses exact-version reproducibility.
- **Pass an overlaid nixpkgs.** `mkFlake` takes `nixpkgs` as an argument, so the flake can
  build `pkgs` with `rust-overlay` applied and hand that in:
  `bloomery.mkFlake { nixpkgs = import nixpkgs { overlays = [ rust-overlay.overlays.default ]; }; ... }`,
  then `[toolchain] rustc = "rust-bin.stable.<ver>.default"`. Keeps pinning; re-introduces
  the rust-overlay input the migration otherwise removes.

This is the one genuine decision the migration surfaces.

### 7. Migration order (low risk first)

1. Land `config.toml` + new `flake.nix` on a scratch branch (`cleanup_and_nix` already
   exists as a candidate) without deleting the old flake.
2. `sync`, then `nix build` and confirm both bins (`geoforge`, `geoforge_v2`) build.
3. `nix flake check`; triage failures (expect the `noDefaultFeatures` /
   `workspaceDependencies` policy checks to be the first casualties — already disabled in
   §2's config).
4. Delete the crane/blueprint inputs and the old flake body in the same commit that makes
   the new one canonical.
5. Do not commit; hand back for review (Bloomery AGENTS.md: commits and pushes are the
   user's).

## Alternatives considered

- **Stay on crane/blueprint.** Viable and already working. Rejected by the user's explicit
  goal of deploying on Bloomery; recorded here so the removed inputs are traceable.
- **Coexist (both flakes).** A parallel `flake.bloomery.nix` would let the two builds be
  compared, but doubles the flake-lock and input surface. Rejected: the migration order in
  §7 already gives a safe rollback (branch + old flake retained until step 4).
- **naersk / buildRustPackage.** Same monolithic-invalidation class as crane; no reason to
  prefer them over what we already have.
- **Pin Rust via rust-overlay (B3).** Kept as an option, not the default (§6): it preserves
  exact version pinning at the cost of one input.

## Open Questions

1. **Toolchain pinning (§6).** nixpkgs `rustc` (drop the pin) vs overlaid nixpkgs (keep the
   pin). Pending.
2. **v2 workspace shape (B2 / Proposal 0001 A3).** Single workspace vs `v2/` sub-workspace
   composed via `[flakes.v2]`. Decides whether §4 is even needed.
3. **Release features.** Does the shipped binary need `export-full` (GDAL), or is GDAL a
   dev-only capability? Decides whether the §3 override is load-bearing or incidental.
4. **`noDefaultFeatures` / `workspaceDependencies`.** When do we tighten these from `false`
   to Bloomery's `true` default?
5. **`crates/libs/Materials/`.** Loose directory with no manifest; clean up or formalize
   with a `Cargo.toml`?
6. **Sub-flake vs monorepo for v2.** If v2 is a sub-flake, does its `bloomery.lock` live in
   `v2/` (a second lock to keep in sync)?