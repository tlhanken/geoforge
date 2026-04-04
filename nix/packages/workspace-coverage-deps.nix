{
  pkgs,
  inputs,
  ...
}: let
  inherit (inputs.self.lib.${pkgs.stdenv.hostPlatform.system}) craneLib commonArgs;
in
  craneLib.buildDepsOnly (
    commonArgs
    // {
      pname = "workspace-coverage";
      RUSTFLAGS = "-Cdebuginfo=2 -Cstrip=none --cfg=tarpaulin -Cdebug-assertions=off -Cinstrument-coverage";
      buildPhaseCargoCommand = "cargo test --release --all-features";
    }
  )
