{
  pkgs,
  inputs,
  perSystem,
  ...
}: let
  inherit (inputs.self.lib.${pkgs.stdenv.hostPlatform.system}) craneLib commonArgs;
in
  craneLib.cargoAudit (
    commonArgs
    // {
      name = "workspace-audit";
      advisory-db = inputs.advisory-db;
      nativeBuildInputs = with perSystem.nixpkgs-unstable; [cargo-audit];
    }
  )
