{
  pkgs,
  inputs,
  perSystem,
  ...
}: let
  inherit (inputs.self.lib.${pkgs.stdenv.hostPlatform.system}) craneLib commonArgs;

  cargoArtifacts = inputs.self.packages.${pkgs.stdenv.hostPlatform.system}.workspace-deps;
in
  craneLib.mkCargoDerivation (
    commonArgs
    // {
      inherit cargoArtifacts;
      pname = "workspace-unify";
      doInstallCargoArtifacts = false;
      nativeBuildInputs = with pkgs; [
        perSystem.nixpkgs-unstable.cargo-rail
        gitMinimal
      ];
      buildPhaseCargoCommand = ''
        git init > /dev/null
        cargo rail unify --check
      '';
    }
  )
