{
  pkgs,
  inputs,
  ...
}: let
  inherit (inputs.self.lib.${pkgs.stdenv.hostPlatform.system}) craneLib commonArgs;

  cargoArtifacts = inputs.self.packages.${pkgs.stdenv.hostPlatform.system}.workspace-deps;
in
  craneLib.cargoDocTest (
    commonArgs
    // {
      inherit cargoArtifacts;
      doCheck = true;
      pname = "workspace";
    }
  )
