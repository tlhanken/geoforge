{
  pkgs,
  inputs,
  ...
}: let
  inherit (inputs.self.lib.${pkgs.stdenv.hostPlatform.system}) craneLib commonArgs package-version;

  cargoArtifacts = inputs.self.packages.${pkgs.stdenv.hostPlatform.system}.workspace-deps;
  pname = "geoforge";
in
  craneLib.buildPackage (commonArgs
    // {
      inherit cargoArtifacts pname;
      version = package-version;
    })
