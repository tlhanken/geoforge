{
  pkgs,
  inputs,
  ...
}: let
  inherit (inputs.self.lib.${pkgs.stdenv.hostPlatform.system}) craneLib commonArgs;

  cargoArtifacts = inputs.self.packages.${pkgs.stdenv.hostPlatform.system}.workspace-deps;
in
  craneLib.cargoNextest (
    commonArgs
    // {
      inherit cargoArtifacts;
      pname = "workspace";
      doCheck = true;
      partitions = 10;
      partitionType = "hash";
      cargoNextestPartitionsExtraArgs = "--no-tests=pass";
    }
  )
