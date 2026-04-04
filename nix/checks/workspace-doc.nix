{
  pkgs,
  inputs,
  ...
}: let
  inherit (inputs.self.lib.${pkgs.stdenv.hostPlatform.system}) craneLib commonArgs;

  cargoArtifacts = inputs.self.packages.${pkgs.stdenv.hostPlatform.system}.workspace-deps;
in
  craneLib.cargoDoc (
    commonArgs
    // {
      inherit cargoArtifacts;
      pname = "workspace";
      env.RUSTDOCFLAGS = "--deny warnings";
    }
  )
