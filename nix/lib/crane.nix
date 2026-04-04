{
  pkgs,
  inputs,
}: let
  versions = import ./versions.nix;
  pkgsWithRust = pkgs.extend inputs.rust-overlay.overlays.default;
  rustToolchain = pkgsWithRust.rust-bin.stable.${versions.rust}.default.override {
    extensions = [
      "rust-src"
      "rust-std"
      "llvm-tools-preview"
    ];
  };
  craneLib = (inputs.crane.mkLib pkgs).overrideToolchain rustToolchain;

  crateSources = import ./crate-sources.nix {
    inherit (pkgs) lib;
    root = ../../.;
  };

  workspaceSrc = pkgs.lib.fileset.toSource {
    root = ../../.;
    fileset = crateSources.workspaceFileset;
  };

  commonArgs = {
    src = workspaceSrc;
    strictDeps = true;
    doCheck = false;
    version = "0.0.1";
    buildInputs =
      [
        pkgs.pkg-config
      ]
      ++ pkgs.lib.optionals pkgs.stdenv.isDarwin [
        pkgs.libiconv
      ];
  };
in {
  inherit craneLib commonArgs;
}
