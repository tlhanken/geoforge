{inputs, ...}: let
  eachSystem = inputs.nixpkgs.lib.genAttrs (import inputs.systems);
in
  eachSystem (
    system: let
      pkgs = inputs.nixpkgs.legacyPackages.${system};
      crane = import ./crane.nix {inherit pkgs inputs;};
      versions = import ./versions.nix;
    in {
      inherit (crane) craneLib commonArgs;
      treefmt = inputs.treefmt-nix.lib.evalModule pkgs (import ./treefmt-config.nix);
      package-version = versions.this;
    }
  )
