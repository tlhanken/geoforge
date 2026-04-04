{
  pkgs,
  flake,
  system,
  perSystem,
  ...
}: let
  inherit (flake.lib.${system}) craneLib;
in
  craneLib.devShell {
    packages = [
      pkgs.just
      # perSystem.agenix.default
      # pkgs.omnix

      # Cargo tools
      pkgs.cargo-nextest
      pkgs.cargo-tarpaulin
      perSystem.nixpkgs-unstable.cargo-audit
      perSystem.nixpkgs-unstable.cargo-rail

      # Agent Persona Manager
      # perSystem.agent-persona-manager.persona-cli

      # Nix Tools
      pkgs.nix-fast-build
      pkgs.nix-output-monitor
    ];
  }
