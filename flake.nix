{
  inputs = {
    # Nixpkgs
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-25.11";
    nixpkgs-unstable.url = "github:NixOS/nixpkgs/nixos-unstable";
    nixpkgs-lib.url = "github:nix-community/nixpkgs.lib";
    systems.url = "github:nix-systems/default";

    # Numtide Utilities
    blueprint = {
      url = "github:numtide/blueprint";
      inputs.nixpkgs.follows = "nixpkgs";
      inputs.systems.follows = "systems";
    };
    treefmt-nix = {
      url = "github:numtide/treefmt-nix";
      inputs.nixpkgs.follows = "nixpkgs";
    };

    # Rust Inputs
    crane.url = "github:ipetkov/crane";
    rust-overlay = {
      url = "github:oxalica/rust-overlay";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    advisory-db = {
      url = "github:rustsec/advisory-db";
      flake = false;
    };

    # Testing Tools
    # agent-persona-manager = {
    #   url = "github:tlhanken/agent-persona-manager";
    #   inputs = {
    #     advisory-db.follows = "advisory-db";
    #     blueprint.follows = "blueprint";
    #     crane.follows = "crane";
    #     nixpkgs.follows = "nixpkgs";
    #     nixpkgs-unstable.follows = "nixpkgs-unstable";
    #     rust-overlay.follows = "rust-overlay";
    #     systems.follows = "systems";
    #     treefmt-nix.follows = "treefmt-nix";
    #   };
    # };
  };

  outputs = inputs: let
    systems = ["x86_64-linux"];
    flake = inputs.blueprint {
      inherit inputs;
      prefix = "nix/";
      inherit systems;
    };
  in
    flake;
}
