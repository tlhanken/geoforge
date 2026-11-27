{
  description = "geoforge: procedural worlds you can edit at any layer";

  inputs = {
    bloomery.url = "github:actionable-work/bloomery";
    nixpkgs.follows = "bloomery/nixpkgs";
  };

  outputs = {
    nixpkgs,
    bloomery,
    ...
  }:
    bloomery.mkFlake {
      inherit nixpkgs;
      root = ./.;
      systems = ["x86_64-linux" "aarch64-linux" "aarch64-darwin"];
    };
}
