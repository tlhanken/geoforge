{
  lib,
  root,
}: let
  commonFiles = lib.fileset.unions [
    (root + /Cargo.toml)
    (root + /Cargo.lock)
  ];

  # Helper to create a source for a crate
  makeCrate = path: {
    inherit path;
    fileset = lib.fileset.unions [
      (path + /Cargo.toml)
      (lib.fileset.fileFilter (file: file.hasExt "rs") path)
    ];
  };

  # Define crates
  bins = let
    dir = root + /crates/bins;
  in {
    geoforge-cli = makeCrate (dir + /geoforge-cli);
    # The viewer embeds its page (web/: HTML, its script and three.js) in the binary.
    geoforge-viewer = let
      crate = makeCrate (dir + /geoforge-viewer);
    in
      crate // {fileset = lib.fileset.unions [crate.fileset (dir + /geoforge-viewer/web)];};
  };

  libs = let
    dir = root + /crates/libs;
  in {
    geoforge-core = makeCrate (dir + /geoforge-core);
    geoforge-names = makeCrate (dir + /geoforge-names);
    geoforge-surface = makeCrate (dir + /geoforge-surface);
    geoforge-analysis = makeCrate (dir + /geoforge-analysis);
    geoforge-cosmos = makeCrate (dir + /geoforge-cosmos);
    geoforge-store = makeCrate (dir + /geoforge-store);
    geoforge-almanac = makeCrate (dir + /geoforge-almanac);
    geoforge-pipeline = makeCrate (dir + /geoforge-pipeline);
  };

  allCrates = bins // libs;

  # Combined fileset for the workspace
  workspaceFileset = lib.fileset.unions (
    [commonFiles] ++ (lib.mapAttrsToList (_: crate: crate.fileset) allCrates)
  );
in {
  inherit
    bins
    libs
    allCrates
    workspaceFileset
    ;
}
