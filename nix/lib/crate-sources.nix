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
    geoforge = makeCrate (dir + /geoforge);
  };

  libs = {
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
