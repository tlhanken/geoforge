{
  pkgs,
  inputs,
  ...
}: let
  versions = import ../lib/versions.nix;
  cargoToml = builtins.fromTOML (builtins.readFile (inputs.self + "/Cargo.toml"));
  cargoMsrv = cargoToml.workspace.package.rust-version;
  cargoVersion = cargoToml.workspace.package.version;
in
  pkgs.runCommand "version-sync-check" {} ''
    # Verify MSRV
    if [ "${versions.msrv}" != "${cargoMsrv}" ]; then
      echo "ERROR: MSRV mismatch!"
      echo "  nix/lib/versions.nix: ${versions.msrv}"
      echo "  Cargo.toml:           ${cargoMsrv}"
      exit 1
    fi

    # Verify Version
    if [ "${versions.this}" != "${cargoVersion}" ]; then
      echo "ERROR: Version mismatch!"
      echo "  nix/lib/versions.nix: ${versions.this}"
      echo "  Cargo.toml:           ${cargoVersion}"
      exit 1
    fi

    echo "Versions are in sync:"
    echo "  MSRV:    ${versions.msrv}"
    echo "  Version: ${versions.this}"
    touch $out
  ''
