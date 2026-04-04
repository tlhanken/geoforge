{
  pkgs,
  inputs,
  ...
}: let
  inherit (inputs.self.lib.${pkgs.stdenv.hostPlatform.system}) craneLib commonArgs;
in
  craneLib.buildDepsOnly (
    commonArgs
    // {
      pname = "workspace";
    }
  )
