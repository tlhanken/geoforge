{
  pkgs,
  inputs,
  ...
}: let
  inherit (inputs.self.lib.${pkgs.stdenv.hostPlatform.system}) craneLib commonArgs;

  cargoArtifacts = inputs.self.packages.${pkgs.stdenv.hostPlatform.system}.workspace-coverage-deps;
in
  # TODO: Evaluate using LlvmCov instead of Tarpaulin
  craneLib.cargoTarpaulin (
    commonArgs
    // {
      inherit cargoArtifacts;
      pname = "workspace";
      cargoTarpaulinExtraArgs = "--skip-clean --out lcov --output-dir $out --engine llvm";

      # TODO: Re-enable this when we have a baseline set of tests
      # cargoTarpaulinExtraArgs = "--skip-clean --out lcov --output-dir $out --fail-under 60 --engine llvm";

      # Use to get the rust flags for the dependencies
      # cargoTarpaulinExtraArgs = "--print-rust-flags --skip-clean --out lcov --output-dir $out --fail-under 60 --engine llvm";
    }
  )
