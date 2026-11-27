# The viewer embeds its page (web/: HTML, its script and three.js) in the binary.
{lib, ...}: {
  fileset = lib.fileset.unions [
    ./Cargo.toml
    ./src
    ./web
  ];
}
