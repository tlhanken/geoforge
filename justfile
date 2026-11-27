# geoforge tasks. `just` lists them; `just <task> <args>` runs one.

world_dir := "worlds"

# List the tasks.
default:
    @just --list

# Build the geoforge command and the viewer (target/release/geoforge, geoforge-viewer).
build:
    cargo build --release -p geoforge-cli -p geoforge-viewer

# Look around a universe in a browser: just view worlds/world-42
view dir port="8137": build
    ./target/release/geoforge-viewer {{dir}} --port {{port}}

# Run every test.
test:
    cargo test --release --workspace

# Lint and check formatting.
check:
    cargo clippy --release --workspace --all-targets
    cargo fmt --all --check

# The Nix checks (each crate's test, clippy, doc and doctest; licences; formatting).
nix-check:
    bloomery check

# Reconcile Cargo.lock and bloomery.lock after changing dependencies or crates.
sync:
    bloomery sync

# Format the code.
fmt:
    cargo fmt --all

# Make a world from a seed and generate every layer with previews: just world 42
world seed="42" *flags="": build
    ./target/release/geoforge new {{world_dir}}/world-{{seed}} --seed {{seed}} {{flags}}
    ./target/release/geoforge generate {{world_dir}}/world-{{seed}} --render

# Describe a world's system and focus: just show worlds/world-42
show dir: build
    ./target/release/geoforge show {{dir}}

# What is at a place: just place worlds/world-42 45 10
place dir lat lon: build
    ./target/release/geoforge place {{dir}} --lat {{lat}} --lon {{lon}}

# Two weeks of weather at a place from a day: just weather worlds/world-42 45 10 100
weather dir lat lon day="0": build
    ./target/release/geoforge weather {{dir}} --lat {{lat}} --lon {{lon}} --day {{day}}

# The weather map of a day: just weather-map worlds/world-42 100
weather-map dir day="0": build
    ./target/release/geoforge weather-map {{dir}} --day {{day}}

# The night sky on a day, with a chart in the world's renders folder: just sky worlds/world-42 30
sky dir day="0": build
    ./target/release/geoforge sky {{dir}} --day {{day}} --to {{dir}}/renders

# A report page with every map and its key, places, weather and sky: just report worlds/world-42
report dir *flags="": build
    ./target/release/geoforge report {{dir}} {{flags}}

# Export a layer to edit: just edit worlds/world-42 rough_elevation
edit dir product: build
    ./target/release/geoforge export {{dir}} {{product}} --caps

# Import an edited layer and regenerate what follows from it: just apply worlds/world-42 rough_elevation
apply dir product: build
    ./target/release/geoforge import {{dir}} {{product}}
    ./target/release/geoforge generate {{dir}} --render
