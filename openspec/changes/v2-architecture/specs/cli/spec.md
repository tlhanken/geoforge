## ADDED Requirements

### Requirement: Command surface
The `geoforge` binary MUST provide the commands `init`, `generate all`, `generate layer`,
`insert`, `export`, `export all`, `list`, `show`, and `seed`.

#### Scenario: every command is invocable
- **WHEN** the binary is run with any listed command and valid arguments
- **THEN** the command executes and exits with code 0

### Requirement: Random seed when omitted
A generating command invoked without `--seed` MUST draw a random seed and print it, so the run
is reproducible afterwards.

#### Scenario: generate without a seed
- **WHEN** `geoforge generate all` is run without `--seed`
- **THEN** a random seed is drawn, printed to the terminal, and used for the run

### Requirement: Missing upstream inputs error
`generate layer <name>` MUST fail, listing the missing inputs, when a declared input layer is
absent from the world. It MUST NOT silently generate upstream layers.

#### Scenario: missing upstream
- **WHEN** `geoforge generate layer climate` runs on a world that has no elevation layer
- **THEN** the command exits non-zero and lists the missing input layers

### Requirement: Overwrite protection
Generation that would replace existing layer data MUST warn and MUST refuse unless `--force`
is given.

#### Scenario: overwrite refused
- **WHEN** a generating command would replace a layer that already exists
- **THEN** it warns and leaves the existing layer unchanged

#### Scenario: forced overwrite
- **WHEN** the same command is run with `--force`
- **THEN** the existing layer is replaced

### Requirement: Layer range selection
`generate all` MUST accept `--from <layer>` and `--to <layer>` to select a contiguous range of
the roster.

#### Scenario: regenerate from a layer
- **WHEN** `geoforge generate all --from climate` is run
- **THEN** climate and every downstream layer are regenerated and earlier layers are left as-is

### Requirement: Machine-readable output
Every command MUST accept `--json` and emit machine-readable output when it is set.

#### Scenario: json output
- **WHEN** any command is run with `--json`
- **THEN** its output is emitted as JSON on stdout

### Requirement: Exit codes
The CLI MUST exit 0 on success, 1 on a generation or validation error, and 2 on a usage error.

#### Scenario: validation error exit
- **WHEN** a command fails validation
- **THEN** the process exits with code 1
