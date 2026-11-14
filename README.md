# h2-data-gen

Synthetic data generation for PHREEQC models of hydrogen storage in rock.
`datagen.py` samples the uncertain inputs of a model (temperature, porosity,
water saturation, brine composition, gas composition, mineral fractions),
runs PHREEQC once per sample and collects the results in one folder per run.

## Requirements

- Python 3.9 or newer with `numpy`, `scipy` and `matplotlib`, plus `tzdata`
  on Windows (run folders are named in Istanbul time, and Windows Python has
  no time-zone data of its own):

  ```bash
  pip install numpy scipy matplotlib tzdata
  ```

- PHREEQC 3 (the command-line program) and its `phreeqc.dat` database.
  Set their paths in each model's `config.cfg` (see [`[phreeqc]`](#phreeqc)).

## Usage

```bash
python datagen.py
```

This runs every model in `models/` whose config has `enabled = true`, one
model after another. At the end it prints one line per model: `done`,
`smoke test failed` or `disabled`.

To run only some models, set `enabled = false` in the others.

## Repository layout

```
datagen.py                  the script
models/
    <name>/
        <name>.phr          PHREEQC template with {PLACEHOLDERS}
        config.cfg          settings and sampled parameters for this model
        runs/               created by datagen.py (not in git)
```

A model's folder name and its template name must match: `models/trona/`
holds `trona.phr`.

Current models: `calcite`, `dolomite`, `mica-carbonate-schist`, `mica-schist`,
`quartzite-schist`, `sandstone`, `trona`.

## What a run does

For each enabled model, `datagen.py`:

1. Reads `config.cfg` and the `.phr` template.
2. Draws all `n_simulations` samples up front from one seed.
3. Creates `models/<name>/runs/<YYYYMMDD-HHMMSS>/` (Istanbul time) and copies the config and
   template into it, so the run records exactly what produced it.
4. Writes one input file per sample to `input/`, before anything runs.
5. Runs the smoke test, if `smoke_test = true` (see below).
6. Runs the simulations in parallel on `workers` processes. For each one it
   fills the template's placeholders, runs PHREEQC, checks the result and
   saves it.
7. Writes `report.txt`.

### Run folder

```
models/<name>/runs/<timestamp>/
    input/00001_input.txt     the placeholder values of each sample
    output/00001_output.txt   PHREEQC's selected output of each usable run
    failed/00001_input.txt    input files of failed, timed-out and flagged runs
    report.txt                settings, seed, counts, parameter ranges, failures
    config.cfg                copy of the config used
    <name>.phr                copy of the template used
    run.log, logs/, plot/     only with debug = true
```

- **Input files** have one `{PLACEHOLDER} value` line per entry of
  `input_file` in `[output]`. The values are the exact text PHREEQC was given.
- **Output files** are tab-separated, with a header line and one row per
  kinetic step. The columns are the template's `USER_PUNCH -headings`, minus
  any in `drop_columns`.
- **`logs/`** (debug only) holds the filled-in PHREEQC input (`.pqi`) and
  PHREEQC's output (`.pqo`) of every simulation, plus a traceback when the
  script itself fails on a run.
- **`plot/`** (debug only) holds one figure per simulation, laid out by the
  `[plot.*]` sections.

### Result status

Every simulation ends in one of four states:

| Status | Meaning | Where it goes |
|---|---|---|
| `success` | Ran to the end and passed every check | `output/` |
| `flagged` | Ran to the end but failed a check in `[checks]` | `failed/` (or `output/` with `include_flagged = true`) |
| `failed` | PHREEQC reported an error, exited abnormally or stopped early | `failed/` |
| `timeout` | Took longer than `timeout_s` | `failed/` |

`report.txt` lists every run that isn't `success`, with the reason.

### Smoke test

With `smoke_test = true`, the script first runs one simulation with every
sampled parameter at the median of its distribution. If that run isn't a
`success`, the model is skipped: its run folder keeps the input files and a
report saying why, and the script moves on to the next model. If it passes,
the log shows its runtime and a rough estimate for the whole batch.

## The config file

Each model has a `config.cfg` in INI format. Lists are comma-separated and
may continue on indented lines. Comments start with `#`; on the same line as
a value, put a space before the `#`.

### `[run]`

| Key | Meaning |
|---|---|
| `enabled` | `false` skips this model. |
| `model_encoding` | Encoding of the `.phr` file, e.g. `utf-8`, `cp1252`. |
| `n_simulations` | Number of samples to draw and run. |
| `debug` | `true` also writes `run.log`, `logs/` and `plot/`. |
| `workers` | Number of PHREEQC runs in parallel. |
| `timeout_s` | Time limit per simulation, in seconds. |
| `smoke_test` | Run the median sample once before the batch. |

### `[phreeqc]`

| Key | Meaning |
|---|---|
| `executable` | PHREEQC program: a name on `PATH` or a full path. |
| `database` | Full path to the database, normally `phreeqc.dat`. |

### `[sampling]`

| Key | Meaning |
|---|---|
| `method` | `lhs` (Latin hypercube), `sobol` (use a power of 2 for `n_simulations`) or `random`. |
| `seed` | `auto` or an integer. The seed used is written to `report.txt`, so `auto` runs can be repeated. |

### `[model]`

| Key | Meaning |
|---|---|
| `volume_m3` | Rock volume represented by one simulation. |
| `total_pressure_atm` | Total gas pressure. Must match `-pressure` in the template. |
| `h2_pressure_atm` | H2 partial pressure. Must match `Hdg(g)` in the template. |
| `normalize_mineral_fractions` | `true` scales the mineral fractions of each sample to sum to 1. |

### `[parameter.<name>]`

One section per input parameter. Sampled parameters are drawn in the order
of their sections. Each section has a `distribution` and the keys that
distribution needs:

| `distribution` | Keys | Value |
|---|---|---|
| `uniform` | `low`, `high` | Uniform between the bounds. |
| `loguniform` | `low`, `high` | Uniform in log10 between the bounds. |
| `triangular` | `low`, `mode`, `high` | Triangular. |
| `normal` | `low`, `high`, `mean`, `std` | Normal, truncated to the bounds. |
| `lognormal` | `low`, `high`, `median`, `sigma_log10` | Lognormal, truncated to the bounds. `sigma_log10` is the standard deviation of log10(x). |
| `constant` | `value` | Not sampled: always `value`. |
| `proportional` | `to`, `base`, `base_to` | Not sampled: `<to> * base / base_to`, keeping the base model's ratio to another parameter (e.g. Cl to Na). |

Add `placeholder = {NAME}` to put the value straight into the template. For
dissolved species the values are in mol/kgw, and for gases they are partial
pressures in atm.

Every model needs these parameters, because the script derives the volumes,
CH4 and the mineral amounts from them:

| Parameter | Used for |
|---|---|
| `porosity` | Pore volume, and the solid volume for the minerals. |
| `water_fraction` | Share of the pore volume filled with water. |
| `pCO2`, `pH2S`, `pN2` | Gas partial pressures; CH4 makes up the rest. |
| one fraction parameter per mineral | Named by `fraction_parameter` in `[mineral.*]`. |

### `[mineral.<name>]`

One section per kinetic mineral in the template.

| Key | Meaning |
|---|---|
| `fraction_parameter` | The `[parameter.*]` holding this mineral's volume fraction of the solid. |
| `moles_placeholder`, `area_placeholder` | Template placeholders for its initial moles and reactive surface area. |
| `diameter_m` | Grain diameter (m). |
| `density_kg_m3` | Density (kg/m³). |
| `molar_mass_kg_mol` | Molar mass (kg/mol). |
| `reactive_fraction` | Share of the geometric surface that reacts. Also applied to the moles. |
| `moles_multiplier`, `area_multiplier` | Extra correction factors; 1 for none. |

### Derived values

These placeholders are filled from the parameters above, not from a
`[parameter.*]` section:

| Placeholder | Value |
|---|---|
| `{WATER_VOLUME}` | `volume_m3 × porosity × 1000 × water_fraction` (L) |
| `{GAS_VOLUME}` | `volume_m3 × porosity × 1000 − {WATER_VOLUME}` (L) |
| `{CH4}` | `total_pressure_atm − pCO2 − pH2S − pN2 − h2_pressure_atm` (atm) |
| mineral moles | `mass / molar_mass × reactive_fraction × moles_multiplier` |
| mineral area | `6 / diameter × mass / density × reactive_fraction × area_multiplier` (m²) |

Here `mass = density × volume_m3 × (1 − porosity) × fraction`. All values are
written in scientific notation with four decimals (e.g. `1.2345e-03`).

### `[checks]`

Checks on each completed run. A run that fails one is `flagged`.

| Key | Meaning |
|---|---|
| `time_column` | Output column with the time in days. |
| `end_time_tolerance` | The run must reach this share of the `-steps` time in the template, or it counts as `failed`. |
| `ph_column`, `ph_range` | pH column and its allowed range. |
| `nonnegative` | Columns that must never go below 0. |
| `must_change` | Columns that must differ between the first kinetic step and the end. Catches reactions that never run. Leave empty to turn this check off. |
| `change_atol` | Smallest change that counts for `must_change`. |

### `[output]`

| Key | Meaning |
|---|---|
| `input_file` | Placeholders written to each input file, in this order. Must cover every placeholder in the template. |
| `drop_columns` | Output columns left out of the output files. |
| `include_flagged` | `true` saves flagged runs in `input/` and `output/` instead of `failed/`. |
| `keep_initial_row` | `true` also keeps PHREEQC's t = 0 row. That row is the initial-solution calculation, done before the gas phase and kinetics exist, so it is dropped by default. |

### `[plot.<title>]`

One section per panel of the debug figure, in order, two panels per row.

| Key | Meaning |
|---|---|
| `columns` | Output columns to plot against time. |
| `scale` | Factor on the values, e.g. `1000` for mol to mmol. |
| `ylabel` | Y-axis label. |

## Adding a model

1. Create `models/<name>/` and save the template as `<name>.phr`. Replace each
   value to sample with a `{PLACEHOLDER}`, and use `{WATER_VOLUME}`,
   `{GAS_VOLUME}` and `{CH4}` for the derived values.
2. Copy the `config.cfg` of the most similar model and edit it: the
   parameters and their ranges, the minerals, `input_file`, and the checks
   and plots for the template's output columns.
3. Set `model_encoding` to the template's encoding. The existing templates
   are UTF-8.
4. Try it with a small `n_simulations` and `debug = true` before a full run.
