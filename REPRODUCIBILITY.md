# Reproducibility Guide

This document describes the expected runtime environment and the workflow needed to reproduce the published results.

## Runtime Environment

- MATLAB with Python support compatible with your selected interpreter
- Python for BattMo JSON validation and the optional BPX / PyBaMM workflow
- BattMo clone containing `startupBattMo.m`

## Required External Dependency

BattMo is not vendored in this repository. You must provide a local BattMo clone containing `startupBattMo.m`.

Recommended setup:

```powershell
$env:BATTMO_DIR = '/path/to/BattMo'
$env:HYDRA_PYTHON_EXECUTABLE = '/path/to/python'
```

## Bootstrap

From the repository root:

```matlab
startup
```

`startup.m` adds repository paths, initializes BattMo, loads the MRST modules required by the scripts, and tries to configure Python automatically.

Environment variables take priority over path guessing:

- `BATTMO_DIR`
- `HYDRA_PYTHON_EXECUTABLE`
- `PYTHON_EXECUTABLE`
- `BATTMO_PYTHON_EXECUTABLE`

## Preferred Publication Runner

The preferred publication entry point is the BattMo-first runner:

```powershell
.\run-publication.ps1
```

This executes the publication-facing BattMo workflow and generates the primary BattMo-versus-experiment outputs.
It also exports BattMo-side `.fig` and `.png` files for paper figures 12, 13, and 14 together with per-run supporting voltage curves and state dashboards in `figures/supporting/`.

To additionally run the optional BPX / PyBaMM interoperability workflow:

```powershell
.\run-publication.ps1 -IncludeBpx
```

## Full Workflow

Run the scripts in this order:

```matlab
run(fullfile('scripts', 'low-rate-calibration', 'runEquilibriumCalibration.m'));
run(fullfile('scripts', 'high-rate-calibration', 'runHighRateCalibration.m'));
run(fullfile('scripts', 'runValidation.m'));
```

Expected canonical outputs:

- `parameters/equilibrium-calibration-parameters.json`
- `parameters/high-rate-calibration-parameters.json`

Additional generated artifacts can include:

- packed simulation output in `output/`
- MATLAB diary logs beginning with `_diary-`
- optional figures exported from the scripts

## Optional BPX / PyBaMM Workflow

Python dependencies are listed in `requirements.txt`.

Create the BPX export:

```powershell
python scripts\export_bpx.py
```

Verify that the BPX file runs in PyBaMM:

```powershell
python scripts\verify_bpx.py
```

Export the BattMo validation reference curves from MATLAB:

```matlab
startup
run(fullfile('scripts', 'exportValidationReference.m'));
```

Compare BattMo and PyBaMM:

```powershell
python scripts\compare_battmo_pybamm.py
```

This workflow generates:

- `parameters/INP5-70-120-H0B_graphite-lnmo_schmitt-2026_validation.bpx.json`
- `figures/battmo-validation-reference.json`

The exported BPX file is based on the exact BattMo validation inputs:

- `parameters/h0b-base.json`
- `parameters/equilibrium-calibration-parameters.json`
- `parameters/high-rate-calibration-parameters.json`

Known limitation:

- BattMo's graphite negative-electrode `j0(soc)` table is approximated by a single BPX reaction-rate constant because standard BPX does not provide an exact equivalent field.
- BattMo's calibrated interface surface areas are exported to BPX through geometry-consistent surface-area values with compensating reaction-rate scaling because PyBaMM assumes `a = 3 * eps_s / r`.
- The exported cathode OCP table includes a linearly extrapolated boundary point below the first raw tabulated stoichiometry to reproduce BattMo's evaluated 100% SOC voltage more closely.

## Headless Verification Command

For headless verification, run:

```powershell
matlab -batch "startup; run(fullfile('scripts','runValidation.m'));"
```

## Troubleshooting

### BattMo is not found

Set `BATTMO_DIR` explicitly to the BattMo root that contains `startupBattMo.m`.

### Python is not configured

If BattMo JSON validation fails, set Python explicitly before running the workflow:

```matlab
pyenv('Version', '/path/to/python');
startup
```

### MATLAB rejects the installed Python version

Python support depends on the MATLAB release. If you are using an older MATLAB version, install a Python version supported by that MATLAB release or upgrade MATLAB.

MathWorks compatibility table:

- https://www.mathworks.com/support/requirements/python-compatibility.html

## What To Archive For Publication

For a publication snapshot, archive at minimum:

- all source code in this repository
- the raw input data in `raw-data/`
- the canonical calibrated parameter files in `parameters/`
- the linked-data exports in `linked-data/`
- the exact software citation metadata in `CITATION.cff`

Transient cache directories and local logs do not need to be archived as part of the canonical source release.
