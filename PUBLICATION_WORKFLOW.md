# Publication Workflow

This repository is organized around a BattMo-first publication workflow.

## Primary Entry Point

From the repository root:

```powershell
.\run-publication.ps1
```

This runs the primary publication workflow in BattMo and generates:

- `figures/battmo-validation-reference.json`
- `figures/figure-12-cell-balancing-under-equilibrium-assumption.fig`
- `figures/figure-12-cell-balancing-under-equilibrium-assumption.png`
- `figures/figure-13-high-rate-calibration-at-2C.fig`
- `figures/figure-13-high-rate-calibration-at-2C.png`
- `figures/figure-14-experimental-voltages-and-p2d-results.fig`
- `figures/figure-14-experimental-voltages-and-p2d-results.png`
- `figures/supporting/...`
- `figures/publication/INP5-70-120-H0B_graphite-lnmo_schmitt-2026_battmo-vs-experiment-summary.json`
- `figures/publication/INP5-70-120-H0B_graphite-lnmo_schmitt-2026_battmo-vs-experiment.png`
- `figures/rate-study/battmo-rate-study-reference.json`

To refresh the interactive GitHub Pages assets after regenerating the figures:

```powershell
python scripts\prepare_docs_site.py
```

To build the site locally:

```powershell
python -m pip install -r requirements-docs.txt
python -m mkdocs serve
```

The primary publication story is:

1. run the calibrated BattMo validation simulations
2. compare BattMo with the experimental discharge data
3. export the paper figures and per-run supporting BattMo dashboards
4. export the BattMo reference curves used elsewhere in the repository

## Optional FAIR / Interoperability Path

To additionally run the BPX / PyBaMM interoperability workflow:

```powershell
.\run-publication.ps1 -IncludeBpx
```

This keeps BattMo as the main research workflow while still producing the FAIR-data interoperability outputs:

- `parameters/INP5-70-120-H0B_graphite-lnmo_schmitt-2026_validation.bpx.json`

## Numerical Resolution Used In The Publication Runner

For the controlled BattMo / PyBaMM rate-study benchmark, both tools use the same requested electrochemical discretisation:

- negative electrode points: `10`
- separator points: `10`
- positive electrode points: `10`
- negative particle points: `10`
- positive particle points: `10`

The same nominal time grid is also requested in both tools:

`numTimesteps = max(240, ceil(8 / DRate))`

Internal nonlinear/adaptive solver substeps remain tool-specific.

## MATLAB-Only Alternative

If you want to stay fully inside MATLAB, run:

```matlab
startup
run(fullfile('scripts', 'exportValidationReference.m'));
run(fullfile('scripts', 'exportRateStudyReference.m'));
run(fullfile('scripts', 'exportPublicationFigures.m'));
```

Then generate the primary validation figure with:

```powershell
python scripts\plot_battmo_validation.py
```
