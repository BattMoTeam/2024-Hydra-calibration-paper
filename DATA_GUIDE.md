# Data Guide

This document summarizes the repository structure, publication naming, and the canonical outputs of the BattMo workflow.

## Repository Structure

- [`raw-data`](./raw-data): experimental cell data used for calibration and validation
- [`parameters`](./parameters): base model parameters and calibrated parameter outputs
- [`linked-data`](./linked-data): JSON-LD representations of the published parameter set
- [`figures/publication`](./figures/publication): primary publication-facing BattMo validation outputs
- [`figures/supporting`](./figures/supporting): per-run voltage curves and BattMo state dashboards
- [`scripts/low-rate-calibration`](./scripts/low-rate-calibration): equilibrium calibration workflow
- [`scripts/high-rate-calibration`](./scripts/high-rate-calibration): high-rate calibration workflow
- [`scripts/runValidation.m`](./scripts/runValidation.m): end-to-end validation against the experimental dataset

## Naming Convention

Publication-facing artifacts use descriptive names built from:

- `INP5-70-120-H0B`: cell model
- `graphite-lnmo`: electrode chemistry
- `schmitt-2026`: paper reference

Examples:

- [`parameters/INP5-70-120-H0B_graphite-lnmo_schmitt-2026_validation.bpx.json`](./parameters/INP5-70-120-H0B_graphite-lnmo_schmitt-2026_validation.bpx.json)
- [`linked-data/INP5-70-120-H0B_graphite-lnmo_schmitt-2026_optimization.jsonld`](./linked-data/INP5-70-120-H0B_graphite-lnmo_schmitt-2026_optimization.jsonld)
- [`figures/publication/INP5-70-120-H0B_graphite-lnmo_schmitt-2026_battmo-vs-experiment.png`](./figures/publication/INP5-70-120-H0B_graphite-lnmo_schmitt-2026_battmo-vs-experiment.png)

Legacy internal BattMo file names such as `h0b-base.json` are kept where they are coupled directly into the original model workflow.

## Canonical BattMo Outputs

Low-rate calibration:

- [`parameters/equilibrium-calibration-parameters.json`](./parameters/equilibrium-calibration-parameters.json)

High-rate calibration:

- [`parameters/high-rate-calibration-parameters.json`](./parameters/high-rate-calibration-parameters.json)

Validation and publication outputs:

- [`figures/battmo-validation-reference.json`](./figures/battmo-validation-reference.json)
- [`figures/figure-12-cell-balancing-under-equilibrium-assumption.fig`](./figures/figure-12-cell-balancing-under-equilibrium-assumption.fig)
- [`figures/figure-12-cell-balancing-under-equilibrium-assumption.png`](./figures/figure-12-cell-balancing-under-equilibrium-assumption.png)
- [`figures/figure-13-high-rate-calibration-at-2C.fig`](./figures/figure-13-high-rate-calibration-at-2C.fig)
- [`figures/figure-13-high-rate-calibration-at-2C.png`](./figures/figure-13-high-rate-calibration-at-2C.png)
- [`figures/figure-14-experimental-voltages-and-p2d-results.fig`](./figures/figure-14-experimental-voltages-and-p2d-results.fig)
- [`figures/figure-14-experimental-voltages-and-p2d-results.png`](./figures/figure-14-experimental-voltages-and-p2d-results.png)
- [`figures/publication/INP5-70-120-H0B_graphite-lnmo_schmitt-2026_battmo-vs-experiment-summary.json`](./figures/publication/INP5-70-120-H0B_graphite-lnmo_schmitt-2026_battmo-vs-experiment-summary.json)
- [`figures/publication/INP5-70-120-H0B_graphite-lnmo_schmitt-2026_battmo-vs-experiment.png`](./figures/publication/INP5-70-120-H0B_graphite-lnmo_schmitt-2026_battmo-vs-experiment.png)
- [`figures/rate-study/battmo-rate-study-reference.json`](./figures/rate-study/battmo-rate-study-reference.json)

## Data Policy

- Treat [`raw-data`](./raw-data) as immutable research input.
- Treat the calibrated parameter files in [`parameters`](./parameters) as canonical derived outputs.
- Treat transient logs, caches, and local debugging artifacts as non-publication material.
