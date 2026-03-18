# FAIR Data and Interoperability

This repository includes FAIR-data exports alongside the primary BattMo workflow.

The BattMo JSON files remain the canonical model inputs for the published workflow. The JSON-LD and BPX files are included to support metadata publication, machine readability, and cross-platform interoperability.

## Linked Data

The linked-data exports are:

- [`linked-data/INP5-70-120-H0B_graphite-lnmo_schmitt-2026_optimization.jsonld`](./linked-data/INP5-70-120-H0B_graphite-lnmo_schmitt-2026_optimization.jsonld)

This file expresses the optimization and calibration metadata in JSON-LD using EMMO-oriented terms and units.

## BPX Export

The publication-facing BPX file is:

- [`parameters/INP5-70-120-H0B_graphite-lnmo_schmitt-2026_validation.bpx.json`](./parameters/INP5-70-120-H0B_graphite-lnmo_schmitt-2026_validation.bpx.json)

Create the BPX export:

```powershell
python scripts\export_bpx.py
```

Verify that the exported BPX file runs in PyBaMM:

```powershell
python scripts\verify_bpx.py
```

## Optional BattMo / PyBaMM Comparison

Export the BattMo validation reference curves from MATLAB:

```matlab
startup
run(fullfile('scripts', 'exportValidationReference.m'));
```

Then compare BattMo and PyBaMM:

```powershell
python scripts\compare_battmo_pybamm.py
```

The PyBaMM path is intended as an interoperability check rather than the primary published simulation workflow.

## BPX Approximation Notes

Some BattMo-specific features cannot be represented exactly in standard BPX / PyBaMM:

- the graphite negative-electrode `j0(soc)` table
- independently calibrated volumetric surface area and active-material volume fraction
- cathode OCP boundary behavior outside the first tabulated stoichiometry point

Accordingly:

- the graphite `j0(soc)` data is approximated by a single BPX reaction-rate constant
- BattMo surface-area scaling is folded into exported reaction-rate constants where needed
- the cathode OCP table includes a boundary extrapolation to better align the initialized voltage
