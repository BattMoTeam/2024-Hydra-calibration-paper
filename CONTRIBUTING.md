# Contributing

This repository is primarily a research artifact, so contributions should favor traceability, reproducibility, and scientific clarity over convenience.

## Principles

- Preserve the scientific record.
- Keep raw input data immutable.
- Document every behavior-changing modification.
- Prefer reproducible scripts over manual steps.
- Update documentation whenever setup, outputs, or assumptions change.

## Before Opening a Change

1. State the scientific or technical motivation clearly.
2. Identify whether the change affects inputs, calibration, validation, or publication metadata.
3. Re-run the smallest relevant verification workflow.

## Expected Checks

- `startup` completes with BattMo available.
- Any changed calibration or validation script still runs.
- If parameter files change, explain why and record which script produced them.
- If publication-facing metadata changes, update `README.md`, `CITATION.cff`, and any affected guidance documents together.

## Data Handling

- Do not overwrite or silently edit files in `raw-data/`.
- Treat `parameters/equilibrium-calibration-parameters.json` and `parameters/high-rate-calibration-parameters.json` as canonical derived outputs.
- Do not commit transient simulation caches, diary logs, or machine-specific environment files unless they are intentionally part of a documented release artifact.

## Style

- Use explicit paths and deterministic scripts where practical.
- Prefer environment variables over user-specific hard-coded paths.
- Keep MATLAB scripts runnable from a clean checkout after `startup`.
- Add comments only where they clarify non-obvious scientific or technical intent.
