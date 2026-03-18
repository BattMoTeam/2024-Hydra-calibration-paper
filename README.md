# Graphite/LNMO BattMo Dataset and Workflow

This repository accompanies the paper "Comprehensive parameter and electrochemical dataset for a 1 Ah graphite/LNMO battery cell for physical modelling as a blueprint for data reporting in battery research" by Schmitt et al.

The primary scientific workflow is a BattMo workflow for:

- equilibrium calibration
- high-rate calibration
- validation against the experimental discharge dataset
- publication of the calibrated parameter set

## Get Started

1. Clone this repository.
2. Clone BattMo locally.
3. Set `BATTMO_DIR` to the BattMo root containing `startupBattMo.m`.
4. Optionally set `HYDRA_PYTHON_EXECUTABLE` for BattMo JSON validation.
5. Start MATLAB in the repository root and run:

```matlab
startup
```

To reproduce the publication-facing BattMo outputs, run:

```powershell
.\run-publication.ps1
```

## Main Entry Points

- [`run-publication.ps1`](./run-publication.ps1): BattMo-first publication workflow
- [`PUBLICATION_WORKFLOW.md`](./PUBLICATION_WORKFLOW.md): what the publication runner produces
- [`REPRODUCIBILITY.md`](./REPRODUCIBILITY.md): environment and rerun guidance
- [`DATA_GUIDE.md`](./DATA_GUIDE.md): repository structure, naming, and canonical outputs
- [`FAIR_DATA.md`](./FAIR_DATA.md): linked-data JSON-LD and optional BPX / PyBaMM interoperability
- [`docs`](./docs): GitHub Pages-ready interactive site

## Core Contents

- [`raw-data`](./raw-data): experimental input data
- [`parameters`](./parameters): BattMo base and calibrated parameter files
- [`linked-data`](./linked-data): JSON-LD publication metadata
- [`figures`](./figures): publication figures and supporting simulation outputs

## Citation

Citation metadata is in [`CITATION.cff`](./CITATION.cff). The accompanying paper is available at `https://arxiv.org/abs/2601.10507`.

## License

This repository is distributed under the GNU General Public License v3.0 or later. See [`COPYING`](./COPYING).
