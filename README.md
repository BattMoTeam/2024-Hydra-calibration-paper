# 2024-Hydra-calibration-paper

This repository contains scripts accompying the paper.

To run the scripts, first install BattMo according to the instructions
on <https://github.com/BattMoTeam/BattMo>. Then modify the `startup.m`
file in this repository to point to the directory where you installed
BattMo.

There are three main scripts:

- `runEquilibriumCalibration.m` for performing the calibration under equilibrium assumption. This script saves the calibrated parameters to disk as a json file for inspection and later use.

- `runHighRateCalibration.m` performs the high-rate calibration. This will take a few minutes to run. It requires the parameters from `runEquilibriumCalibration`. The script also saves the calibrated parameters to disk in json format.

- `runValidation.m` uses both sets of calibrated parameters, and runs and compares the P2D model with the experimental results at the provided rates.

## Linked Data Implementation

The goal of the linked data implementation is to express the parameters in a way that:

- enables Semantic Web integration
- facilitates cross-platform interoperability
- adheres to the FAIR data principles