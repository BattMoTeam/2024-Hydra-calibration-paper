from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pybamm
from bpx import parse_bpx_file

from publication_names import PUBLICATION_BPX_PATH


ROOT = Path(__file__).resolve().parents[1]
PARAMETERS_DIR = ROOT / "parameters"


def run_verification(bpx_path: Path) -> list[dict]:
    bpx = parse_bpx_file(bpx_path)
    parameter_values = pybamm.ParameterValues.create_from_bpx(bpx_path, target_soc=1.0)

    model = pybamm.lithium_ion.DFN({"thermal": "isothermal"})
    solver = pybamm.CasadiSolver(mode="safe")

    lower_cutoff = bpx.parameterisation.cell.lower_voltage_cutoff
    results = []

    for name, validation in bpx.validation.items():
        current = float(np.mean(validation.current))
        duration_hours = float(validation.time[-1] / 3600.0)
        experiment = pybamm.Experiment(
            [f"Discharge at {current:.10g} A for {duration_hours:.10g} hours"],
            termination=f"{lower_cutoff} V",
        )
        simulation = pybamm.Simulation(
            model,
            parameter_values=parameter_values,
            experiment=experiment,
            solver=solver,
        )
        solution = simulation.solve()
        time_s = solution["Time [s]"].entries
        voltage_v = solution["Voltage [V]"].entries
        results.append(
            {
                "case": name,
                "current_a": current,
                "duration_h": duration_hours,
                "points": int(len(time_s)),
                "final_time_s": float(time_s[-1]),
                "initial_voltage_v": float(voltage_v[0]),
                "final_voltage_v": float(voltage_v[-1]),
                "min_voltage_v": float(np.min(voltage_v)),
            }
        )

    return results


def main() -> None:
    parser = argparse.ArgumentParser(description="Verify that a BPX parameter set runs in PyBaMM.")
    parser.add_argument(
        "--input",
        type=Path,
        default=PUBLICATION_BPX_PATH,
        help="BPX file to verify",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional JSON summary file",
    )
    args = parser.parse_args()

    results = run_verification(args.input)
    if args.output is not None:
        args.output.write_text(json.dumps(results, indent=2), encoding="utf-8")
        print(f"Wrote {args.output}")

    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
