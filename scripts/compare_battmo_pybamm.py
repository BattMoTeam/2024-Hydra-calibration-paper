from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pybamm

from publication_names import PUBLICATION_BPX_PATH


ROOT = Path(__file__).resolve().parents[1]
FIGURES_DIR = ROOT / "figures"
CODEX_FIGURES_DIR = ROOT / "codex" / "figures"
PARAMETERS_DIR = ROOT / "parameters"


def load_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as f:
        return json.load(f)


def run_pybamm_case(parameter_values: pybamm.ParameterValues, current_a: float, duration_h: float) -> tuple[np.ndarray, np.ndarray]:
    model = pybamm.lithium_ion.DFN({"thermal": "isothermal"})
    solver = pybamm.CasadiSolver(mode="safe")
    experiment = pybamm.Experiment([f"Discharge at {current_a:.16g} A for {duration_h:.16g} hours"])
    simulation = pybamm.Simulation(
        model,
        parameter_values=parameter_values,
        experiment=experiment,
        solver=solver,
    )
    solution = simulation.solve()
    time_s = solution["Time [s]"].entries
    voltage_v = solution["Voltage [V]"].entries
    return np.asarray(time_s, dtype=float), np.asarray(voltage_v, dtype=float)


def compute_metrics(reference_time: np.ndarray, reference_voltage: np.ndarray,
                    comparison_time: np.ndarray, comparison_voltage: np.ndarray) -> dict:
    common_end = min(reference_time[-1], comparison_time[-1])
    mask = reference_time <= common_end
    reference_time = reference_time[mask]
    reference_voltage = reference_voltage[mask]
    comparison_interp = np.interp(reference_time, comparison_time, comparison_voltage)
    diff = comparison_interp - reference_voltage
    return {
        "common_duration_s": float(common_end),
        "rmse_v": float(np.sqrt(np.mean(diff**2))),
        "mae_v": float(np.mean(np.abs(diff))),
        "max_abs_v": float(np.max(np.abs(diff))),
        "final_voltage_diff_v": float(comparison_interp[-1] - reference_voltage[-1]),
        "initial_voltage_diff_v": float(comparison_interp[0] - reference_voltage[0]),
    }


def make_plot(cases: list[dict], output_path: Path) -> None:
    fig, axes = plt.subplots(3, 2, figsize=(11, 12), constrained_layout=True)
    axes = axes.ravel()

    for idx, case in enumerate(cases):
        ax = axes[idx]
        current_a = case["current_a"]
        exp = case["experimental"]
        battmo = case["battmo"]
        pybamm_case = case["pybamm"]

        exp_capacity = np.asarray(exp["time_s"], dtype=float) * current_a / 3600.0
        battmo_capacity = np.asarray(battmo["time_s"], dtype=float) * current_a / 3600.0
        pybamm_capacity = np.asarray(pybamm_case["time_s"], dtype=float) * current_a / 3600.0

        ax.plot(exp_capacity, exp["voltage_v"], "k--", linewidth=1.5, label="Experiment")
        ax.plot(battmo_capacity, battmo["voltage_v"], color="#1f77b4", linewidth=1.8, label="BattMo")
        ax.plot(pybamm_capacity, pybamm_case["voltage_v"], color="#d62728", linewidth=1.6, label="PyBaMM BPX")
        ax.set_title(f"{case['case_name']} ({current_a:.3g} A)")
        ax.set_xlabel("Capacity / Ah")
        ax.set_ylabel("Voltage / V")
        ax.grid(True, alpha=0.25)
        ax.set_ylim(3.0, 4.95)

    axes[-1].axis("off")
    handles, labels = axes[0].get_legend_handles_labels()
    axes[-1].legend(handles, labels, loc="center")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def compare(reference_json: Path, bpx_path: Path, summary_path: Path, figure_path: Path) -> list[dict]:
    reference = load_json(reference_json)
    parameter_values = pybamm.ParameterValues.create_from_bpx(bpx_path, target_soc=1.0)

    cases_out = []
    for case in reference["cases"]:
        current_a = float(case["current_a"])
        duration_h = float(case["experimental"]["time_s"][-1]) / 3600.0
        pybamm_time, pybamm_voltage = run_pybamm_case(parameter_values, current_a, duration_h)

        battmo_time = np.asarray(case["battmo"]["time_s"], dtype=float)
        battmo_voltage = np.asarray(case["battmo"]["voltage_v"], dtype=float)
        metrics = compute_metrics(battmo_time, battmo_voltage, pybamm_time, pybamm_voltage)

        case_out = {
            "case_name": case["case_name"],
            "current_a": current_a,
            "drate": float(case["drate"]),
            "metrics": metrics,
            "experimental": case["experimental"],
            "battmo": case["battmo"],
            "pybamm": {
                "time_s": pybamm_time.tolist(),
                "voltage_v": pybamm_voltage.tolist(),
            },
        }
        cases_out.append(case_out)

    summary = {
        "reference_file": str(reference_json.relative_to(ROOT)),
        "bpx_file": str(bpx_path.relative_to(ROOT)),
        "summary_metrics": [
            {
                "case_name": case["case_name"],
                "current_a": case["current_a"],
                **case["metrics"],
            }
            for case in cases_out
        ],
    }

    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    make_plot(cases_out, figure_path)
    return summary["summary_metrics"]


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare BattMo validation curves with PyBaMM BPX simulations.")
    parser.add_argument(
        "--reference",
        type=Path,
        default=FIGURES_DIR / "battmo-validation-reference.json",
        help="BattMo validation export JSON",
    )
    parser.add_argument(
        "--bpx",
        type=Path,
        default=PUBLICATION_BPX_PATH,
        help="BPX parameter file",
    )
    parser.add_argument(
        "--summary-output",
        type=Path,
        default=CODEX_FIGURES_DIR / "battmo-vs-pybamm-bpx-summary.json",
        help="JSON file for comparison metrics",
    )
    parser.add_argument(
        "--figure-output",
        type=Path,
        default=CODEX_FIGURES_DIR / "battmo-vs-pybamm-bpx.png",
        help="Figure file for the curve comparison",
    )
    args = parser.parse_args()

    metrics = compare(args.reference, args.bpx, args.summary_output, args.figure_output)
    print(json.dumps(metrics, indent=2))
    print(f"Wrote {args.summary_output}")
    print(f"Wrote {args.figure_output}")


if __name__ == "__main__":
    main()

"""
  Copyright 2021-2026 SINTEF Industry, Sustainable Energy Technology
  and SINTEF Digital, Mathematics & Cybernetics.

  This file is part of The Battery Modeling Toolbox BattMo

  BattMo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  BattMo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
"""
