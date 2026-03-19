from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from publication_names import (
    PUBLICATION_BATTMO_EXPERIMENT_FIGURE_PATH,
    PUBLICATION_BATTMO_EXPERIMENT_SUMMARY_PATH,
)


ROOT = Path(__file__).resolve().parents[1]
FIGURES_DIR = ROOT / "figures"


def load_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as f:
        return json.load(f)


def compute_metrics(
    experimental_time: np.ndarray,
    experimental_voltage: np.ndarray,
    battmo_time: np.ndarray,
    battmo_voltage: np.ndarray,
) -> dict:
    common_end = min(experimental_time[-1], battmo_time[-1])
    mask = experimental_time <= common_end
    experimental_time = experimental_time[mask]
    experimental_voltage = experimental_voltage[mask]
    battmo_interp = np.interp(experimental_time, battmo_time, battmo_voltage)
    diff = battmo_interp - experimental_voltage
    return {
        "common_duration_s": float(common_end),
        "rmse_v": float(np.sqrt(np.mean(diff**2))),
        "mae_v": float(np.mean(np.abs(diff))),
        "max_abs_v": float(np.max(np.abs(diff))),
        "initial_voltage_diff_v": float(diff[0]),
        "final_voltage_diff_v": float(diff[-1]),
    }


def make_plot(cases: list[dict], output_path: Path) -> None:
    fig, axes = plt.subplots(3, 2, figsize=(11, 12), constrained_layout=True)
    axes = axes.ravel()

    for idx, case in enumerate(cases):
        ax = axes[idx]
        current_a = case["current_a"]
        experimental = case["experimental"]
        battmo = case["battmo"]

        experimental_capacity = np.asarray(experimental["time_s"], dtype=float) * current_a / 3600.0
        battmo_capacity = np.asarray(battmo["time_s"], dtype=float) * current_a / 3600.0

        ax.plot(experimental_capacity, experimental["voltage_v"], "k--", linewidth=1.6, label="Experiment")
        ax.plot(battmo_capacity, battmo["voltage_v"], color="#1f77b4", linewidth=2.0, label="BattMo")
        ax.set_title(f"{case['case_name']} ({current_a:.3g} A)")
        ax.set_xlabel("Capacity / Ah")
        ax.set_ylabel("Voltage / V")
        ax.grid(True, alpha=0.25)
        ax.set_ylim(3.0, 4.95)

    axes[-1].axis("off")
    handles, labels = axes[0].get_legend_handles_labels()
    axes[-1].legend(handles, labels, loc="center", frameon=False)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def summarise(reference_path: Path, summary_output: Path, figure_output: Path) -> dict:
    reference = load_json(reference_path)

    cases_out = []
    for case in reference["cases"]:
        experimental_time = np.asarray(case["experimental"]["time_s"], dtype=float)
        experimental_voltage = np.asarray(case["experimental"]["voltage_v"], dtype=float)
        battmo_time = np.asarray(case["battmo"]["time_s"], dtype=float)
        battmo_voltage = np.asarray(case["battmo"]["voltage_v"], dtype=float)

        cases_out.append(
            {
                "case_name": case["case_name"],
                "current_a": float(case["current_a"]),
                "drate": float(case["drate"]),
                "metrics": compute_metrics(
                    experimental_time,
                    experimental_voltage,
                    battmo_time,
                    battmo_voltage,
                ),
                "experimental": case["experimental"],
                "battmo": case["battmo"],
            }
        )

    summary = {
        "reference_file": str(reference_path.relative_to(ROOT)),
        "summary_metrics": [
            {
                "case_name": case["case_name"],
                "current_a": case["current_a"],
                **case["metrics"],
            }
            for case in cases_out
        ],
    }

    summary_output.parent.mkdir(parents=True, exist_ok=True)
    summary_output.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    make_plot(cases_out, figure_output)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot BattMo validation curves against the experimental data.")
    parser.add_argument(
        "--reference",
        type=Path,
        default=FIGURES_DIR / "battmo-validation-reference.json",
        help="BattMo validation reference JSON",
    )
    parser.add_argument(
        "--summary-output",
        type=Path,
        default=PUBLICATION_BATTMO_EXPERIMENT_SUMMARY_PATH,
        help="Summary JSON output",
    )
    parser.add_argument(
        "--figure-output",
        type=Path,
        default=PUBLICATION_BATTMO_EXPERIMENT_FIGURE_PATH,
        help="Figure output",
    )
    args = parser.parse_args()

    summary = summarise(args.reference, args.summary_output, args.figure_output)
    print(json.dumps(summary["summary_metrics"], indent=2))
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
