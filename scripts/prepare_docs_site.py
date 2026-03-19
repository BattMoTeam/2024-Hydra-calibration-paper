from __future__ import annotations

import json
import shutil
from pathlib import Path

from plotly.offline import get_plotlyjs


ROOT = Path(__file__).resolve().parents[1]
DOCS_DIR = ROOT / "docs"
DOCS_ASSETS_DIR = DOCS_DIR / "assets"
DATA_DIR = DOCS_ASSETS_DIR / "data"
IMAGES_DIR = DOCS_ASSETS_DIR / "images"
VENDOR_DIR = DOCS_ASSETS_DIR / "vendor"


def copy_file(source: Path, target: Path) -> str:
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)
    return str(target.relative_to(DOCS_DIR)).replace("\\", "/")


def write_data_script(key: str, payload: object, target: Path) -> str:
    target.parent.mkdir(parents=True, exist_ok=True)
    script = (
        "window.__DOCS_DATA__ = window.__DOCS_DATA__ || {};\n"
        f"window.__DOCS_DATA__[{json.dumps(key)}] = {json.dumps(payload, ensure_ascii=False)};\n"
    )
    target.write_text(script, encoding="utf-8")
    return str(target.relative_to(DOCS_DIR)).replace("\\", "/")


def write_plotly_bundle() -> str:
    target = VENDOR_DIR / "plotly.min.js"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(get_plotlyjs(), encoding="utf-8")
    return str(target.relative_to(DOCS_DIR)).replace("\\", "/")


def load_json(path: Path) -> object:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    IMAGES_DIR.mkdir(parents=True, exist_ok=True)
    VENDOR_DIR.mkdir(parents=True, exist_ok=True)

    plotly_bundle = write_plotly_bundle()

    def prepare_json_data(key: str, source: Path, target_name: str) -> dict[str, str]:
        raw_rel = copy_file(source, DATA_DIR / target_name)
        data = load_json(source)
        script_rel = write_data_script(key, data, DATA_DIR / f"{Path(target_name).stem}.js")
        return {"path": raw_rel, "script": script_rel, "data_key": key}

    def prepare_text_data(key: str, source: Path, target_name: str) -> dict[str, str]:
        raw_rel = copy_file(source, DATA_DIR / target_name)
        text = load_text(source)
        script_rel = write_data_script(key, text, DATA_DIR / f"{Path(target_name).stem}.js")
        return {"path": raw_rel, "script": script_rel, "data_key": key}

    data_map = {
        "validation_reference": prepare_json_data(
            "validation_reference",
            ROOT / "figures" / "battmo-validation-reference.json",
            "battmo-validation-reference.json",
        ),
        "validation_summary": prepare_json_data(
            "validation_summary",
            ROOT
            / "figures"
            / "publication"
            / "INP5-70-120-H0B_graphite-lnmo_schmitt-2026_battmo-vs-experiment-summary.json",
            "battmo-vs-experiment-summary.json",
        ),
        "figure12": prepare_json_data(
            "figure12",
            ROOT / "figures" / "figure-12-cell-balancing-under-equilibrium-assumption.json",
            "figure-12-cell-balancing-under-equilibrium-assumption.json",
        ),
        "figure13": prepare_json_data(
            "figure13",
            ROOT / "figures" / "figure-13-high-rate-calibration-at-2C.json",
            "figure-13-high-rate-calibration-at-2C.json",
        ),
        "figure14": prepare_json_data(
            "figure14",
            ROOT / "figures" / "figure-14-experimental-voltages-and-p2d-results.json",
            "figure-14-experimental-voltages-and-p2d-results.json",
        ),
        "supporting_states": prepare_json_data(
            "supporting_states",
            ROOT / "figures" / "supporting" / "battmo-validation-states.json",
            "battmo-validation-states.json",
        ),
        "rate_study_reference": prepare_json_data(
            "rate_study_reference",
            ROOT / "figures" / "rate-study" / "battmo-rate-study-reference.json",
            "battmo-rate-study-reference.json",
        ),
    }

    publication_gallery = [
        {
            "title": "Figure 12",
            "description": "Cell balancing under equilibrium assumption.",
            "image": copy_file(
                ROOT / "figures" / "figure-12-cell-balancing-under-equilibrium-assumption.png",
                IMAGES_DIR / "publication" / "figure-12-cell-balancing-under-equilibrium-assumption.png",
            ),
        },
        {
            "title": "Figure 13",
            "description": "High-rate calibration at 2C from two different negative-electrode diffusion-coefficient initializations.",
            "image": copy_file(
                ROOT / "figures" / "figure-13-high-rate-calibration-at-2C.png",
                IMAGES_DIR / "publication" / "figure-13-high-rate-calibration-at-2C.png",
            ),
        },
        {
            "title": "Figure 14",
            "description": "Experimental voltages and P2D model results over different discharge rates.",
            "image": copy_file(
                ROOT / "figures" / "figure-14-experimental-voltages-and-p2d-results.png",
                IMAGES_DIR / "publication" / "figure-14-experimental-voltages-and-p2d-results.png",
            ),
        },
        {
            "title": "BattMo vs experiment",
            "description": "Publication-facing validation panel across all five discharge rates.",
            "image": copy_file(
                ROOT
                / "figures"
                / "publication"
                / "INP5-70-120-H0B_graphite-lnmo_schmitt-2026_battmo-vs-experiment.png",
                IMAGES_DIR / "publication" / "battmo-vs-experiment.png",
            ),
        },
    ]

    supporting_gallery = []
    for case_dir in sorted((ROOT / "figures" / "supporting").glob("discharge-rate-*")):
        case_name = case_dir.name.replace("-", " ").title()
        supporting_gallery.append(
            {
                "title": f"{case_name} voltage curve",
                "description": "BattMo voltage curve against the experimental discharge trace.",
                "image": copy_file(
                    case_dir / f"{case_dir.name}-voltage.png",
                    IMAGES_DIR / "supporting" / case_dir.name / f"{case_dir.name}-voltage.png",
                ),
            }
        )
        supporting_gallery.append(
            {
                "title": f"{case_name} state dashboard",
                "description": "Electrolyte, potential, and particle-stoichiometry contour dashboard for the BattMo run.",
                "image": copy_file(
                    case_dir / f"{case_dir.name}-state-dashboard.png",
                    IMAGES_DIR / "supporting" / case_dir.name / f"{case_dir.name}-state-dashboard.png",
                ),
            }
        )

    fair_documents = [
        {
            "label": "Linked-data optimization JSON-LD",
            "description": "Semantic-web description of the optimization and calibration metadata.",
            **prepare_text_data(
                "fair_optimization_jsonld",
                ROOT / "linked-data" / "INP5-70-120-H0B_graphite-lnmo_schmitt-2026_optimization.jsonld",
                "INP5-70-120-H0B_graphite-lnmo_schmitt-2026_optimization.jsonld",
            ),
        },
        {
            "label": "Publication BPX",
            "description": "Publication-facing BPX export of the validation parameter set.",
            **prepare_text_data(
                "fair_bpx",
                ROOT / "parameters" / "INP5-70-120-H0B_graphite-lnmo_schmitt-2026_validation.bpx.json",
                "INP5-70-120-H0B_graphite-lnmo_schmitt-2026_validation.bpx.json",
            ),
        },
        {
            "label": "BattMo base parameters",
            "description": "Canonical BattMo base parameter file used by the model workflow.",
            **prepare_text_data("fair_h0b_base", ROOT / "parameters" / "h0b-base.json", "h0b-base.json"),
        },
        {
            "label": "BattMo equilibrium calibration parameters",
            "description": "Low-rate calibrated BattMo parameters.",
            **prepare_text_data(
                "fair_equilibrium_parameters",
                ROOT / "parameters" / "equilibrium-calibration-parameters.json",
                "equilibrium-calibration-parameters.json",
            ),
        },
        {
            "label": "BattMo high-rate calibration parameters",
            "description": "High-rate calibrated BattMo parameters used in validation.",
            **prepare_text_data(
                "fair_high_rate_parameters",
                ROOT / "parameters" / "high-rate-calibration-parameters.json",
                "high-rate-calibration-parameters.json",
            ),
        },
        {
            "label": "Software citation metadata",
            "description": "Machine-readable citation metadata for the repository.",
            **prepare_text_data("fair_citation_cff", ROOT / "CITATION.cff", "CITATION.cff"),
        },
    ]

    manifest = {
        "summary": {
            "validation_cases": 5,
            "publication_figures": 3,
            "supporting_dashboards": len(supporting_gallery),
            "fair_documents": len(fair_documents),
        },
        "vendor": {"plotly_bundle": plotly_bundle},
        "data": data_map,
        "publication_gallery": publication_gallery,
        "supporting_gallery": supporting_gallery,
        "fair_documents": fair_documents,
    }

    manifest_json = DATA_DIR / "site-manifest.json"
    manifest_json.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    manifest_js = DATA_DIR / "site-manifest.js"
    write_data_script("site_manifest", manifest, manifest_js)
    print(f"Wrote {manifest_json}")
    print(f"Wrote {manifest_js}")


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
