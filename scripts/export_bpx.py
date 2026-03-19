from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

import numpy as np
from bpx import parse_bpx_obj
from scipy.io import loadmat

from publication_names import PUBLICATION_BPX_PATH


ROOT = Path(__file__).resolve().parents[1]
PARAMETERS_DIR = ROOT / "parameters"
RAW_DATA_DIR = ROOT / "raw-data"


def load_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as f:
        return json.load(f)


def merge_dicts(base: dict, update: dict) -> dict:
    merged = dict(base)
    for key, value in update.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = merge_dicts(merged[key], value)
        else:
            merged[key] = value
    return merged


def parse_matlab_table(path: Path) -> np.ndarray:
    text = path.read_text(encoding="utf-8")
    match = re.search(r"=\s*\[(.*?)\];", text, flags=re.DOTALL)
    if match is None:
        raise ValueError(f"Could not find MATLAB table in {path}")
    rows = []
    for line in match.group(1).splitlines():
        line = line.strip()
        if not line:
            continue
        rows.append([float(item) for item in line.split()])
    return np.array(rows, dtype=float)


def extend_table_to_unit_interval(table: np.ndarray) -> np.ndarray:
    x = table[:, 0]
    y = table[:, 1]
    rows = []
    if x[0] > 0.0:
        slope = (y[1] - y[0]) / (x[1] - x[0])
        rows.append([0.0, y[0] + slope * (0.0 - x[0])])
    rows.extend(table.tolist())
    if x[-1] < 1.0:
        slope = (y[-1] - y[-2]) / (x[-1] - x[-2])
        rows.append([1.0, y[-1] + slope * (1.0 - x[-1])])
    return np.array(rows, dtype=float)


def interpolate_table(table: np.ndarray, x_value: float) -> float:
    return float(np.interp(x_value, table[:, 0], table[:, 1]))


def electrolyte_diffusivity_nyman2008(concentration: np.ndarray) -> np.ndarray:
    scaled = concentration / 1000.0
    return 8.794e-11 * scaled**2 - 3.972e-10 * scaled + 4.862e-10


def electrolyte_conductivity_nyman2008(concentration: np.ndarray) -> np.ndarray:
    scaled = concentration / 1000.0
    return 0.1297 * scaled**3 - 2.51 * scaled**1.5 + 3.329 * scaled


def compute_active_fraction_within_solids(
    active_mass_fraction: float,
    active_density: float,
    binder_mass_fraction: float,
    binder_density: float,
    additive_mass_fraction: float,
    additive_density: float,
) -> float:
    active_specific_volume = active_mass_fraction / active_density
    binder_specific_volume = binder_mass_fraction / binder_density
    additive_specific_volume = additive_mass_fraction / additive_density
    solid_specific_volume = active_specific_volume + binder_specific_volume + additive_specific_volume
    return active_specific_volume / solid_specific_volume


def compute_active_material_volume_fraction(
    total_solid_volume_fraction: float,
    active_mass_fraction: float,
    active_density: float,
    binder_mass_fraction: float,
    binder_density: float,
    additive_mass_fraction: float,
    additive_density: float,
) -> float:
    active_fraction_within_solids = compute_active_fraction_within_solids(
        active_mass_fraction,
        active_density,
        binder_mass_fraction,
        binder_density,
        additive_mass_fraction,
        additive_density,
    )
    return total_solid_volume_fraction * active_fraction_within_solids


def compute_porosity(total_solid_volume_fraction: float) -> float:
    return 1.0 - total_solid_volume_fraction


def compute_effective_conductivity(conductivity: float, total_solid_volume_fraction: float, bruggeman: float) -> float:
    return conductivity * total_solid_volume_fraction**bruggeman


def transport_efficiency(porosity: float, bruggeman: float) -> float:
    return porosity**bruggeman


def compute_geometric_surface_area_per_unit_volume(active_volume_fraction: float, particle_radius: float) -> float:
    return 3.0 * active_volume_fraction / particle_radius


def convert_battmo_positive_reaction_rate_constant_to_bpx(
    k_battmo: float, c_max: float, c_e0: float, surface_area_scale: float
) -> float:
    return k_battmo * c_max * math.sqrt(c_e0) * surface_area_scale


def fit_bpx_negative_reaction_rate_constant(
    j0_table: np.ndarray, theta_min: float, theta_max: float, surface_area_scale: float
) -> float:
    faraday = 96485.33212
    soc = j0_table[:, 0]
    j0 = j0_table[:, 1] * 1e4 * surface_area_scale
    sto = theta_min + soc * (theta_max - theta_min)
    mask = (sto > 0.0) & (sto < 1.0)
    basis = faraday * np.sqrt(sto[mask] * (1.0 - sto[mask]))
    return float(np.dot(basis, j0[mask]) / np.dot(basis, basis))


def load_validation_data() -> dict:
    experiment = loadmat(RAW_DATA_DIR / "TE_1473.mat", squeeze_me=True, struct_as_record=False)["experiment"]
    validation = {}
    for idx, (time_h, current, voltage) in enumerate(zip(experiment.time, experiment.current, experiment.voltage), start=1):
        time_s = np.asarray(time_h, dtype=float) * 3600.0
        current_a = -np.asarray(current, dtype=float)
        voltage_v = np.asarray(voltage, dtype=float)
        validation[f"Discharge rate {idx}"] = {
            "Time [s]": time_s.tolist(),
            "Current [A]": current_a.tolist(),
            "Voltage [V]": voltage_v.tolist(),
        }
    return validation


def build_bpx_dict() -> dict:
    base = load_json(PARAMETERS_DIR / "h0b-base.json")
    geom_3d = load_json(PARAMETERS_DIR / "h0b-geometry-3d.json")
    eq = load_json(PARAMETERS_DIR / "equilibrium-calibration-parameters.json")
    high = load_json(PARAMETERS_DIR / "high-rate-calibration-parameters.json")
    params = merge_dicts(base, eq)
    params = merge_dicts(params, high)

    neg = params["NegativeElectrode"]["Coating"]
    pos = params["PositiveElectrode"]["Coating"]
    sep = params["Separator"]
    elyte = params["Electrolyte"]
    ctrl = params["Control"]
    geometry = geom_3d["Geometry"]

    neg_total_solid = neg["volumeFraction"]
    pos_total_solid = pos["volumeFraction"]
    neg_active_share = compute_active_fraction_within_solids(
        neg["ActiveMaterial"]["massFraction"],
        neg["ActiveMaterial"]["density"],
        neg["Binder"]["massFraction"],
        neg["Binder"]["density"],
        neg["ConductingAdditive"]["massFraction"],
        neg["ConductingAdditive"]["density"],
    )
    pos_active_share = compute_active_fraction_within_solids(
        pos["ActiveMaterial"]["massFraction"],
        pos["ActiveMaterial"]["density"],
        pos["Binder"]["massFraction"],
        pos["Binder"]["density"],
        pos["ConductingAdditive"]["massFraction"],
        pos["ConductingAdditive"]["density"],
    )
    neg_active = neg_total_solid * neg_active_share
    pos_active = pos_total_solid * pos_active_share
    neg_porosity = compute_porosity(neg_total_solid)
    pos_porosity = compute_porosity(pos_total_solid)

    neg_transport = transport_efficiency(
        neg_porosity, high["Electrolyte"]["regionBruggemanCoefficients"]["NegativeElectrode"]
    )
    pos_transport = transport_efficiency(
        pos_porosity, high["Electrolyte"]["regionBruggemanCoefficients"]["PositiveElectrode"]
    )
    sep_transport = transport_efficiency(
        sep["porosity"], high["Electrolyte"]["regionBruggemanCoefficients"]["Separator"]
    )

    neg_conductivity = compute_effective_conductivity(
        neg["electronicConductivity"], neg_total_solid, neg["bruggemanCoefficient"]
    )
    pos_conductivity = compute_effective_conductivity(
        pos["electronicConductivity"], pos_total_solid, pos["bruggemanCoefficient"]
    )

    neg_ocp = extend_table_to_unit_interval(parse_matlab_table(PARAMETERS_DIR / "computeOCPanodeH0b.m"))
    pos_ocp = extend_table_to_unit_interval(parse_matlab_table(PARAMETERS_DIR / "computeOCPcathodeH0b.m"))
    neg_j0 = parse_matlab_table(PARAMETERS_DIR / "computeJ0anodeH0b.m")

    neg_radius = neg["ActiveMaterial"]["SolidDiffusion"]["particleRadius"]
    pos_radius = pos["ActiveMaterial"]["SolidDiffusion"]["particleRadius"]
    neg_surface_area_geom = compute_geometric_surface_area_per_unit_volume(neg_active, neg_radius)
    pos_surface_area_geom = compute_geometric_surface_area_per_unit_volume(pos_active, pos_radius)
    neg_surface_area_scale = neg["ActiveMaterial"]["Interface"]["volumetricSurfaceArea"] / neg_surface_area_geom
    pos_surface_area_scale = pos["ActiveMaterial"]["Interface"]["volumetricSurfaceArea"] / pos_surface_area_geom

    neg_theta_min = neg["ActiveMaterial"]["Interface"]["guestStoichiometry0"]
    neg_theta_max = neg["ActiveMaterial"]["Interface"]["guestStoichiometry100"]
    pos_theta_min = pos["ActiveMaterial"]["Interface"]["guestStoichiometry100"]
    pos_theta_max = pos["ActiveMaterial"]["Interface"]["guestStoichiometry0"]
    ocv_0 = interpolate_table(pos_ocp, pos_theta_max) - interpolate_table(neg_ocp, neg_theta_min)
    ocv_100 = interpolate_table(pos_ocp, pos_theta_min) - interpolate_table(neg_ocp, neg_theta_max)

    c_e0 = elyte["species"]["nominalConcentration"]
    neg_k_bpx = fit_bpx_negative_reaction_rate_constant(
        neg_j0, neg_theta_min, neg_theta_max, neg_surface_area_scale
    )
    pos_k_bpx = convert_battmo_positive_reaction_rate_constant_to_bpx(
        pos["ActiveMaterial"]["Interface"]["reactionRateConstant"],
        pos["ActiveMaterial"]["Interface"]["saturationConcentration"],
        c_e0,
        pos_surface_area_scale,
    )

    concentration_grid = np.linspace(0.0, 3000.0, 301)
    validation = load_validation_data()

    electrode_pair_area = geometry["width"] * geometry["length"]
    stack_thickness = geometry["nLayers"] * (
        params["NegativeElectrode"]["CurrentCollector"]["thickness"]
        + neg["thickness"]
        + sep["thickness"]
        + pos["thickness"]
        + params["PositiveElectrode"]["CurrentCollector"]["thickness"]
    )
    external_surface_area = 2.0 * (
        geometry["width"] * geometry["length"]
        + geometry["width"] * stack_thickness
        + geometry["length"] * stack_thickness
    )
    cell_volume = geometry["width"] * geometry["length"] * stack_thickness

    return {
        "Header": {
            "BPX": 1.0,
            "Title": "HYDRA graphite/LNMO validation parameter set",
            "Description": (
                "BPX export of the BattMo parameter set used in runValidation.m. "
                "Negative-electrode kinetics are approximated by a single BPX reaction-rate constant "
                "fitted to the BattMo j0(soc) table because BPX does not support a direct j0(soc) field."
            ),
            "References": "Schmitt et al., arXiv:2601.10507; BattMo parameter export",
            "Model": "DFN",
        },
        "Parameterisation": {
            "Cell": {
                "Electrode area [m2]": electrode_pair_area,
                "External surface area [m2]": external_surface_area,
                "Volume [m3]": cell_volume,
                "Number of electrode pairs connected in parallel to make a cell": geometry["nLayers"],
                "Lower voltage cut-off [V]": ctrl["lowerCutoffVoltage"],
                "Upper voltage cut-off [V]": 4.9,
                "Nominal cell capacity [A.h]": 1.0,
                "Ambient temperature [K]": params["initT"],
                "Initial temperature [K]": params["initT"],
                "Reference temperature [K]": params["initT"],
            },
            "Electrolyte": {
                "Initial concentration [mol.m-3]": c_e0,
                "Cation transference number": elyte["species"]["transferenceNumber"],
                "Diffusivity [m2.s-1]": {
                    "x": concentration_grid.tolist(),
                    "y": electrolyte_diffusivity_nyman2008(concentration_grid).tolist(),
                },
                "Conductivity [S.m-1]": {
                    "x": concentration_grid.tolist(),
                    "y": electrolyte_conductivity_nyman2008(concentration_grid).tolist(),
                },
            },
            "Negative electrode": {
                "Thickness [m]": neg["thickness"],
                "Porosity": neg_porosity,
                "Transport efficiency": neg_transport,
                "Conductivity [S.m-1]": neg_conductivity,
                "Minimum stoichiometry": neg_theta_min,
                "Maximum stoichiometry": neg_theta_max,
                "Maximum concentration [mol.m-3]": neg["ActiveMaterial"]["Interface"]["saturationConcentration"],
                "Particle radius [m]": neg_radius,
                "Surface area per unit volume [m-1]": neg_surface_area_geom,
                "Diffusivity [m2.s-1]": neg["ActiveMaterial"]["SolidDiffusion"]["referenceDiffusionCoefficient"],
                "Diffusivity activation energy [J.mol-1]": neg["ActiveMaterial"]["SolidDiffusion"]["activationEnergyOfDiffusion"],
                "OCP [V]": {"x": neg_ocp[:, 0].tolist(), "y": neg_ocp[:, 1].tolist()},
                "Entropic change coefficient [V.K-1]": 0.0,
                "Reaction rate constant [mol.m-2.s-1]": neg_k_bpx,
            },
            "Positive electrode": {
                "Thickness [m]": pos["thickness"],
                "Porosity": pos_porosity,
                "Transport efficiency": pos_transport,
                "Conductivity [S.m-1]": pos_conductivity,
                "Minimum stoichiometry": pos_theta_min,
                "Maximum stoichiometry": pos_theta_max,
                "Maximum concentration [mol.m-3]": pos["ActiveMaterial"]["Interface"]["saturationConcentration"],
                "Particle radius [m]": pos_radius,
                "Surface area per unit volume [m-1]": pos_surface_area_geom,
                "Diffusivity [m2.s-1]": pos["ActiveMaterial"]["SolidDiffusion"]["referenceDiffusionCoefficient"],
                "Diffusivity activation energy [J.mol-1]": pos["ActiveMaterial"]["SolidDiffusion"]["activationEnergyOfDiffusion"],
                "OCP [V]": {"x": pos_ocp[:, 0].tolist(), "y": pos_ocp[:, 1].tolist()},
                "Entropic change coefficient [V.K-1]": 0.0,
                "Reaction rate constant [mol.m-2.s-1]": pos_k_bpx,
                "Reaction rate constant activation energy [J.mol-1]": pos["ActiveMaterial"]["Interface"]["activationEnergyOfReaction"],
            },
            "Separator": {
                "Thickness [m]": sep["thickness"],
                "Porosity": sep["porosity"],
                "Transport efficiency": sep_transport,
            },
            "User-defined": {
                "Negative electrode total solid volume fraction": neg_total_solid,
                "Positive electrode total solid volume fraction": pos_total_solid,
                "Negative electrode active material share within solids": neg_active_share,
                "Positive electrode active material share within solids": pos_active_share,
                "Negative electrode active material volume fraction": neg_active,
                "Positive electrode active material volume fraction": pos_active,
                "Negative electrode inactive solid volume fraction": neg_total_solid - neg_active,
                "Positive electrode inactive solid volume fraction": pos_total_solid - pos_active,
                "Open-circuit voltage at 0% SOC [V]": ocv_0,
                "Open-circuit voltage at 100% SOC [V]": ocv_100,
                "BattMo negative electrode volumetric surface area [m-1]": neg["ActiveMaterial"]["Interface"]["volumetricSurfaceArea"],
                "BattMo positive electrode volumetric surface area [m-1]": pos["ActiveMaterial"]["Interface"]["volumetricSurfaceArea"],
                "BPX negative electrode surface area per unit volume [m-1]": neg_surface_area_geom,
                "BPX positive electrode surface area per unit volume [m-1]": pos_surface_area_geom,
                "BattMo negative electrode exchange-current density j0 [A.m-2]": {
                    "x": neg_j0[:, 0].tolist(),
                    "y": (neg_j0[:, 1] * 1e4).tolist(),
                },
                "BattMo negative electrode BPX-fitted reaction rate constant [mol.m-2.s-1]": neg_k_bpx,
            },
        },
        "Validation": validation,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Export the validated BattMo parameter set to BPX.")
    parser.add_argument(
        "--output",
        type=Path,
        default=PUBLICATION_BPX_PATH,
        help="Output BPX file path",
    )
    args = parser.parse_args()

    bpx_dict = build_bpx_dict()
    bpx_model = parse_bpx_obj(bpx_dict).model_dump(by_alias=True, exclude_none=True)
    serialised = json.dumps(bpx_model, indent=2)
    args.output.write_text(serialised, encoding="utf-8")
    print(f"Wrote {args.output}")


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
