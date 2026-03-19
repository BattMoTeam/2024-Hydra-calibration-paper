from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PARAMETERS_DIR = ROOT / "parameters"
PUBLICATION_FIGURES_DIR = ROOT / "figures" / "publication"

PUBLICATION_PREFIX = "INP5-70-120-H0B_graphite-lnmo_schmitt-2026"

PUBLICATION_BPX_FILENAME = f"{PUBLICATION_PREFIX}_validation.bpx.json"
PUBLICATION_BPX_PATH = PARAMETERS_DIR / PUBLICATION_BPX_FILENAME

PUBLICATION_BATTMO_EXPERIMENT_BASENAME = f"{PUBLICATION_PREFIX}_battmo-vs-experiment"
PUBLICATION_BATTMO_EXPERIMENT_FIGURE_PATH = PUBLICATION_FIGURES_DIR / f"{PUBLICATION_BATTMO_EXPERIMENT_BASENAME}.png"
PUBLICATION_BATTMO_EXPERIMENT_SUMMARY_PATH = PUBLICATION_FIGURES_DIR / f"{PUBLICATION_BATTMO_EXPERIMENT_BASENAME}-summary.json"


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
