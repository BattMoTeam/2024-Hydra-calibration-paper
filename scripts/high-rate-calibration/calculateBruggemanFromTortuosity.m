function bman = calculateBruggemanFromTortuosity(model, jsonstructEC, tortuosityRef)

    am    = 'ActiveMaterial';
    itf   = 'Interface';
    pe    = 'PositiveElectrode';
    ne    = 'NegativeElectrode';
    co    = 'Coating';
    sd    = 'SolidDiffusion';
    ctrl  = 'Control';
    elyte = 'Electrolyte';
    sep   = 'Separator';

    bruggeman = @(vf, tau) -log(tau)/log(vf);

    bman = struct(pe, bruggeman(jsonstructEC.(pe).(co).volumeFraction, tortuosityRef.(pe)), ...
                  ne, bruggeman(jsonstructEC.(ne).(co).volumeFraction, tortuosityRef.(ne)), ...
                  sep, bruggeman(1 - model.(sep).porosity, tortuosityRef.(sep)));

end


%{
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
%}
