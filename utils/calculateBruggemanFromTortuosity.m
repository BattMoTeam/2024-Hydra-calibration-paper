function json = calculateBruggemanFromTortuosity(model, jsonstructEC, tortuosityRef, bruggemandef)

    am    = 'ActiveMaterial';
    itf   = 'Interface';
    pe    = 'PositiveElectrode';
    ne    = 'NegativeElectrode';
    co    = 'Coating';
    sd    = 'SolidDiffusion';
    ctrl  = 'Control';
    elyte = 'Electrolyte';
    sep   = 'Separator';

    if nargin == 3
        bruggemandef = 'landesfeind';
    else
        bruggemandef = 'newman';
    end

    switch bruggemandef
      case 'landesfeind'
        bruggeman = @(vf, tau) -log(tau)/log(vf);
      otherwise
        bruggeman = @(vf, tau) 1 - log(tau)/log(vf);
    end

    % Get volume fraction from provided jsonstruct or model
    if isempty(jsonstructEC)
        obj = model;
    else
        obj = jsonstructEC;
    end

    ne_bg = bruggeman(obj.(ne).(co).volumeFraction, tortuosityRef.(ne));
    pe_bg = bruggeman(obj.(pe).(co).volumeFraction, tortuosityRef.(pe));
    sep_bg = bruggeman(1 - model.(sep).porosity, tortuosityRef.(sep));

    json = struct(ne, ...
                  struct(co, ...
                         struct('bruggemanCoefficient', ne_bg)), ...
                  pe, ...
                  struct(co, ...
                         struct('bruggemanCoefficient',  pe_bg)), ...
                  sep, ...
                  struct('bruggemanCoefficient', sep_bg));

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
