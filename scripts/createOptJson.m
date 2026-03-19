%%

clear all
close all

filename      = fullfile(getHydra0Dir(), 'parameters', 'equilibrium-calibration-parameters.json');
jsonstructEC  = parseBattmoJson(filename);
filename      = fullfile(getHydra0Dir(), 'parameters', 'high-rate-calibration-parameters.json');
jsonstructHRC = parseBattmoJson(filename);

input = struct('lowRateParams', jsonstructEC, ...
               'highRateParams', jsonstructHRC, ...
               'include_current_collectors', true, ...
               'useRegionBruggemanCoefficients', true);

output = runHydra(input, 'runSimulation', false);

jsonOpt = output.jsonstruct;

jsonOpt = rmfield(jsonOpt, 'Control');

% Compare
filename = fullfile(getHydra0Dir(), 'parameters', 'h0b-base.json');
jsonBase = parseBattmoJson(filename);
json00 = mergeJsonStructs({jsonstructEC, jsonstructHRC, jsonBase});

jsonDiff(jsonOpt, json00, 'opt', '00');

% An effective P2D model electronic conductivity is used, but we don't
% save it to file
jsonOpt = json00;

% Remove all comment fields from jsonOpt before saving
jsonOpt = removeJsonComments(jsonOpt);

dosave = true;
if dosave
    filename = fullfile(getHydra0Dir(), 'parameters', 'h0b-opt.json');
    writeJsonStruct(jsonOpt, filename);
end

function jsonOut = removeJsonComments(jsonIn)

    jsonOut = jsonIn;
    if isstruct(jsonIn)
        fnames = fieldnames(jsonIn);
        for k = 1:numel(fnames)
            fname = fnames{k};
            if strcmp(fname, 'comment')
                jsonOut = rmfield(jsonOut, fname);
            else
                jsonOut.(fname) = removeJsonComments(jsonIn.(fname));
            end
        end
    elseif iscell(jsonIn)
        for k = 1:numel(jsonIn)
            jsonOut{k} = removeJsonComments(jsonIn{k});
        end
    end

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
