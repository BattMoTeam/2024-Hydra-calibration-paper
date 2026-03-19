% Export BattMo discharge curves for a controlled low-rate study.

clear all
close all

mrstDebug(0);

ctrl = 'Control';

getTime = @(states) cellfun(@(s) s.time, states);
getE = @(states) cellfun(@(s) s.(ctrl).E, states);

rates = [1/30, 1/10, 1/5, 1/3];
rateLabels = {'C/30', 'C/10', 'C/5', 'C/3'};

filename      = fullfile(getHydra0Dir(), 'parameters', 'equilibrium-calibration-parameters.json');
jsonstructEC  = parseBattmoJson(filename);
filename      = fullfile(getHydra0Dir(), 'parameters', 'high-rate-calibration-parameters.json');
jsonstructHRC = parseBattmoJson(filename);

inputCap  = struct('lowRateParams', jsonstructEC, ...
                   'include_current_collectors', true);
outputCap = runHydra(inputCap, 'runSimulation', false);
cap       = computeCellCapacity(outputCap.model);
cutoffV   = outputCap.model.(ctrl).lowerCutoffVoltage;

cases = cell(numel(rates), 1);

for k = 1:numel(rates)

    DRate = rates(k);
    totalTime = 1.2 * hour / DRate;
    numTimesteps = max(240, ceil(8 / DRate));

    input = struct('DRate'                         , DRate            , ...
                   'totalTime'                     , totalTime        , ...
                   'numTimesteps'                  , numTimesteps     , ...
                   'lowRateParams'                 , jsonstructEC     , ...
                   'highRateParams'                , jsonstructHRC    , ...
                   'useRegionBruggemanCoefficients', true             , ...
                   'include_current_collectors'    , true);

    output = runHydra(input, 'clearSimulation', true);

    simTime = getTime(output.states);
    simVoltage = getE(output.states);

    cutoffIndex = find(simVoltage <= cutoffV + 1e-6, 1, 'first');
    if isempty(cutoffIndex)
        cutoffIndex = numel(simTime);
    end

    simTime = simTime(1:cutoffIndex);
    simVoltage = simVoltage(1:cutoffIndex);

    currentA = DRate * cap / hour;

    cases{k} = struct('case_name'     , sprintf('Discharge %s', rateLabels{k}) , ...
                      'rate_label'    , rateLabels{k}                           , ...
                      'rate_c'        , DRate                                   , ...
                      'current_a'     , currentA                                , ...
                      'time_s'        , simTime(:)'                             , ...
                      'voltage_v'     , simVoltage(:)'                          , ...
                      'capacity_ah'   , simTime(:)' * currentA / hour          , ...
                      'final_voltage' , simVoltage(end)                         , ...
                      'cutoff_voltage', cutoffV);

end

exportStruct = struct('generated_at'    , datestr(now, 'yyyy-mm-ddTHH:MM:SS')                     , ...
                      'source'          , 'scripts/exportRateStudyReference.m'                    , ...
                      'matlab_release'  , version('-release')                                     , ...
                      'capacity_as'     , cap                                                      , ...
                      'capacity_ah'     , cap / hour                                               , ...
                      'parameter_files' , {{'parameters/h0b-base.json', ...
                                             'parameters/equilibrium-calibration-parameters.json', ...
                                             'parameters/high-rate-calibration-parameters.json'}} , ...
                      'cases'           , {cases});

outdir = fullfile(getHydra0Dir(), 'figures', 'rate-study');
if ~isfolder(outdir)
    mkdir(outdir);
end

outfile = fullfile(outdir, 'battmo-rate-study-reference.json');
jsonstr = jsonencode(exportStruct, 'PrettyPrint', true);
fid = fopen(outfile, 'w');
assert(fid > 0, 'Could not open %s for writing', outfile);
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s', jsonstr);

fprintf('Wrote %s\n', outfile);

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
