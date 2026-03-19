% Export BattMo validation curves to JSON for cross-tool comparison.

clear all
close all

mrstDebug(0);

am   = 'ActiveMaterial';
itf  = 'Interface';
pe   = 'PositiveElectrode';
ne   = 'NegativeElectrode';
co   = 'Coating';
sd   = 'SolidDiffusion';
ctrl = 'Control';

getTime = @(states) cellfun(@(s) s.time, states);
getE = @(states) cellfun(@(s) s.(ctrl).E, states);

% Experimental data
datafilename = fullfile(getHydra0Dir(), 'raw-data', 'TE_1473.mat');
saveddata    = load(datafilename);
dataraw      = saveddata.experiment;

% Calibrated parameters
filename      = fullfile(getHydra0Dir(), 'parameters', 'equilibrium-calibration-parameters.json');
jsonstructEC  = parseBattmoJson(filename);
filename      = fullfile(getHydra0Dir(), 'parameters', 'high-rate-calibration-parameters.json');
jsonstructHRC = parseBattmoJson(filename);

% Find model capacity used to define the validation rates
inputCap  = struct('lowRateParams', jsonstructEC, ...
                   'include_current_collectors', true);
outputCap = runHydra(inputCap, 'runSimulation', false);
cap       = computeCellCapacity(outputCap.model);

cases = cell(numel(dataraw.time), 1);

for k = 1:numel(dataraw.time)

    expdata = struct('time', dataraw.time{k} * hour, ...
                     'U'   , dataraw.voltage{k}    , ...
                     'I'   , abs(mean(dataraw.current{k})));

    DRate = expdata.I / cap * hour;

    input = struct('DRate'                         , DRate            , ...
                   'totalTime'                     , expdata.time(end), ...
                   'lowRateParams'                 , jsonstructEC     , ...
                   'highRateParams'                , jsonstructHRC    , ...
                   'useRegionBruggemanCoefficients', true             , ...
                   'include_current_collectors'    , true);

    output = runHydra(input, 'clearSimulation', true);

    simTime = getTime(output.states);
    simVoltage = getE(output.states);

    cases{k} = struct('case_name'    , sprintf('Discharge rate %d', k)     , ...
                      'current_a'    , expdata.I                           , ...
                      'drate'        , DRate                               , ...
                      'experimental' , struct('time_s'   , expdata.time(:)', ...
                                              'voltage_v', expdata.U(:)')  , ...
                      'battmo'       , struct('time_s'   , simTime(:)'     , ...
                                              'voltage_v', simVoltage(:)'));

end

exportStruct = struct('generated_at'      , datestr(now, 'yyyy-mm-ddTHH:MM:SS')                        , ...
                      'source'            , 'scripts/exportValidationReference.m'                      , ...
                      'matlab_release'    , version('-release')                                        , ...
                      'capacity_ah'       , cap                                                         , ...
                      'parameter_files'   , {{'parameters/h0b-base.json', ...
                                               'parameters/equilibrium-calibration-parameters.json', ...
                                               'parameters/high-rate-calibration-parameters.json'}}    , ...
                      'cases'             , {cases});

outdir = fullfile(getHydra0Dir(), 'figures');
if ~isfolder(outdir)
    mkdir(outdir);
end

outfile = fullfile(outdir, 'battmo-validation-reference.json');
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
