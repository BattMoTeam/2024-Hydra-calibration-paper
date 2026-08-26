% Validate the effect of restoring each calibrated parameter to its base value

clear all
close all

mrstDebug(0);

set(0, 'defaultlinelinewidth', 2)
set(0, 'defaulttextfontsize', 15);
set(0, 'defaultaxesfontsize', 15);

getTime = @(states) cellfun(@(s) s.time, states);
getE = @(states) cellfun(@(s) s.Control.E, states);

%% Fetch experimental data and parameter sets

datafilename = fullfile(getHydra0Dir(), 'raw-data', 'TE_1473.mat');
saveddata    = load(datafilename);
dataraw      = saveddata.experiment;

parameterdir = fullfile(getHydra0Dir(), 'parameters');
jsonstruct0 = parseBattmoJson(fullfile(parameterdir, 'h0b-base.json'));

filename      = fullfile(parameterdir, 'equilibrium-calibration-parameters.json');
jsonstructEQC = parseBattmoJson(filename);
filename      = fullfile(parameterdir, 'high-rate-calibration-parameters.json');
jsonstructHRC = parseBattmoJson(filename);

% Each case starts from the complete calibrated parameter set and restores
% exactly one scalar parameter to the corresponding value in h0b-base.
parameterCases = buildParameterCases(jsonstructEQC, jsonstructHRC);

%% Check the nominal rates using the fully calibrated model

input = struct('lowRateParams', jsonstructEQC, ...
               'include_current_collectors', true);
outputCap = runHydra(input, 'runSimulation', false);
calibratedCap = computeCellCapacity(outputCap.model);

rates = [0.05, 0.2, 0.5, 1, 2];
assert(numel(dataraw.time) == numel(rates), ...
       'Expected one experimental data set for each nominal rate.');

for k = 1:numel(dataraw.time)
    experimentalCurrent = abs(mean(dataraw.current{k}));
    DRate = experimentalCurrent / calibratedCap * hour;
    assert(abs(DRate - rates(k))/rates(k) < 0.1, ...
           'DRate %g does not match expected rate %g', DRate, rates(k));
end

%% Validate one restored parameter at a time

for icase = 1:numel(parameterCases)

    parameterCase = parameterCases(icase);
    jsonstructEQC_i = jsonstructEQC;
    jsonstructHRC_i = jsonstructHRC;

    baseValue = getNestedValue(jsonstruct0, parameterCase.path);

    switch parameterCase.parameterSet
      case 'EQC'
        jsonstructEQC_i = setNestedValue(jsonstructEQC_i, ...
                                         parameterCase.path, baseValue);
      case 'HRC'
        jsonstructHRC_i = setNestedValue(jsonstructHRC_i, ...
                                         parameterCase.path, baseValue);
      otherwise
        error('Unknown parameter set %s', parameterCase.parameterSet);
    end

    % EQC parameters affect cell capacity. Recompute it for every case so
    % that each simulation uses the measured experimental current.
    input = struct('lowRateParams', jsonstructEQC_i, ...
                   'include_current_collectors', true);
    outputCap = runHydra(input, 'runSimulation', false);
    cap = computeCellCapacity(outputCap.model);

    fullParameterName = sprintf('%s.%s', parameterCase.parameterSet, ...
                                strjoin(parameterCase.path, '.'));
    parameterName = sprintf('%s.%s', parameterCase.parameterSet, ...
                            getShortParameterPath(parameterCase.path));
    parameterTitle = sprintf('%s: %.4g -> %.4g', ...
                             parameterName, baseValue, ...
                             parameterCase.calibratedValue);
    fprintf('Validating %s (calibrated: %g, base: %g)\n', ...
            fullParameterName, parameterCase.calibratedValue, baseValue);

    fig = figure('Units', 'inches', ...
                 'Position', [0.1, 0.1, 8, 6], ...
                 'Name', parameterTitle, ...
                 'NumberTitle', 'off');
    axPlot = axes('Parent', fig);
    hold(axPlot, 'on');
    colors = lines(numel(dataraw.time));

    RMSE = nan(size(rates));
    hp2d = gobjects(1, numel(rates));
    hp = gobjects(1, 2);

    for k = 1:numel(dataraw.time)

        expdata = struct('time', dataraw.time{k} * hour, ...
                         'U'   , dataraw.voltage{k}    , ...
                         'I'   , abs(mean(dataraw.current{k})));

        DRate = expdata.I / cap * hour;

        input = struct('DRate'                         , DRate             , ...
                       'totalTime'                     , expdata.time(end) , ...
                       'lowRateParams'                 , jsonstructEQC_i   , ...
                       'highRateParams'                , jsonstructHRC_i   , ...
                       'useRegionBruggemanCoefficients', true              , ...
                       'include_current_collectors'    , true);

        output = runHydra(input, 'clearSimulation', false);

        RMSE(k) = l2error(expdata.time, expdata.U, ...
                          getTime(output.states), getE(output.states), ...
                          'extrap', true);

        plot(axPlot, expdata.time/hour * expdata.I, expdata.U, '--', ...
             'color', colors(k,:));
        hp2d(k) = plot(axPlot, getTime(output.states)/hour * expdata.I, ...
                       getE(output.states), 'color', colors(k,:));
        drawnow

    end

    xlabel(axPlot, 'Capacity  /  Ah')
    ylabel(axPlot, 'Voltage  /  V')
    title(axPlot, parameterTitle, 'Interpreter', 'none')
    axis(axPlot, 'tight')
    ylim(axPlot, [3.45, 4.9])

    hp(1) = plot(axPlot, nan, nan, 'k', 'linestyle', '--');
    hp(2) = plot(axPlot, nan, nan, 'k', 'linestyle', '-');
    havg = plot(axPlot, nan, nan, 'linestyle', 'none', 'marker', 'none');
    legend(axPlot, hp, {'exp', 'P2D'});

    legtxt = cell(1, numel(rates) + 1);
    for k = 1:numel(rates)
        legtxt{k} = sprintf('%1.2gC RMSE=%2.1f mV', ...
                            rates(k), RMSE(k)/milli);
    end
    legtxt{end} = sprintf('Avg RMSE=%2.1f mV', mean(RMSE)/milli);

    ax = axes('Parent', fig, ...
              'position', get(axPlot, 'position'), ...
              'visible', 'off');
    legend(ax, [hp2d, havg], legtxt, 'location', 'sw');

end


function parameterCases = buildParameterCases(jsonstructEQC, jsonstructHRC)

    parameterSets = {'EQC', jsonstructEQC; ...
                     'HRC', jsonstructHRC};
    pathsBySet = {getScalarPaths(jsonstructEQC, {}); ...
                  getScalarPaths(jsonstructHRC, {})};
    numCases = sum(cellfun(@numel, pathsBySet));

    emptyCase = struct('parameterSet', '', ...
                       'path', {{}}, ...
                       'calibratedValue', []);
    parameterCases = repmat(emptyCase, 1, numCases);
    icase = 0;

    for iset = 1:size(parameterSets, 1)
        parameterSet = parameterSets{iset, 1};
        jsonstruct = parameterSets{iset, 2};
        paths = pathsBySet{iset};

        for ipath = 1:numel(paths)
            icase = icase + 1;
            path = paths{ipath};
            calibratedValue = getNestedValue(jsonstruct, path);
            assert(isnumeric(calibratedValue) && isscalar(calibratedValue), ...
                   'Expected %s.%s to be a numeric scalar.', ...
                   parameterSet, strjoin(path, '.'));

            parameterCases(icase) = struct( ...
                'parameterSet', parameterSet, ...
                'path', {path}, ...
                'calibratedValue', calibratedValue);
        end
    end

end


function paths = getScalarPaths(jsonstruct, parentPath)

    paths = {};
    fields = fieldnames(jsonstruct);

    for ifield = 1:numel(fields)
        field = fields{ifield};
        value = jsonstruct.(field);
        path = [parentPath, {field}];

        if isstruct(value)
            childPaths = getScalarPaths(value, path);
            paths = [paths, childPaths]; %#ok<AGROW>
        else
            paths{end + 1} = path; %#ok<AGROW>
        end
    end

end


function value = getNestedValue(jsonstruct, path)

    field = path{1};
    assert(isfield(jsonstruct, field), ...
           'Parameter %s is not present in h0b-base.json.', ...
           strjoin(path, '.'));

    if isscalar(path)
        value = jsonstruct.(field);
    else
        value = getNestedValue(jsonstruct.(field), path(2:end));
    end

end


function jsonstruct = setNestedValue(jsonstruct, path, value)

    field = path{1};
    assert(isfield(jsonstruct, field), ...
           'Parameter %s is not present in the calibrated parameter set.', ...
           strjoin(path, '.'));

    if isscalar(path)
        jsonstruct.(field) = value;
    else
        jsonstruct.(field) = setNestedValue(jsonstruct.(field), ...
                                             path(2:end), value);
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
