%% Script to calibrate parameters using high-rate data

clear all
close all

diary(sprintf('_diary-%s-%s.txt', mfilename, datetime('now', 'Format', 'yyyyMMdd-HHmmss')));

set(0, 'defaultlinelinewidth', 2)
set(0, 'defaulttextfontsize', 15);
set(0, 'defaultaxesfontsize', 15);

am    = 'ActiveMaterial';
itf   = 'Interface';
pe    = 'PositiveElectrode';
ne    = 'NegativeElectrode';
co    = 'Coating';
sd    = 'SolidDiffusion';
ctrl  = 'Control';
elyte = 'Electrolyte';
sep   = 'Separator';

% mrstDebug(0);

doplot = true;
debug = false;

getTime = @(states) cellfun(@(s) s.time, states);
getE = @(states) cellfun(@(s) s.(ctrl).E, states);
printer = @(s) disp(jsonencode(s, 'PrettyPrint', true));

%% Fetch experimental data

datafilename = fullfile(getHydra0Dir(), 'raw-data', 'TE_1473.mat');
saveddata    = load(datafilename);
dataraw      = saveddata.experiment;

% Highest DRate is last
k = numel(dataraw.time);
expdata = struct('time', dataraw.time{k} * hour, ...
                 'U'   , dataraw.voltage{k}    , ...
                 'I'   , abs(mean(dataraw.current{k})));

%% Initial guess using equilibrium calibration data

filename     = fullfile(getHydra0Dir(), 'parameters', 'equilibrium-calibration-parameters.json');
jsonstructEC = parseBattmoJson(filename);

shortnames = {'ne_vsa', 'pe_vsa', 'ne_D', 'pe_D', 'elyte_bgfactor'};
disp('shortnames:');
printer(shortnames);
useRegionBruggemanCoefficients = any(strcmp(shortnames, 'elyte_bgfactor'));

numTimesteps = 100; % 400

input0 = struct('I'                             , expdata.I, ...
                'totalTime'                     , expdata.time(end)             , ...
                'numTimesteps'                  , numTimesteps                  , ...
                'lowRateParams'                 , jsonstructEC                  , ...
                'useRegionBruggemanCoefficients', useRegionBruggemanCoefficients, ...
                'include_current_collectors'    , true);
output0 = runHydra(input0, 'clearSimulation', false);

% Avoid adaptive time stepping to ensure adjoint consistency
output0.nls.timeStepSelector = SimpleTimeStepSelector();
output0.nls.maxTimestepCuts = 0;

if debug
    % Check how exp and initial guess compare
    figure; hold on; grid on;
    plot(expdata.time/hour, expdata.U, 'k--');
    plot(getTime(output0.states)/hour, getE(output0.states));
    xlabel('time / h')
    ylabel('potential / V')
    title('initial guess')
    drawnow
end

%% Setup optimization

% Evaluate experimental data at simulation times (allow for
% extrapolation since expdata.time(end) is very close to
% output.states{end}.time)
simtimes = getTime(output0.states);
assert(expdata.time(1) <= simtimes(1));
assert(abs(expdata.time(end) - simtimes(end)) < 1e-11);

Evals     = interp1(expdata.time, expdata.U, simtimes, 'linear', 'extrap');
statesExp = cell(numel(output0.states), 1);

for k = 1:numel(output0.states)
    statesExp{k}.time     = simtimes(k);
    statesExp{k}.(ctrl).E = Evals(k);
end

if debug
    % Check that the extracted values are the same as the raw values
    figure; hold on; grid on;
    plot(expdata.time/hour, expdata.U, 'k--');
    plot(getTime(statesExp)/hour, getE(statesExp));
    xlabel('Time / h')
    ylabel('Potential / V')
    title('statesExp')
    drawnow
end

simulatorSetup = SimulationSetup(struct('model'          , output0.model   , ...
                                        'schedule'       , output0.schedule, ...
                                        'initstate'      , output0.initstate, ...
                                        'NonLinearSolver', output0.nls     , ...
                                        'OutputMinisteps', false));

% Setup parameters to be calibrated
HRC = HighRateCalibration(simulatorSetup, 'shortnames', shortnames);
parameters = HRC.getParams();

% Objective function
lsq = @(simsetup, states, varargin) leastSquaresEI(simsetup, states, statesExp, varargin{:});
v = lsq(simulatorSetup, output0.states);
scaling = sum([v{:}]);
objective = @(p, varargin) evalObjectiveBattmo(p, lsq, simulatorSetup, parameters, ...
                                               'objScaling', scaling, varargin{:});

% Compute and classify sensitivities at the initial parameter values.
X0 = getScaledParameterVector(simulatorSetup, parameters);
initialParameterValues = cellfun( ...
    @(parameter) parameter.getParameterValue(simulatorSetup), parameters);
sensitivityReport = computeSensitivities( ...
    X0, objective, HRC.shortnames, scaling);

if debug
    % The least squares function evaluated at the experimental values
    % should be zero
    v = lsq(simulatorSetup, statesExp);
    assert(norm([v{:}]) == 0.0);

    % Compare gradients calculated using adjoints and finite
    % difference approximation
    disp('Gradient comparison at initial parameters:');
    compareAdjointAndFiniteDifferenceGradients(X0, objective, HRC.shortnames);

    %return
end

%% Run optimization

v0 = sensitivityReport.objectiveValue;

callbackfunc = @(history, it) callbackplot(history, it, simulatorSetup, parameters, expdata, ...
                                           'plotEveryIt', 10     , ...
                                           'objScaling' , scaling, ...
                                           'doplot'     , doplot);

gradTol = 1e-3;
objChangeTol = 1e-10;
maxit = 100;
[vopt, Xopt, history] = unitBoxBFGS(X0, objective, ...
                                    'gradTol'         , gradTol     , ...
                                    'objChangeTol'    , objChangeTol, ...
                                    'lineSearchMaxIt' , 10          , ...
                                    'maxInitialUpdate', 0.02        , ...
                                    'maximize'        , false       , ...
                                    'maxit'           , maxit       , ...
                                    'logPlot'         , true        , ...
                                    'callbackfunc'    , callbackfunc, ...
                                    'plotEvolution'   , doplot);

setupOpt = updateSetupFromScaledParameters(simulatorSetup, parameters, Xopt);

fprintf('obj val=%1.2f (%1.2f), iter=%d\n', vopt, v0, numel(history.val));
reasonStr = getReasonStr(history, ...
                         'gradTol'     , gradTol     , ...
                         'objChangeTol', objChangeTol, ...
                         'maxit'       , maxit);
disp(reasonStr);

if debug && numel(history.val) >= 2 && ...
        abs(history.val(end) - history.val(end-1)) < objChangeTol

    % Calculate fd and adjoint gradients at final point
    disp('Gradient comparison at optimized parameters:');
    compareAdjointAndFiniteDifferenceGradients(Xopt, objective, HRC.shortnames);
end

%% Extract parameters

jsonstructHRC = HRC.export(setupOpt);
filename = fullfile(getHydra0Dir(), 'parameters', 'high-rate-calibration-parameters.json');
writeStruct(jsonstructHRC, filename);
printer(jsonstructHRC);

Dne = output0.model.(ne).(co).(am).(sd).referenceDiffusionCoefficient;
Dpe = output0.model.(pe).(co).(am).(sd).referenceDiffusionCoefficient;
filename = fullfile(getHydra0Dir(), 'parameters', sprintf('high-rate-calibration-parameters-%g-%g.json', Dne, Dpe));
writeStruct(jsonstructHRC, filename);

%% Run model with calibrated parameters

inputOpt = struct('I'                             , expdata.I                     , ...
                  'totalTime'                     , expdata.time(end)             , ...
                  'numTimesteps'                  , numTimesteps                  , ...
                  'lowRateParams'                 , jsonstructEC                  , ...
                  'highRateParams'                , jsonstructHRC                 , ...
                  'useRegionBruggemanCoefficients', useRegionBruggemanCoefficients, ...
                  'include_current_collectors'    , true);
outputOpt = runHydra(inputOpt, 'clearSimulation', false);

assert(outputOpt.model.Electrolyte.bgfactor == setupOpt.model.Electrolyte.bgfactor);

%% Quantify differences

vfinal = lsq(simulatorSetup, outputOpt.states);

expdataUinterp1 = @(t) interp1(expdata.time, expdata.U, t, 'linear', 'extrap');
RMSE = l2error(getTime(outputOpt.states), getE(outputOpt.states), expdata.time, expdata.U, 'extrap', true);

fprintf('Final least squares values:\n');
fprintf('vopt: %g\n', vopt);
fprintf('Sum of squares: %g\n', sum([vfinal{:}]));
fprintf('RMSE: %g mV\n', RMSE/milli);

if doplot
    % plot differences
    figure; hold on; grid on;
    plot(getTime(outputOpt.states), (getE(outputOpt.states) - expdataUinterp1(getTime(outputOpt.states))).^2, 'displayname', '|E_{sim} - E_{exp}|^2');
    plot(getTime(outputOpt.states), [vfinal{:}], 'displayname', 'vfinal');
end

%% Plot

if doplot
    colors = lines(2);
    fig = figure('Units', 'inches', 'Position', [0.1, 0.1, 8, 6]);
    hold on;
    plot(expdata.time/hour, expdata.U, 'k--', 'displayname', 'Experiment 2C');
    plot(getTime(output0.states)/hour, getE(output0.states), 'color', colors(1,:), 'displayname', 'Initial guess')
    plot(getTime(outputOpt.states)/hour, getE(outputOpt.states), 'color', colors(2,:), 'displayname', 'Calibrated');
    xlabel('Time  /  h')
    ylabel('E  /  V')
    legend('location', 'sw')
    axis tight
    ylim([3.45, 4.9])

    dosave = true;
    if dosave
        exportgraphics(fig, 'high-rate-calibration.png', 'resolution', 300)
    end
end

%% Quantify difference between experiment and calibrated

RMSE = l2error(getTime(outputOpt.states), getE(outputOpt.states), expdata.time, expdata.U, 'extrap', true);
fprintf('RMSE after calibration: %g mV\n', RMSE/milli);

%% Print

disp('Results HRC');
printer(jsonstructHRC);

% Print tortuosities using bgfactor: we have eff cond = bgfactor *
% cond * poro^bman = cond * poro^(bman*lg(bgfactor)) = cond *
% poro^(eff bman). From this we can calculate the tortuosities as tau
% = poro^-bman.

model = outputOpt.model;
rbc = model.(elyte).regionBruggemanCoefficients;
bgfactor = model.(elyte).bgfactor;
effbman = @(poro, bman) bman * log(bgfactor) / log(poro);
effbmen = struct();
effbmen.(ne) = effbman(1-model.(ne).(co).volumeFraction, rbc.(ne));
effbmen.(pe) = effbman(1-model.(pe).(co).volumeFraction, rbc.(pe));
effbmen.(sep) = effbman(model.(sep).porosity, rbc.(sep));
disp('Pseudo Bruggeman coefficients:');
printer(effbmen);

% Convert to tortuosities
tortuosity = @(vf, bman) vf.^(-bman);
tau = struct();
tau.(ne) = tortuosity(1-model.(ne).(co).volumeFraction, effbmen.(ne));
tau.(pe) = tortuosity(1-model.(pe).(co).volumeFraction, effbmen.(pe));
tau.(sep) = tortuosity(model.(sep).porosity, effbmen.(sep));
disp('Derived tortuosities');
printer(tau);

% % Postprocess: Report effective electrode conductivities and
% % electrolyte tortuosities
% regionTags = outputOpt.model.(elyte).regionTags;
% vfs = struct();
% vfs.(ne)  = unique(outputOpt.model.(elyte).volumeFraction(regionTags == 1));
% vfs.(pe)  = unique(outputOpt.model.(elyte).volumeFraction(regionTags == 2));
% vfs.(sep) = unique(outputOpt.model.(elyte).volumeFraction(regionTags == 3));
% assert(isscalar(vfs.(ne)));
% assert(isscalar(vfs.(pe)));
% assert(isscalar(vfs.(sep)));

% if useRegionBruggemanCoefficients
%     effectiveRegionBruggeman = outputOpt.rbc;
% else
%     bg = outputOpt.model.(elyte).bruggemanCoefficient;
%     assert(isscalar(bg));
%     effectiveRegionBruggeman = struct(ne, bg, pe, bg, sep, bg);
% end

% tortuosityParams = struct();
% tortuosityParams.(elyte).regionBruggemanCoefficients = effectiveRegionBruggeman;
% tau = calculateTortuosityFromBruggeman(vfs, tortuosityParams);
% disp('Effective electrolyte region Bruggeman coefficients:');
% printer(effectiveRegionBruggeman);
% disp('Tortuosities:');
% printer(tau);

effCond = struct(pe, outputOpt.model.(pe).(co).effectiveElectronicConductivity, ...
                 ne, outputOpt.model.(ne).(co).effectiveElectronicConductivity);
disp('Effective electronic conductivities:');
printer(effCond);

% For testing: print NE volumetric surface area
fprintf('Initial diffusion Dne=%g volumetricsurfacearea=%g\n', ...
        Dne, jsonstructHRC.(ne).(co).(am).(itf).volumetricSurfaceArea);
fprintf('RMSE %g mV\n', RMSE/milli);
disp(reasonStr)

% Print initial and final vals, plus sensitivities
finalParameterValues = cellfun( ...
    @(parameter) parameter.getParameterValue(setupOpt), parameters);
sensitivitySummary = table( ...
    HRC.shortnames(:), initialParameterValues(:), finalParameterValues(:), ...
    sensitivityReport.initialGroup(:), ...
    'VariableNames', ...
    {'Shortname', 'InitialValue', 'FinalValue', 'InitialGroup'});
fprintf('\nInitial sensitivity classification and calibration results:\n');
disp(sensitivitySummary);

diary off;



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
