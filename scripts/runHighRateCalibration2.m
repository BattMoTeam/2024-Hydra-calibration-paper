% Script to calibrate parameters using high-rate data

clear all
close all

mrstModule add ad-core optimization mpfa

%mrstDebug(99);

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

getTime = @(states) cellfun(@(s) s.time, states);
getE = @(states) cellfun(@(s) s.(ctrl).E, states);
printer = @(s) disp(jsonencode(s, 'PrettyPrint', true));

debug = true;

%% Fetch experimental data

% Original data
datafilename = fullfile(getHydra0Dir(), 'rawData', 'TE_1473.mat');
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

% Find capacity
input     = struct('lowRateParams', jsonstructEC);
outputCap = runHydra(input, 'runSimulation', false);
cap       = computeCellCapacity(outputCap.model);


% baseJson = parseBattmoJson(fullfile(getHydra0Dir(), 'parameters', 'h0b-base.json'));

% baseJson.(elyte).useRegionBruggemanCoefficients = true;
% baseJson.(elyte).regionBruggemanCoefficients = struct(ne, 4., ...
%                                                       pe, 6., ...
%                                                       sep, 1.5);
% baseJson.(elyte) = rmfield(baseJson.(elyte), 'bruggemanCoefficient');

% baseJsonName = fullfile('/tmp', 'baseJson.json');
% writeJsonStruct(baseJson, baseJsonName);

% Initial guess
input = struct('DRate'        , expdata.I / cap * hour, ...
               'totalTime'    , expdata.time(end)     , ...
               'lowRateParams', jsonstructEC);
output0 = runHydra(input);

plotDashboard(output0.model, output0.states, 'step', 22)
return

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

return

%% Setup optimization

% Evaluate experimental data at simulation times (allow for
% extrapolation since expdata.time(end) is very close to
% output.states{end}.time)
simtimes = getTime(output0.states);
assert(expdata.time(1) <= simtimes(1));
assert(abs(expdata.time(end) - simtimes(end)) < 1e-11, 'expdata end time %g but sim end time %g', expdata.time(end)/hour, simtimes(end)/hour);

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
    xlabel('time / h')
    ylabel('potential / V')
    title('statesExp')
    drawnow
end

simulatorSetup = struct('model', output0.model, ...
                        'schedule', output0.schedule, ...
                        'state0', output0.initstate);

% % Setup parameters to be calibrated
% HRC = HighRateCalibration(simulatorSetup);
% parameters = HRC.getParams();

paramsetter  = ParamSetter();

getValues = @(model, notused) paramsetter.getValues(model);
setValues = @(model, notused, v) paramsetter.setValues(model, v);

boxLims = [[1e4, 1e7];
           [1e4, 1e7];
           [1.5, 20];
           [1.5, 20];
           [1e-15, 1e-11];
           [1e-15, 1e-11];
           [1.0, 10];
           [1.0, 10]
          ];

parameters{1} = ModelParameter(simulatorSetup            , ...
                               'name'     , 'lengthPoros', ...
                               'belongsTo', 'model'      , ...
                               'location' , {''}         , ...
                               'boxLims'  , boxLims      , ...
                               'getfun'   , getValues    , ...
                               'setfun'   , setValues);



% Objective function
lsq = @(model, states, schedule, varargin) leastSquaresEI(model, states, statesExp, schedule, varargin{:});
v = lsq(simulatorSetup.model, output0.states, simulatorSetup.schedule);
scaling = sum([v{:}]);

objective = @(p, varargin) evalObjectiveBattmo(p, lsq, simulatorSetup, parameters, ...
                                               'objScaling', scaling, varargin{:});

if debug
    % The least squares function evaluated at the experimental values
    % should be zero
    v = lsq(output0.model, statesExp, simulatorSetup.schedule);
    assert(norm([v{:}]) == 0.0);

    % Compare gradients calculated using adjoints and finite
    % difference approximation
    Xtmp = getScaledParameterVector(simulatorSetup, parameters);

    [vad, gad] = evalObjectiveBattmo(Xtmp, lsq, simulatorSetup, parameters, ...
                                     'gradientMethod', 'AdjointAD', ...
                                     'objScaling', scaling);

    [vnum, gnum] = evalObjectiveBattmo(Xtmp, lsq, simulatorSetup, parameters, ...
                                       'gradientMethod', 'PerturbationADNUM', ...
                                       'PerturbationSize', 1e-7, ...
                                       'objScaling', scaling);
    assert(abs(vad - vnum) < eps);
    assert(all(abs(gad) > 0));
    assert(all(abs(gnum) > 0));
    assert(norm((gad-gnum)./gnum, 'inf') < 1e-5);
end

%% Run optimization

X0 = getScaledParameterVector(simulatorSetup, parameters);
v0 = objective(X0);

callbackfunc = @(history, it) callbackplot(history, it, simulatorSetup, parameters, statesExp, ...
                                           'plotEveryIt', 1, ...
                                           'objScaling', scaling);

[vopt, Xopt, history] = unitBoxBFGS(X0, objective, ...
                                    'objChangeTol', 1e-8 , ...
                                    'maximize'    , false, ...
                                    'maxit'       , 1  , ...
                                    'logPlot'     , true, ...
                                    'callbackfunc', callbackfunc);

setupOpt = updateSetupFromScaledParameters(simulatorSetup, parameters, Xopt);

fprintf('obj val=%1.2f (%1.2f), iter=%d\n', vopt, v0, numel(history.val));

%% Extract parameters

vals0 = parameters{1}.getParameter(simulatorSetup);
paramsetter.print(vals0);
valsOpt = parameters{1}.getParameter(setupOpt);
paramsetter.print(valsOpt);

return


jsonstructHRC = HRC.export(setupOpt);
filename = fullfile(getHydra0Dir(), 'parameters', 'high-rate-calibration-parameters.json');
writeJsonStruct(jsonstructHRC, filename);
disp(jsonencode(jsonstructHRC, 'PrettyPrint', true));

%% Run model with calibrated parameters

input = struct('DRate'         , expdata.I / cap * hour    , ...
               'totalTime'     , expdata.time(end)         , ...
               'lowRateParams' , jsonstructEC, ...
               'highRateParams', jsonstructHRC);
outputOpt = runHydra(input);

% %% run with original bruggeman (from report)
% jsonstructHRCtmp = jsonstructHRC;
% jsonstructHRCtmp.(pe).(co).bruggemanCoefficient = 3.0;
% jsonstructHRCtmp.(ne).(co).bruggemanCoefficient = 3.46;
% input = struct('DRate'         , expdata.I / cap * hour    , ...
%                'totalTime'     , expdata.time(end)         , ...
%                'lowRateParams' , jsonstructEC, ...
%                'highRateParams', jsonstructHRCtmp);
% outputBman = runHydra(input);

%% run with cc

% input = struct('DRate'                     , expdata.I / cap * hour, ...
%                'totalTime'                 , expdata.time(end)     , ...
%                'lowRateParams' , jsonstructEC, ...
%                'highRateParams', jsonstructHRC, ...
%                'include_current_collectors', true);
% outputCC = runHydra(input);


%% Plot

colors = lines(4);
fig = figure;
hold on;
plot(expdata.time/hour, expdata.U, 'k--', 'displayname', 'experiment 2 C');
plot(getTime(output0.states)/hour, getE(output0.states), 'color', colors(1,:), 'displayname', 'initial guess')
plot(getTime(outputOpt.states)/hour, getE(outputOpt.states), 'color', colors(2,:), 'displayname', 'calibrated');
%plot(getTime(outputBman.states)/hour, getE(outputBman.states), 'color', colors(3,:), 'displayname', 'org Bman');
%plot(getTime(outputCC.states)/hour, getE(outputCC.states), 'color', colors(4,:), 'displayname', 'with CC');

xlabel('time / h')
ylabel('voltage / V')
legend('location', 'sw')
axis tight
ylim([3.45, 4.9])

dosave = false;
if dosave
    exportgraphics(fig, 'high-rate-calibration-without-bruggeman-with-CC.png', 'resolution', 300)
end


return

%% Play with Bruggemen

bman_pe = linspace(3.46, 7, 5);
bman_ne = linspace(3, 4, 5);

for k = 1:numel(bman_pe)
    for m = 1:numel(bman_ne)
        jsonstruct = outputOpt.jsonstruct;

        jsonstruct.(pe).(co).bruggemanCoefficient = bman_pe(k);
        jsonstruct.(ne).(co).bruggemanCoefficient = bman_ne(m);

        name = tempname;
        writeJsonStruct(jsonstruct, name);

        input = struct('DRate'         , expdata.I / cap * hour    , ...
                       'totalTime'     , expdata.time(end)         , ...
                       'lowRateParams' , jsonstructEC, ...
                       'highRateParams', jsonstructHRC, ...
                       'baseJson', name);

        outputTmp = runHydra(input);

        % plot
        figure(fig); hold on;
        s = sprintf('PE %1.2f NE %1.2f', jsonstruct.(pe).(co).bruggemanCoefficient, jsonstruct.(ne).(co).bruggemanCoefficient);
        plot(getTime(outputTmp.states)/hour, getE(outputTmp.states), 'displayname', s);
        drawnow

    end
end

%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
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
