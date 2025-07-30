% Script to calibrate parameters under equilibrium assumptions

clear all
close all

mrstModule add ad-core optimization mpfa

mrstDebug(0);

set(0, 'defaultlinelinewidth', 2)
set(0, 'defaulttextfontsize', 15);
set(0, 'defaultaxesfontsize', 15);

am   = 'ActiveMaterial';
itf  = 'Interface';
pe   = 'PositiveElectrode';
ne   = 'NegativeElectrode';
co   = 'Coating';
sd   = 'SolidDiffusion';
ctrl = 'Control';

getTime = @(states) cellfun(@(s) s.time, states);
getE = @(states) cellfun(@(s) s.(ctrl).E, states);
printer = @(s) disp(jsonencode(s, 'PrettyPrint', true));

%% Fetch experimental data

% Smooth voltages
datafilename = fullfile(getHydra0Dir(), 'rawData', 'TE_1473smooth');
saveddata    = load(datafilename);
datasmooth   = saveddata.expsmooth;

% Original data (for currents)
datafilename = fullfile(getHydra0Dir(), 'rawData', 'TE_1473.mat');
saveddata    = load(datafilename);
dataraw      = saveddata.experiment;

% Pack data with lowest rate (first item)
expdata = struct('time' , datasmooth.time{1} * hour                           , ...
                 'E'    , datasmooth.voltage{1}                               , ...
                 'cap'  , abs(trapz(dataraw.time{1}*hour, dataraw.current{1})), ...
                 'I'    , abs(mean(dataraw.current{1}))                       , ...
                 'DRate', datasmooth.crate{1});

%% Initial guess simulation

input = struct('DRate'    , expdata.DRate, ...
               'totalTime', expdata.time(end));
output0 = runHydra(input);

%% Setup and run optimization

EC = EquilibriumCalibration(output0.model, expdata);
[fexp, fcomp] = EC.setupFunctions();

X0 = EC.getDefaultValue();
objective = @(X) EC.objective(X);
v0 = objective(X0);

linIneq.A = -eye(4);
linIneq.b = -0.01*ones(4, 1);

[~, Xopt, history] = unitBoxBFGS(X0, objective, ...
                                 'maximize'    , false  , ...
                                 'linIneq'     , linIneq, ...
                                 'objChangeTol', 1e-11  , ...
                                 'maxit'       , 100    , ...
                                 'logPlot'     , true);
vopt = objective(Xopt);

%% Print

fprintf('\t\tInitial \t Optimized\n');
elde = {'PE', 'NE'};
for e = 1:2
    fprintf('%s\n', elde{e});
    fprintf('theta100\t %1.5f \t %1.5f\n', X0(2*e-1), Xopt(2*e-1));
    fprintf('alpha    \t %1.5f \t %1.5f\n', X0(2*e), Xopt(2*e));
end

fprintf('obj val=%1.2f (%1.2f), iter=%d\n', vopt, v0, numel(history.val));

%% Extract parameters

jsonstructEC = EC.extractAlpha(output0.model, Xopt);
filename = fullfile(getHydra0Dir(), 'parameters', 'equilibrium-calibration-parameters.json');
writeJsonStruct(jsonstructEC, filename);
printer(jsonstructEC);

%% Simulation with calibrated parameters

input = struct('DRate'        , expdata.I * hour / expdata.cap, ...
               'totalTime'    , expdata.time(end)             , ...
               'lowRateParams', jsonstructEC);
outputOpt = runHydra(input);

%% Plot

colors = lines(4);
fig = figure;
hold on
plot(expdata.time/hour, expdata.E, 'k--', 'displayname', 'experiment 0.05 C');
plot(expdata.time/hour, fcomp(expdata.time, X0), 'color', colors(3,:), 'displayname', 'EA initial guess');
plot(expdata.time/hour, fcomp(expdata.time, Xopt), 'color', colors(4,:), 'displayname', 'EA calibrated guess');
plot(getTime(output0.states)/hour, getE(output0.states), 'color', colors(1,:), 'displayname', 'initial guess');
plot(getTime(outputOpt.states)/hour, getE(outputOpt.states), 'color', colors(2,:), 'displayname', 'calibrated');
xlabel 'time / h';
ylabel 'cell voltage / V';
legend('location', 'sw')
axis tight
ylim([3.45, 4.9])

dosave = false;
if dosave
    exportgraphics(fig, 'equilibriumcalibration.png', 'resolution', 300);
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
