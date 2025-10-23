% Script to validate the calibrated P2D model at different rates

clear all
close all

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

%% Fetch experimental data

% Original data
datafilename = fullfile(getHydra0Dir(), 'rawData', 'TE_1473.mat');
saveddata    = load(datafilename);
dataraw      = saveddata.experiment;

% Calibrated params
filename      = fullfile(getHydra0Dir(), 'parameters', 'equilibrium-calibration-parameters.json');
jsonstructEC  = parseBattmoJson(filename);
filename      = fullfile(getHydra0Dir(), 'parameters', 'high-rate-calibration-parameters.json');
jsonstructHRC = parseBattmoJson(filename);

% Find capacity
input     = struct('lowRateParams', jsonstructEC);
outputCap = runHydra(input, 'runSimulation', false);
cap       = computeCellCapacity(outputCap.model);

fig = figure('Units', 'inches', 'Position', [0.1, 0.1, 8, 6]);
hold on;
colors = lines(numel(dataraw.time));
expname = @(c) sprintf('exp %1.2g C', c);
p2dname = @(c) sprintf('P2D %1.2g C', c);

rates = [0.05, 0.2, 0.5, 1, 2];

for k = 1:numel(dataraw.time)

    expdata = struct('time', dataraw.time{k} * hour, ...
                     'U'   , dataraw.voltage{k}    , ...
                     'I'   , abs(mean(dataraw.current{k})));

    DRate = expdata.I / cap * hour;

    input = struct('lowRateParams' , jsonstructEC , ...
                   'highRateParams', jsonstructHRC, ...
                   'DRate'         , DRate        , ...
                   'totalTime'     , expdata.time(end));

    output = runHydra(input, 'clearSimulation', false);

    figure(fig);
    plot(expdata.time/hour * expdata.I, expdata.U, '--', 'color', colors(k,:));
    hp2d(k) = plot(getTime(output.states)/hour * expdata.I, getE(output.states), 'color', colors(k,:)); %#ok
    drawnow

end

xlabel('C  /  Ah')
ylabel('E  /  V')
axis tight
ylim([3.45, 4.9])

hp(1) = plot(nan, nan, 'k', 'linestyle', '--');
hp(2) = plot(nan, nan, 'k', 'linestyle', '-');
legend(gca(), hp, {'exp', 'P2D'});

legtxt = arrayfun(@(r) {sprintf('%1.2g C', r)}, rates);
ax = axes('position', get(gca(), 'position'), 'visible', 'off');
legend(ax, hp2d, legtxt, 'location', 'sw');

dosave = false;
if dosave
    exportgraphics(fig, 'validation.png', 'resolution', 300);
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
