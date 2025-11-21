%% Script to plot the results from the high-rate calibration

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
printer = @(s) disp(jsonencode(s, 'PrettyPrint', true));

%% Load

% Original data
datafilename = fullfile(getHydra0Dir(), 'rawData', 'TE_1473.mat');
saveddata    = load(datafilename);
dataraw      = saveddata.experiment;

% Highest DRate is last
k = numel(dataraw.time);
expdata = struct('time', dataraw.time{k} * hour, ...
                 'U'   , dataraw.voltage{k}    , ...
                 'I'   , abs(mean(dataraw.current{k})));

neD = 1e-14;
peD = 1e-14;
output1 = load(sprintf('high-rate-calibrated-outputOpt-%g-%g.mat', neD, peD));

neD = 1e-13;
peD = 1e-14;
output2 = load(sprintf('high-rate-calibrated-outputOpt-%g-%g.mat', neD, peD));

%% Plot

colors = lines(4);
fig = figure('Units', 'inches', 'Position', [0.1, 0.1, 8, 6]);
hold on;
plot(expdata.time/hour, expdata.U, 'k--', 'displayname', 'Experiment 2 C');

neD = output1.output0.model.(ne).(co).(am).(sd).referenceDiffusionCoefficient;
peD = output1.output0.model.(pe).(co).(am).(sd).referenceDiffusionCoefficient;
dn = sprintf('Initial guess D_{NE}=%g D_{PE}=%g', neD, peD);
plot(getTime(output1.output0.states)/hour, getE(output1.output0.states), 'color', colors(1,:), 'displayname', dn);

neD = output2.output0.model.(ne).(co).(am).(sd).referenceDiffusionCoefficient;
peD = output2.output0.model.(pe).(co).(am).(sd).referenceDiffusionCoefficient;
dn = sprintf('Initial guess D_{NE}=%g D_{PE}=%g', neD, peD);
plot(getTime(output2.output0.states)/hour, getE(output2.output0.states), 'color', colors(2,:), 'displayname', dn);

plot(getTime(output1.outputOpt.states)/hour, getE(output1.outputOpt.states), 'color', colors(3,:), 'displayname', 'Calibrated');

xlabel('Time / h')
ylabel('E  /  V')

legend()

dosave = true;
if dosave
    exportgraphics(fig, 'high-rate-calibration-results.png', 'Resolution', 300);
end
