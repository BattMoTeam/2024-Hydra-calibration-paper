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
getI = @(states) cellfun(@(s) s.(ctrl).I, states);
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
plot(expdata.time/hour, expdata.U, 'k--', 'displayname', 'Experiment 2C');

neD = output1.output0.model.(ne).(co).(am).(sd).referenceDiffusionCoefficient;
peD = output1.output0.model.(pe).(co).(am).(sd).referenceDiffusionCoefficient;
dn = sprintf('Initial guess D_{NE}=%g D_{PE}=%g', neD, peD);
plot(getTime(output1.output0.states)/hour, getE(output1.output0.states), 'color', colors(1,:), 'displayname', dn);

neD = output2.output0.model.(ne).(co).(am).(sd).referenceDiffusionCoefficient;
peD = output2.output0.model.(pe).(co).(am).(sd).referenceDiffusionCoefficient;
dn = sprintf('Initial guess D_{NE}=%g D_{PE}=%g', neD, peD);
plot(getTime(output2.output0.states)/hour, getE(output2.output0.states), 'color', colors(2,:), 'displayname', dn);

plot(getTime(output1.outputOpt.states)/hour, getE(output1.outputOpt.states), '--', 'color', colors(1,:), 'displayname', 'Calibrated');
plot(getTime(output2.outputOpt.states)/hour, getE(output2.outputOpt.states), '--', 'color', colors(2,:), 'displayname', 'Calibrated');

xlabel('Time  /  h')
ylabel('E  /  V')

legend()

dosave = false;
if dosave
    exportgraphics(fig, '/tmp/high-rate-calibration-results.png', 'Resolution', 300);
end

%% Plot difference between experiment and calibrated

fig = figure('Units', 'inches', 'Position', [0.1, 0.1, 8, 6]);
hold on;

expdataUinterp1 = @(t) interp1(expdata.time, expdata.U, t, 'linear', 'extrap');

neD = output1.output0.model.(ne).(co).(am).(sd).referenceDiffusionCoefficient;
peD = output1.output0.model.(pe).(co).(am).(sd).referenceDiffusionCoefficient;
dn = sprintf('Difference with D_{NE}=%g', neD);
plot(getTime(output1.outputOpt.states)/hour, abs(getE(output1.outputOpt.states) - expdataUinterp1(getTime(output1.outputOpt.states))), 'color', colors(1,:), 'displayname', dn);

neD = output2.output0.model.(ne).(co).(am).(sd).referenceDiffusionCoefficient;
peD = output2.output0.model.(pe).(co).(am).(sd).referenceDiffusionCoefficient;
dn = sprintf('Difference with D_{NE}=%g', neD);
plot(getTime(output2.outputOpt.states)/hour, abs(getE(output2.outputOpt.states) - expdataUinterp1(getTime(output2.outputOpt.states))), 'color', colors(2,:), 'displayname', dn);

xlabel('Time  /  h')
ylabel('|E_{exp} - E_{sim}|  /  V')

legend();

%% Plot dqdv: may need spline approx of OCVs

fig = figure %('Units', 'inches', 'Position', [0.1, 0.1, 8, 6]);
hold on;
expQ = cumtrapz(expdata.time, expdata.I*ones(size(expdata.time))) / hour;

dqdvExp = gradient(expQ, expdata.U);
plot(expdata.U, -dqdvExp, 'k-', 'displayname', 'Experiment 2C');

simQ = @(states) cumtrapz(getTime(states), getI(states)) / hour;

dqdv = gradient(simQ(output1.outputOpt.states), getE(output1.outputOpt.states));
plot(getE(output1.outputOpt.states), -dqdv, 'color', colors(1,:), 'displayname', 'Calibrated D_{NE}=1e-14');

dqdv = gradient(simQ(output2.outputOpt.states), getE(output2.outputOpt.states));
plot(getE(output2.outputOpt.states), -dqdv, 'color', colors(2,:), 'displayname', 'Calibrated D_{NE}=1e-13');

xlabel 'E  /  V'
ylabel '-dQdV  /  Ah V^{-1}'

%% Plot dvdq: may need spline approx of OCVs

fig = figure%('Units', 'inches', 'Position', [0.1, 0.1, 8, 6]);
hold on;
dvdqExp = gradient(expdata.U, expQ);
plot(expQ, -dvdqExp, 'k--', 'displayname', 'Experiment 2C');

dvdq = gradient(getE(output1.outputOpt.states), simQ(output1.outputOpt.states));
plot(simQ(output1.outputOpt.states), -dvdq, 'color', colors(1,:), 'displayname', 'Calibrated D_{NE}=1e-14');
dvdq = gradient(getE(output2.outputOpt.states), simQ(output2.outputOpt.states));
plot(simQ(output2.outputOpt.states), -dvdq, 'color', colors(2,:), 'displayname', 'Calibrated D_{NE}=1e-13');
ylim([0, max(abs(dvdqExp))*0.2]);

xlabel 'Q  /  Ah'
ylabel '-dVdQ  /  V A^{-1} h^{-1}'

%% Calculate tortuosity

tau1init = tortuosity(output1.output0.model);
tau1opt = tortuosity(output1.outputOpt.model);

tau2init = tortuosity(output2.output0.model);
tau2opt = tortuosity(output2.outputOpt.model);

% Add to jsonstructHRC
output1.jsonstructHRC.(ne).tortuosity = tau1opt.(ne);
output1.jsonstructHRC.(pe).tortuosity = tau1opt.(pe);
output2.jsonstructHRC.(ne).tortuosity = tau2opt.(ne);
output2.jsonstructHRC.(pe).tortuosity = tau2opt.(pe);


%% Compare
jsonDiff(output1.jsonstructHRC, ...
         output2.jsonstructHRC, ...
         sprintf('neD %g', output1.output0.model.(ne).(co).(am).(sd).referenceDiffusionCoefficient), ...
         sprintf('neD %g', output2.output0.model.(ne).(co).(am).(sd).referenceDiffusionCoefficient));
