%%

clear all
% close all

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
                 'U'    , datasmooth.voltage{1}                               , ...
                 'cap'  , abs(trapz(dataraw.time{1}*hour, dataraw.current{1})), ...
                 'I'    , abs(mean(dataraw.current{1}))                       , ...
                 'DRate', datasmooth.crate{1});


%% P2D simulation

% Initial
rate = expdata.I / expdata.cap * hour;
input0 = struct('DRate'    , rate, ...
                'totalTime', expdata.time(end));
output0 = runHydra(input0, 'clearSimulation', false);

% Calibrated
jsonstructEC = parseBattmoJson(fullfile(getHydra0Dir(), 'parameters', 'equilibrium-calibration-parameters.json'));
jsonstructHRC = parseBattmoJson(fullfile(getHydra0Dir(), 'parameters', 'high-rate-calibration-parameters.json'));

inputOpt = struct('DRate'        , rate                        , ...
                  'totalTime'    , expdata.time(end)             , ...
                  'lowRateParams', jsonstructEC,  ...
                  'highRateParams', jsonstructHRC);

outputOpt = runHydra(inputOpt, 'clearSimulation', false);

%% Plot half cell OCPs and full cell discharge

%fig = figure('Units', 'inches', 'Position', [0.1, 0.1, 8, 6]);
fig = figure;
hold on

% Full cell
capexp = cumtrapz(expdata.time, expdata.I * ones(size(expdata.time)));
cap0 = cumtrapz(getTime(output0.states), getI(output0.states));
capOpt = cumtrapz(getTime(outputOpt.states), getI(outputOpt.states));

plot(capexp / hour, expdata.U, 'k--', 'DisplayName', 'Exp');
plot(cap0 / hour, getE(output0.states), 'b-', 'DisplayName', 'P2D initial');
plot(capOpt / hour, getE(outputOpt.states), 'r-', 'DisplayName', 'P2D calibrated');
title('Full cell discharge');

% normalize capacity to 1
colors = lines(5);

capmaxes = [capexp(end), cap0(end), capOpt(end)];
[~, imax] = max(capmaxes);
capmax = capmaxes(imax);

msz = 16;

fig = figure;
hold on
plot(capexp / capmax, expdata.U, 'k--', 'DisplayName', 'Exp');
plot(capexp(1) / capmax, expdata.U(1), 'kx', 'MarkerSize', 10, 'DisplayName', 'Exp start');
plot(capexp(end) / capmax, expdata.U(end), 'ko', 'MarkerSize', 10, 'DisplayName', 'Exp end');

plot(cap0 / capmax, getE(output0.states), 'color', colors(1,:), 'DisplayName', 'P2D initial');
plot(cap0(1) / capmax, getE(output0.states(1)), 'x', 'color', colors(1,:), 'MarkerSize', 14, 'DisplayName', 'P2D initial start');
plot(cap0(end) / capmax, getE(output0.states(end)), 'o', 'color', colors(1,:), 'MarkerSize', 14, 'DisplayName', 'P2D initial end');

plot(capOpt / capmax, getE(outputOpt.states), 'color', colors(2,:), 'DisplayName', 'P2D calibrated');
plot(capOpt(1) / capmax, getE(outputOpt.states(1)), 'x', 'color', colors(2,:), 'MarkerSize', 14, 'DisplayName', 'P2D calibrated start');
plot(capOpt(end) / capmax, getE(outputOpt.states(end)), 'o', 'color', colors(2,:), 'MarkerSize', 14, 'DisplayName', 'P2D calibrated end');

% Half cell
json0 = parseBat
tmoJson(fullfile(getHydra0Dir(), 'parameters', 'h0b-base.json'));
jsonstructEC = parseBattmoJson(fullfile(getHydra0Dir(), 'parameters', 'equilibrium-calibration-parameters.json'));
jsonOpt = mergeJsonStructs({jsonstructEC, json0});
jsons = {json0, jsonOpt};
tags = {'Initial', 'Calibrated'};
linestyles = {'--', '-'};

for k = 1:2
    jsonstruct = jsons{k};

    ne_itf = jsonstruct.(ne).(co).(am).(itf);
    ne_cmax = ne_itf.saturationConcentration;
    ne_c = linspace(ne_itf.guestStoichiometry0, ne_itf.guestStoichiometry100, 100) * ne_cmax;
    ne_ocp = computeOCPanodeH0b(ne_c, [], ne_cmax);

    pe_itf = jsonstruct.(pe).(co).(am).(itf);
    pe_cmax = pe_itf.saturationConcentration;
    pe_c = linspace(pe_itf.guestStoichiometry0, pe_itf.guestStoichiometry100, 100) * pe_cmax;
    pe_ocp = computeOCPcathodeH0b(pe_c, [], pe_cmax);

    plot(1-ne_c/ne_cmax, ne_ocp, 'color', colors(3,:), 'linestyle', linestyles{k}, 'DisplayName', sprintf('NE OCP %s', tags{k}));
    plot(1-ne_c(1)/ne_cmax, ne_ocp(1), 'x', 'color', colors(3,:), 'MarkerSize', msz, 'DisplayName', sprintf('NE OCP %s start', tags{k}));
    plot(1-ne_c(end)/ne_cmax, ne_ocp(end), 'o', 'color', colors(3,:), 'MarkerSize', msz, 'DisplayName', sprintf('NE OCP %s end', tags{k}));

    plot(pe_c/pe_cmax, pe_ocp, 'color', colors(4,:), 'linestyle', linestyles{k}, 'DisplayName', sprintf('PE OCP %s', tags{k}));
    plot(pe_c(1)/pe_cmax, pe_ocp(1), 'x', 'color', colors(4,:), 'MarkerSize', msz, 'DisplayName', sprintf('PE OCP %s start', tags{k}));
    plot(pe_c(end)/pe_cmax, pe_ocp(end), 'o', 'color', colors(4,:), 'MarkerSize', msz, 'DisplayName', sprintf('PE OCP %s end', tags{k}));

end

%legend

%% Plot OCP difference after calibration

json0 = parseBattmoJson(fullfile(getHydra0Dir(), 'parameters', 'h0b-base.json'));
jsonstructEC = parseBattmoJson(fullfile(getHydra0Dir(), 'parameters', 'equilibrium-calibration-parameters.json'));
jsonOpt = mergeJsonStructs({jsonstructEC, json0});

jsonstruct = jsonOpt;



fig = figure;
hold on
grid on
colors = lines(3);

capOpt = cumtrapz(getTime(outputOpt.states), getI(outputOpt.states));
plot(capOpt/capOpt(end), getE(outputOpt.states), 'color', 'k', 'DisplayName', 'P2D initial');

ne_itf = jsonstruct.(ne).(co).(am).(itf);
ne_cmax = ne_itf.saturationConcentration;
ne_c = linspace(ne_itf.guestStoichiometry0, ne_itf.guestStoichiometry100, 100) * ne_cmax;
ne_ocp = computeOCPanodeH0b(ne_c, [], ne_cmax);

pe_itf = jsonstruct.(pe).(co).(am).(itf);
pe_cmax = pe_itf.saturationConcentration;
pe_c = linspace(pe_itf.guestStoichiometry0, pe_itf.guestStoichiometry100, 100) * pe_cmax;
pe_ocp = computeOCPcathodeH0b(pe_c, [], pe_cmax);

ocp = pe_ocp - ne_ocp;

x = (pe_c - min(pe_c)) / (max(pe_c)-min(pe_c)); % / pe_cmax;
plot(x, ocp, 'color', colors(2,:), 'DisplayName', 'OCP diff');


%% Plot OCP vs lithiation

json0 = parseBattmoJson(fullfile(getHydra0Dir(), 'parameters', 'h0b-base.json'));
jsonstructEC = parseBattmoJson(fullfile(getHydra0Dir(), 'parameters', 'equilibrium-calibration-parameters.json'));
jsonOpt = mergeJsonStructs({jsonstructEC, json0});

jsonstruct = jsonOpt;
jsonstruct.(ne).(co).(am).(itf).guestStoichiometry0 = -0.126985;

ne_itf = jsonstruct.(ne).(co).(am).(itf);
ne_cmax = ne_itf.saturationConcentration;
ne_theta = linspace(ne_itf.guestStoichiometry0, ne_itf.guestStoichiometry100, 100);
ne_ocp = computeOCPanodeH0b(ne_theta * ne_cmax, [], ne_cmax);

pe_itf = jsonstruct.(pe).(co).(am).(itf);
pe_cmax = pe_itf.saturationConcentration;
pe_theta = linspace(pe_itf.guestStoichiometry0, pe_itf.guestStoichiometry100, 100);
pe_ocp = computeOCPcathodeH0b(pe_theta * pe_cmax, [], pe_cmax);

fig = figure;

tiledlayout(2,2);
nexttile
hold on; grid on
plot(ne_theta, ne_ocp, 'DisplayName', 'NE OCP');
plot(pe_theta, pe_ocp, 'DisplayName', 'PE OCP');
xlabel 'Lithiation \theta';
ylabel 'OCP';
legend

nexttile
hold on; grid on
plot(ne_theta*ne_cmax*litre, ne_ocp, 'DisplayName', 'NE OCP');
plot(pe_theta*pe_cmax*litre, pe_ocp, 'DisplayName', 'PE OCP');
ocp = pe_ocp - ne_ocp;
plot(pe_theta*pe_cmax*litre, ocp, 'k', 'DisplayName', 'Full cell OCP');
xlabel 'Lithiation c (mol/l)';
legend
ylim([0, 5])


sgtitle('after calibration')
